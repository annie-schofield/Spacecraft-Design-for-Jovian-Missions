import pandas as pd
from datetime import datetime, timedelta
import numpy as np
from tudatpy import constants
from tudatpy.interface import spice
from tudatpy import numerical_simulation
from tudatpy.numerical_simulation import environment_setup
from tudatpy.numerical_simulation import propagation_setup
from tudatpy.astro import element_conversion

def run_ganymede_simulation():
    """
    Propagates orbit of spacecraft around Ganymede. 
    Includes solar radiation pressure pertubations.
    Considers 30 day period. 
    

    Returns
    -------
    state_history : contains state_history with epochs and positions.

    """
    #load spice kernals
    spice.load_standard_kernels()
    
    #create environment
    bodies_to_create = ["Sun", "Jupiter", "Ganymede", "Europa", "Callisto"]

    global_frame_origin = "Ganymede"  
    global_frame_orientation = "ECLIPJ2000"
    
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, 
        global_frame_origin, 
        global_frame_orientation
    )

    #define spacecraft
    body_settings.add_empty_settings("Spacecraft")
    body_settings.get("Spacecraft").constant_mass = 1000.0
    
    #radiation pressure settings
    reference_area_radiation = 5.0
    radiation_pressure_coefficient = 1.2
    occulting_bodies = {"Sun": ["Jupiter"]} #source : blockers
    
    radiation_pressure_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
        reference_area_radiation, 
        radiation_pressure_coefficient, 
        occulting_bodies
    )
    
    body_settings.get("Spacecraft").radiation_pressure_target_settings = radiation_pressure_settings

    #create system of bodies
    bodies = environment_setup.create_system_of_bodies(body_settings)

    #define acceleration settings
    accelerations_settings_spacecraft = dict(
        Ganymede=[
            propagation_setup.acceleration.point_mass_gravity()
        ],
        Jupiter=[
            propagation_setup.acceleration.point_mass_gravity()
        ],
        Sun=[
            propagation_setup.acceleration.point_mass_gravity(),
            propagation_setup.acceleration.radiation_pressure()
        ],
        Europa=[
            propagation_setup.acceleration.point_mass_gravity()
        ]
    )

    acceleration_settings = {"Spacecraft": accelerations_settings_spacecraft}
    
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, 
        acceleration_settings, 
        ["Spacecraft"], # bodies_to_propagate
        ["Ganymede"]    # central_bodies
    )
    #define initial orbit
    mu_ganymede = bodies.get("Ganymede").gravitational_parameter
    
    #define orbit parameters
    radius_ganymede = 2634.0 * 1000.0  
    altitude = 500.0 * 1000.0  
    
    sma = radius_ganymede + altitude
    ecc = 0.0
    inc = np.deg2rad(90.0) # polar orbit
    arg_per = 0.0
    raan = 0.0
    true_anom = 0.0

    initial_state = element_conversion.keplerian_to_cartesian_elementwise(
        semi_major_axis=sma,
        eccentricity=ecc,
        inclination=inc,
        argument_of_periapsis=arg_per,
        longitude_of_ascending_node=raan,
        true_anomaly=true_anom,
        gravitational_parameter=mu_ganymede
    )

    #propagation settings
    simulation_start_epoch = 24.0 * constants.JULIAN_YEAR 
    simulation_duration = 2.0 * constants.JULIAN_DAY
    simulation_end_epoch = simulation_start_epoch + simulation_duration

    termination_settings = propagation_setup.propagator.time_termination(
        simulation_end_epoch
    )

    #using RK4 integrator
    integrator_settings = propagation_setup.integrator.runge_kutta_4(
        initial_time_step=10.0
    )

    propagator_settings = propagation_setup.propagator.translational(
        central_bodies=["Ganymede"],
        acceleration_models=acceleration_models,
        bodies_to_integrate=["Spacecraft"],
        initial_states=initial_state,
        initial_time=simulation_start_epoch,
        integrator_settings=integrator_settings,
        termination_settings=termination_settings
    )

    #run simulation
    dynamics_simulator = numerical_simulation.create_dynamics_simulator(
        bodies, 
        propagator_settings
    )
    
    state_history = dynamics_simulator.state_history

    print(f"Simulation finished. Generated {len(state_history)} data points.")
    return state_history

def export_results_for_spenvis(results_dict, output_filename="spenvis_upload_optimized.txt"):
    """
    Generates txt file that can be used in SPENVIS -> coordinate generator  -> spacecraft trajectories  -> upload trajectory file
    This is 30 days of data.
    Reduces to 60 second intervals so it should be within SPENVIS input limits.

    Parameters
    ----------
    results_dict : uses results from run_ganymede_simulation()
    output_filename : the default is "spenvis_upload_optimized.txt".

    Returns
    -------
    None.

    """
    print(f"Generating optimized SPENVIS file: {output_filename} ...")
    
    #J2000 Epoch in julian days
    jd_j2000 = 2451545.0
    
    #gets all keys (times) and values (states)
    all_times = list(results_dict.keys())
    all_states = list(results_dict.values())
    
    #take every 6th point (10s step -> 60s step)
    #this reduces ~17,000 points to ~2,800 points, to meet SPENVIS limits
    times_downsampled = all_times[::6]
    states_downsampled = all_states[::6]
    
    print(f"Downsampling: Reduced {len(all_times)} points to {len(times_downsampled)} points.")
    
    data_lines = []
    
    for i, epoch_seconds in enumerate(times_downsampled):
        state_sc_wrt_gan = states_downsampled[i]
        
        #convert to Julian date
        current_jd = jd_j2000 + (epoch_seconds / 86400.0)
        
        #ganymede position relative to jupiter (J2000 Frame)
        state_gan_wrt_jup = spice.get_body_cartesian_state_at_epoch(
            target_body_name="Ganymede",
            observer_body_name="Jupiter",
            reference_frame_name="J2000", 
            aberration_corrections="None",
            ephemeris_time=epoch_seconds
        )
        
        #convert metres to km
        x_km = (state_sc_wrt_gan[0] + state_gan_wrt_jup[0]) / 1000.0
        y_km = (state_sc_wrt_gan[1] + state_gan_wrt_jup[1]) / 1000.0
        z_km = (state_sc_wrt_gan[2] + state_gan_wrt_jup[2]) / 1000.0
        
        #format: JD, X, Y, Z
        line = f"{current_jd:.9f}, {x_km:.6f}, {y_km:.6f}, {z_km:.6f}"
        data_lines.append(line)
        
    #writes file
    with open(output_filename, "w") as f:
        # header section
        f.write("******************************************************************\n")
        f.write("* Optimized Ganymede Orbiter Trajectory (60s step)\n")
        f.write("******************************************************************\n")
        f.write("Title: Ganymede Orbiter Analysis\n")
        f.write("Planet: Jupiter\n")
        f.write("Coordinates: PEI\n")
        f.write("Columns: JDCT, X, Y, Z\n")
        f.write("Format: CSV\n")
        f.write("******************************************************************\n")
        f.write("$$BEGIN\n")
        for line in data_lines:
            f.write(line + "\n")
        f.write("$$END\n")

    print(f"Success! Saved to '{output_filename}'. Upload this one to SPENVIS.")

def export_to_ccsds_oem(results_dict, output_filename="ganymede_orbiter_safe.oem"):
    """
    Generates oem file that can be used in SPENVIS -> JOREM -> trajectory upload.
    This is 30 days of data.
    Reduces to 300 second intervals so it should be within SPENVIS input limits.

    Parameters
    ----------
    results_dict : uses results from run_ganymede_simulation()
    output_filename : the default is "ganymede_orbiter_safe.oem".

    Returns
    -------
    None.

    """
    print(f"Generating optimized CCSDS OEM file: {output_filename} ...")
    
    # J2000 epoch reference
    base_date = datetime(2000, 1, 1, 12, 0, 0)
    
    all_times = list(results_dict.keys())
    all_states = list(results_dict.values())

    #original step = 10s, we want 300sto reduce file size even more.
    step_size = 30
    
    times_downsampled = all_times[::step_size]
    states_downsampled = all_states[::step_size]
    
    print(f"Downsampling: Reduced {len(all_times)} points to {len(times_downsampled)} points (Safe for JOREM).")
    
    with open(output_filename, "w") as f:
        # HEADER
        f.write("CCSDS_OEM_VERS = 2.0\n")
        f.write(f"CREATION_DATE  = {datetime.now().strftime('%Y-%m-%dT%H:%M:%S')}\n")
        f.write("ORIGINATOR     = TUDAT_PYTHON\n")
        f.write("\n")
        
        # METADATA
        f.write("META_START\n")
        f.write("OBJECT_NAME          = GANYMEDE_ORBITER\n")
        f.write("OBJECT_ID            = 999\n") 
        f.write("CENTER_NAME          = JUPITER\n")
        f.write("REF_FRAME            = EME2000\n")
        f.write("TIME_SYSTEM          = UTC\n")
        
        t_start = base_date + timedelta(seconds=times_downsampled[0])
        t_stop = base_date + timedelta(seconds=times_downsampled[-1])
        f.write(f"START_TIME           = {t_start.strftime('%Y-%m-%dT%H:%M:%S')}\n")
        f.write(f"STOP_TIME            = {t_stop.strftime('%Y-%m-%dT%H:%M:%S')}\n")
        f.write("META_STOP\n")
        f.write("\n")
        
        # DATA LOOP
        for i, epoch_seconds in enumerate(times_downsampled):
            state_sc_wrt_gan = states_downsampled[i]
            
            # Get Ganymede position relative to Jupiter
            state_gan_wrt_jup = spice.get_body_cartesian_state_at_epoch(
                target_body_name="Ganymede",
                observer_body_name="Jupiter",
                reference_frame_name="J2000", 
                aberration_corrections="None",
                ephemeris_time=epoch_seconds
            )
            
            # Position (km)
            x = (state_sc_wrt_gan[0] + state_gan_wrt_jup[0]) / 1000.0
            y = (state_sc_wrt_gan[1] + state_gan_wrt_jup[1]) / 1000.0
            z = (state_sc_wrt_gan[2] + state_gan_wrt_jup[2]) / 1000.0
            
            # Velocity (km/s)
            vx = (state_sc_wrt_gan[3] + state_gan_wrt_jup[3]) / 1000.0
            vy = (state_sc_wrt_gan[4] + state_gan_wrt_jup[4]) / 1000.0
            vz = (state_sc_wrt_gan[5] + state_gan_wrt_jup[5]) / 1000.0
            
            t_current = base_date + timedelta(seconds=epoch_seconds)
            t_str = t_current.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]
            
            f.write(f"{t_str} {x:.5f} {y:.5f} {z:.5f} {vx:.6f} {vy:.6f} {vz:.6f}\n")
            
    print(f"Success! Saved '{output_filename}'. Upload this smaller file to JOREM.")

if __name__ == "__main__":
    results = run_ganymede_simulation()
    export_results_for_spenvis(results)
    export_to_ccsds_oem(results)