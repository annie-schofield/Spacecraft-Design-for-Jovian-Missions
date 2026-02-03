import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

def calculate_shielding(data_file="spenvis_sao.txt", mission_duration_days=30):
    print(f"--- STARTING SHIELDING OPTIMIZATION ({mission_duration_days} DAYS) ---")
    
    #parsing file, uses re.split to read numbers. 
    #extracts energy (MeV) and flux
    valid_rows = []
    
    try:
        with open(data_file, 'r') as f:
            lines = f.readlines()
            
        print(f"Reading {len(lines)} lines from file...")
        
        for i, line in enumerate(lines):
            #splits by comma, space, or tab (flexible delimiters)
            parts = re.split(r'[,\s]+', line.strip())
            
            #filter out empty strings from the split
            parts = [p for p in parts if p]
            
            #need at least 2 numbers (Energy, ElectronFlux)
            if len(parts) >= 2:
                try:
                    #try converting the first 2 columns to floats
                    energy = float(parts[0])
                    flux_electron = float(parts[1])
                    
                    #if valid, append. We assume Col 1 = Energy, Col 2 = Electrons
                    valid_rows.append([energy, flux_electron])
                except ValueError:
                    #this catches lines like "Energy (MeV)" or "Block #1"
                    continue 

        if not valid_rows:
            print("CRITICAL ERROR: No data found. Please open 'spenvis_sao.txt' and check if it contains numbers.")
            return 0.0
            
        print(f"Success! Parsed {len(valid_rows)} data points.")
        
    except FileNotFoundError:
        print(f"ERROR: Could not find file '{data_file}'.")
        print("Please ensure the text file is in the same folder as this script.")
        return 0.0

    #create DataFrame
    data = pd.DataFrame(valid_rows, columns=["Energy_MeV", "ElectronFlux"])
    
    #convert flux to fluence
    #flux is hits per second, fluence is total hits so multiply by mission length
    seconds_in_mission = mission_duration_days * 24 * 3600
    data['ElectronFluence'] = data['ElectronFlux'] * seconds_in_mission
    
    #filter for physics relevance (> 0.04 MeV), ignore low energy
    target_data = data[data['Energy_MeV'] >= 0.04].copy()
    
    # 3. CALCULATE REQUIRED SHIELDING (Weber's Approximation for Al)
    #calculates required shielding (Weber approximation)
    shielding_required = []
    for e_mev in target_data['Energy_MeV']:
        #empirical range formula (g/cm2)
        if e_mev < 2.5:
            #low energy exponential approximation, E < 2.5 MeV
            if e_mev > 0: # Avoid log(0)
                r = 0.412 * (e_mev**(1.265 - 0.0954 * np.log(e_mev)))
            else:
                r = 0.0
        else:
            #high energy linear approximation, E > 2.5 MeV
            r = 0.530 * e_mev - 0.106
        shielding_required.append(r)
        
    target_data['Shielding_g_cm2'] = shielding_required
    
    #convert to thickness using aluminium density
    target_data['Al_Thickness_mm'] = (target_data['Shielding_g_cm2'] / 2.70) * 10.0 
    
    #we need to stop most of the flux. 
    #find the max energy where the fluence is still dangerous (> 1e9 hits)
    
    hazardous_energies = target_data[target_data['ElectronFluence'] > 1e9]
    
    if not hazardous_energies.empty:
        max_hazard_energy = hazardous_energies['Energy_MeV'].max()
        #find thickness for that specific energy
        design_thickness = hazardous_energies[hazardous_energies['Energy_MeV'] == max_hazard_energy]['Al_Thickness_mm'].values[0]
    else:
        #only if environment is chill
        max_hazard_energy = 0.5
        design_thickness = 2.0 #minimum structural wall
        
    #safety factor
    final_thickness = design_thickness * 1.2

    print(f"Max Hazardous Energy detected: {max_hazard_energy} MeV")
    print(f"Raw Required Thickness: {design_thickness:.2f} mm")
    print(f"Recommended Thickness (with 20% Safety Margin): {final_thickness:.2f} mm")
    
    #plot
    plt.figure(figsize=(10, 5))
    plt.plot(target_data['Energy_MeV'], target_data['Al_Thickness_mm'], 'b-', linewidth=2)
    plt.axhline(y=final_thickness, color='r', linestyle='--', label=f'Design: {final_thickness:.2f}mm')
    plt.xlabel("Electron Energy (MeV)")
    plt.ylabel("Required Shielding (mm Al)")
    plt.title(f"Jovian Shielding Analysis (Mission: {mission_duration_days} days)")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.savefig("shielding_result.png")
    print("Plot saved to 'shielding_result.png'")

if __name__ == "__main__":
    calculate_shielding("spenvis_sao.txt")