# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 09:35:17 2025

@author: rotemb


Code 1 (Screening) 
        is a deterministic run.
        It runs each simulation exactly once to get a quick "Ranking"
        of who inhibits and who doesn't.

Code 2 (Sensitivity Analysis)
        is a stochastic (randomized) run.
        It runs each simulation 10 times with random "noise" to calculate
        P-values and Error Bars, proving that the results are not just luck
        or numerical artifacts.

"""

# === CODE 1 ===
# ==== Imports and Configuration ===

import cobra
import pandas as pd
import numpy as np
import os
import glob
import gc
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
from contextlib import contextmanager
import sys
from scipy import stats
import re



# --- 1. CONFIGURATION ---
# Update these paths to match your folders
MODELS_DIR = r"C:\Users\rotemb\OneDrive\Ph.D\binning_90_5\dFBA_no_bz_rib\final_GSMMs"
MEDIUM_FILE = r"C:\Users\rotemb\OneDrive\Ph.D\binning_90_5\dFBA_no_bz_rib\Apple_Medium.csv"
PENICILLIUM_FILE = "PE_gapfilled.xml"
OUTPUT_FILE = "dFBA_Checkpoint_Results.csv" 

# Simulation Parameters
TOTAL_TIME = 72.0       # Hours
DT = 1.0                # Time step (1 hour is stable)
INITIAL_BIOMASS = 0.01  # gDW/L
DEATH_RATE = 0.005      # 1/h
DEFAULT_VMAX = 10.0     # mmol/gDW/h

# Additives to Apple Medium
#ADDITIVES = {'EX_bz_e': 10.0, 'EX_rib__D_e': 10.0}

# Volatiles to Evaporate
VOLATILES = ['EX_etoh_e', 'EX_ac_e', 'EX_acald_e'] 
EVAPORATION_RATE = 0.01

warnings.filterwarnings('ignore')

# --- 2. DEFINE SMALL HELPER LOGIC ---
# (kept as functions only for calculation logic, not execution flow)

def get_medium_dict(csv_path):
    df = pd.read_csv(csv_path)
    medium = pd.Series(df.Flux.values, index=df.Exchange).to_dict()

    return medium

def safe_apply_bounds(model, medium, biomass, dt):
    """Updates bounds. Ensures LowerBound <= UpperBound to prevent crash."""
    for reaction in model.exchanges:
        rid = reaction.id
        
        # 1. Reset Upper Bound (Open) to prevent conflict
        reaction.upper_bound = 1000.0
        
        # 2. Determine Lower Bound (Uptake)
        if rid in medium:
            available = medium[rid]
            # Avoid divide by zero
            avail_limit = available / (max(biomass, 1e-9) * dt)
            limit = min(DEFAULT_VMAX, avail_limit)
            
            new_lb = -limit
            
            # 3. CRITICAL SAFETY CHECK
            # If floating point math makes new_lb > 1000 (unlikely) or existing bound issues
            if new_lb > reaction.upper_bound:
                reaction.upper_bound = new_lb + 10.0
                
            reaction.lower_bound = new_lb
        else:
            reaction.lower_bound = 0.0

def update_medium(medium, fluxes, biomass, dt):
    """Calculates new concentration: C_new = C_old + Flux * X * dt"""
    for rid, val in fluxes.items():
        if not rid.startswith('EX_'): continue
        change = val * biomass * dt
        medium[rid] = max(0.0, medium.get(rid, 0.0) + change)
    
    # Evaporation
    for v in VOLATILES:
        if v in medium: medium[v] *= (1 - EVAPORATION_RATE * dt)
    return medium

# --- 3. CHECKPOINT INIT ---
# Check what is already done
processed_bacteria = []
if os.path.exists(OUTPUT_FILE):
    try:
        existing_df = pd.read_csv(OUTPUT_FILE)
        processed_bacteria = existing_df['Bacteria'].unique().tolist()
        print(f"Found existing results file. Skipping {len(processed_bacteria)} bacteria already processed.")
    except:
        print("Could not read existing file. Will overwrite/append.")
else:
    # Create new file with header
    with open(OUTPUT_FILE, 'w') as f:
        f.write("Bacteria,Penicillium_Final,Control_Ref,Interaction_Score,Type\n")
    print("Created new results file.")

# --- 4. RUN CONTROL (PENICILLIUM ONLY) ---
# We need this value for the score calculation.
# We run it once.
print("\n--- Step A: Calculating Control (Penicillium Monoculture) ---")
pen_path = os.path.join(MODELS_DIR, PENICILLIUM_FILE)
initial_medium = get_medium_dict(MEDIUM_FILE)

# Variables for Control
p_model = cobra.io.read_sbml_model(pen_path)
p_med = initial_medium.copy()
p_bio = INITIAL_BIOMASS

# Control Loop
for t in np.arange(0, TOTAL_TIME, DT):
    safe_apply_bounds(p_model, p_med, p_bio, DT)
    try:
        sol = p_model.optimize()
        g = sol.objective_value if sol.status == 'optimal' else 0.0
        flx = sol.fluxes if sol.status == 'optimal' else {}
    except:
        g = 0.0; flx = {}
        
    p_bio += (g - DEATH_RATE) * p_bio * DT
    p_med = update_medium(p_med, flx, p_bio, DT)

final_control_biomass = p_bio
print(f"Control Final Biomass: {final_control_biomass:.4f}")

# Cleanup Control Memory
del p_model, p_med, p_bio
gc.collect()

# --- 5. RUN COMPETITION LOOP (FLAT) ---
bact_files = glob.glob(os.path.join(MODELS_DIR, "GSMM_*.xml"))
print(f"\n--- Step B: Starting Competition for {len(bact_files)} Bacteria ---")

for i, b_file in enumerate(bact_files):
    b_name = os.path.basename(b_file).replace('.xml', '')
    
    # CHECKPOINT: Skip if done
    if b_name in processed_bacteria:
        continue
        
    print(f"Processing [{i+1}/{len(bact_files)}]: {b_name} ...", end=" ")
    
    try:
        # A. LOAD FRESH MODELS (Vital for stability)
        model_b = cobra.io.read_sbml_model(b_file)
        model_p = cobra.io.read_sbml_model(pen_path)
        
        # B. INIT STATE
        curr_med = initial_medium.copy()
        bio_b = INITIAL_BIOMASS
        bio_p = INITIAL_BIOMASS
        
        # C. TIME LOOP (72h)
        for t in np.arange(0, TOTAL_TIME, DT):
            # Update Bounds (Shared Medium)
            safe_apply_bounds(model_b, curr_med, bio_b, DT)
            safe_apply_bounds(model_p, curr_med, bio_p, DT)
            
            # Solve Bacteria
            try:
                sol_b = model_b.optimize()
                g_b = sol_b.objective_value if sol_b.status == 'optimal' else 0.0
                f_b = sol_b.fluxes if sol_b.status == 'optimal' else {}
            except: g_b=0.0; f_b={}
            
            # Solve Penicillium
            try:
                sol_p = model_p.optimize()
                g_p = sol_p.objective_value if sol_p.status == 'optimal' else 0.0
                f_p = sol_p.fluxes if sol_p.status == 'optimal' else {}
            except: g_p=0.0; f_p={}
            
            # Update Biomass
            bio_b += (g_b - DEATH_RATE) * bio_b * DT
            bio_p += (g_p - DEATH_RATE) * bio_p * DT
            
            # Clamp (Prevent negative/zero)
            bio_b = max(1e-9, bio_b)
            bio_p = max(1e-9, bio_p)
            
            # Update Medium (Aggregate)
            curr_med = update_medium(curr_med, f_b, bio_b, DT)
            curr_med = update_medium(curr_med, f_p, bio_p, DT)
            
        # D. SAVE RESULT
        score = (bio_p - final_control_biomass) / final_control_biomass
        
        if score < -0.1: i_type = "Inhibition"
        elif score > 0.1: i_type = "Support"
        else: i_type = "Neutral"
        
        # Write single line to CSV
        with open(OUTPUT_FILE, 'a') as f:
            line = f"{b_name},{bio_p:.6f},{final_control_biomass:.6f},{score:.6f},{i_type}\n"
            f.write(line)
            
        print(f"Saved. ({i_type})")
        
        # E. CLEANUP (Crucial!)
        del model_b, model_p, sol_b, sol_p, curr_med
        gc.collect()
        
    except Exception as e:
        print(f"Error with {b_name}: {e}")
        # Try to clean up anyway to save the kernel
        gc.collect()

print("\n--- Script Finished ---")
print(f"Results are in: {OUTPUT_FILE}")



# === Plot the results ===

import matplotlib.colors as mcolors

# --- 1. Load Data ---
csv_file = 'dFBA_Checkpoint_Results.csv'
df = pd.read_csv(csv_file)
    
# --- 2. Rename Bacteria (Standardizing Names) ---
# We map internal IDs to Display Names
rename_map = {
    'GSMM_11_gapfilled_final': 'GSMM_11_gapfilled_(AP)',
    'GSMM_11': 'GSMM_11_(AP)',
    'GSMM_9_gapfilled': 'GSMM_9_gapfilled_(PJ)',
    'GSMM_9': 'GSMM_9_(PJ)',
    'GSMM_23_gapfilled': 'GSMM_23_gapfilled_(PE)',
    'GSMM_23': 'GSMM_23_(PE)',
    'GSMM_8_gapfilled': 'GSMM_8_gapfilled_(KC)',
    'GSMM_8': 'GSMM_8_(KC)',
    'GSMM_19_gapfilled': 'GSMM_19_gapfilled_(CT)',
    'GSMM_19': 'GSMM_19_(CT)',
    'GSMM_15_gapfilled_final': 'GSMM_15_gapfilled_final_(BN)',
    'GSMM_15': 'GSMM_15_(BN)',
    'GSMM_20_gapfilled': 'GSMM_20_gapfilled_(PV)',
    'GSMM_20': 'GSMM_20_(PV)',
    'GSMM_5_gapfilled_final': 'GSMM_5_gapfilled_final_(SY)',
    'GSMM_5': 'GSMM_5_(SY)',
    'GSMM_MR_gapfilled': 'GSMM_MR_gapfilled', # Keep consistency if needed
    'GSMM_MR': 'GSMM_MR'
}

df['Bacteria_Display'] = df['Bacteria'].map(rename_map).fillna(df['Bacteria'])

# --- 3. Define Colors ---
# Helper to lighten color
def lighten_color(color, amount=0.4):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    """
    try:
        c = mcolors.cnames[color]
    except:
        c = color
    c = mcolors.to_rgb(c)
    # Convert to HLS, increase L, convert back
    import colorsys
    h, l, s = colorsys.rgb_to_hls(*c)
    # Make it lighter: increase lightness
    # new_l = l + amount * (1 - l) # amount 0 = no change, 1 = white
    # Alternative simple blending with white
    white = np.array([1, 1, 1])
    c_arr = np.array(c)
    new_c = c_arr + (white - c_arr) * amount
    return mcolors.to_hex(new_c)

# Base Palette (9 distinct colors from Tab10 or similar)
# 1. Blue, 2. Orange, 3. Green, 4. Red, 5. Purple, 6. Brown, 7. Pink, 8. Gray, 9. Olive/Yellow
# Base Palette (9 distinct colors from Tab10)
colors_hex = [
    '#E6194B', # Red (AP)
    '#3CB44B', # Green (MR)
    '#FFE119', # Yellow (PJ)
    '#4363D8', # Blue (PE)
    '#F58231', # Orange (KC)
    '#911EB4', # Purple (CT)
    '#42D4F4', # Cyan (BN)
    '#F032E6', # Magenta (PV)
    '#BFEF45'  # Lime (SY)
    ]


# Map Strains to Base Colors
# We need to link the strain ID (e.g. "11", "MR") to a color index.
strain_order = ['11', 'MR', '9', '23', '8', '19', '15', '20', '5']
strain_color_map = {s: colors_hex[i] for i, s in enumerate(strain_order)}

# Generate the specific color map for the plot
final_color_map = {}

# We iterate through the dataframe to assign colors
for idx, row in df.iterrows():
    disp_name = row['Bacteria_Display']
    
    # Identify the strain key
    strain_key = None
    for s in strain_order:
        # Check if strain identifier is in the name (e.g. GSMM_11...)
        # Special case: MR matches GSMM_MR
        # Numbers match GSMM_11
        if f"_{s}_" in disp_name or f"_{s}" in disp_name or (s == 'MR' and 'MR' in disp_name):
            strain_key = s
            break
    
    if strain_key:
        base = strain_color_map[strain_key]
        if 'gapfilled' in disp_name:
            final_color_map[disp_name] = base # Gapfilled gets Base
        else:
            final_color_map[disp_name] = lighten_color(base, 0.4) # Original gets Lighter
    else:
        final_color_map[disp_name] = '#333333' # Fallback black

# --- 4. Plot ---
df_sorted = df.sort_values('Interaction_Score')
df_sorted['Color'] = df_sorted['Bacteria_Display'].map(final_color_map)

plt.figure(figsize=(12, 14))
bars = plt.barh(df_sorted['Bacteria_Display'], df_sorted['Interaction_Score'], color=df_sorted['Color'])

plt.axvline(x=0, color='black', linestyle='--', linewidth=1.0)
plt.xlabel('Interaction Score\n(Negative = Inhibition, Positive = Support)', fontsize=12)
plt.ylabel('Bacteria Model', fontsize=12)
plt.title('dFBA Interaction Analysis: Bacteria vs Penicillium (72h)', fontsize=14)
plt.grid(axis='x', linestyle='--', alpha=0.5)

# Add Labels
offset = 0.005
for bar in bars:
    width = bar.get_width()
    label_x_pos = width + (offset if width >= 0 else -offset)
    align = 'left' if width >= 0 else 'right'
    plt.text(label_x_pos, bar.get_y() + bar.get_height()/2, f'{width:.2f}', 
             va='center', ha=align, fontsize=9, color='black')

plt.tight_layout()
output_img = 'Interaction_Summary_Plot.tiff'
plt.savefig(output_img, dpi=300)
print(f"Plot saved to {output_img}")


# =============================================================================
# =============================================================================
# === CODE 2 ===


# --- 1. CONFIGURATION ---
OUTPUT_STATS_FILE = "dFBA_Sensitivity_Analysis_All.csv"

# Sensitivity Parameters
N_ITERATIONS = 10       # Runs per pair (10 runs * 33 models = 330 total simulations)
NOISE_LEVEL = 0.10      # 10% variability in Vmax
BASE_VMAX = 10.0        # Mean uptake rate
TOTAL_TIME = 72.0
DT = 1.0
INITIAL_BIOMASS = 0.01
DEATH_RATE = 0.005

# Additives
VOLATILES = ['EX_etoh_e', 'EX_ac_e', 'EX_acald_e'] 
EVAPORATION_RATE = 0.01

warnings.filterwarnings('ignore')

# --- 2. HELPER FUNCTIONS ---

@contextmanager
def suppress_output():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout; old_stderr = sys.stderr
        sys.stdout = devnull; sys.stderr = devnull
        try: yield
        finally: sys.stdout = old_stdout; sys.stderr = old_stderr

def get_medium_dict(csv_path):
    df = pd.read_csv(csv_path)
    medium = pd.Series(df.Flux.values, index=df.Exchange).to_dict()

    return medium

def safe_apply_bounds(model, medium, biomass, dt, current_vmax):
    """Applies bounds using a VARIABLE Vmax (noisy)."""
    for reaction in model.exchanges:
        rid = reaction.id
        reaction.upper_bound = 1000.0
        
        if rid in medium:
            available = medium[rid]
            # Limit based on availability
            avail_limit = available / (max(biomass, 1e-9) * dt)
            
            # Limit based on physiology (Includes random noise)
            limit = min(current_vmax, avail_limit)
            
            new_lb = -limit
            if new_lb > reaction.upper_bound: reaction.upper_bound = new_lb + 10.0
            reaction.lower_bound = new_lb
        else:
            reaction.lower_bound = 0.0

def update_medium(medium, fluxes, biomass, dt):
    for rid, val in fluxes.items():
        if not rid.startswith('EX_'): continue
        change = val * biomass * dt
        medium[rid] = max(0.0, medium.get(rid, 0.0) + change)
    for v in VOLATILES:
        if v in medium: medium[v] *= (1 - EVAPORATION_RATE * dt)
    return medium

# --- 3. IDENTIFY ALL TARGETS ---
print("--- Step 1: Finding All Bacterial Models ---")
bact_files = glob.glob(os.path.join(MODELS_DIR, "GSMM_*.xml"))
all_bacteria = [os.path.basename(f).replace('.xml', '') for f in bact_files]

print(f"Found {len(all_bacteria)} models to analyze.")

# --- 4. RUN SENSITIVITY ANALYSIS ---

# Store all results here: {'Control': [val1...], 'Bacteria_A': [val1...]}
biomass_distributions = {}

# A. Control (Penicillium Only) - Run N times
print(f"\n--- Step 2: Running Control Analysis (N={N_ITERATIONS}) ---")
control_biomasses = []
pen_path = os.path.join(MODELS_DIR, PENICILLIUM_FILE)
base_medium = get_medium_dict(MEDIUM_FILE)

for n in range(N_ITERATIONS):
    print(f"   Iter {n+1}/{N_ITERATIONS}...", end="\r")
    
    # Generate Noisy Parameter
    noisy_vmax = np.random.normal(loc=BASE_VMAX, scale=BASE_VMAX * NOISE_LEVEL)
    noisy_vmax = max(1.0, noisy_vmax) 
    
    try:
        with suppress_output():
            p_model = cobra.io.read_sbml_model(pen_path)
        
        p_med = base_medium.copy()
        p_bio = INITIAL_BIOMASS
        
        for t in np.arange(0, TOTAL_TIME, DT):
            safe_apply_bounds(p_model, p_med, p_bio, DT, noisy_vmax)
            try:
                sol = p_model.optimize()
                g = sol.objective_value if sol.status == 'optimal' else 0.0
                flx = sol.fluxes if sol.status == 'optimal' else {}
            except: g=0; flx={}
            
            p_bio += (g - DEATH_RATE) * p_bio * DT
            p_med = update_medium(p_med, flx, p_bio, DT)
        
        control_biomasses.append(p_bio)
        del p_model, p_med
        gc.collect()
        
    except Exception as e:
        print(f"Error in control iter {n}: {e}")

biomass_distributions['Control'] = control_biomasses
print(f"   Done. Mean Control: {np.mean(control_biomasses):.4f} +/- {np.std(control_biomasses):.4f}")


# B. Co-cultures (All Models)
print(f"\n--- Step 3: Running Co-culture Analysis on {len(all_bacteria)} models ---")

stats_results = []

for idx, b_name in enumerate(all_bacteria):
    print(f"[{idx+1}/{len(all_bacteria)}] Analyzing {b_name} ...")
    
    b_file = os.path.join(MODELS_DIR, f"{b_name}.xml")
    
    coculture_biomasses = []
    
    for n in range(N_ITERATIONS):
        print(f"   Iter {n+1}/{N_ITERATIONS}...", end="\r")
        
        noisy_vmax = np.random.normal(loc=BASE_VMAX, scale=BASE_VMAX * NOISE_LEVEL)
        noisy_vmax = max(1.0, noisy_vmax)
        
        try:
            with suppress_output():
                model_b = cobra.io.read_sbml_model(b_file)
                model_p = cobra.io.read_sbml_model(pen_path)
            
            curr_med = base_medium.copy()
            bio_b = INITIAL_BIOMASS
            bio_p = INITIAL_BIOMASS
            
            for t in np.arange(0, TOTAL_TIME, DT):
                safe_apply_bounds(model_b, curr_med, bio_b, DT, noisy_vmax)
                safe_apply_bounds(model_p, curr_med, bio_p, DT, noisy_vmax)
                
                try:
                    sol_b = model_b.optimize()
                    g_b = sol_b.objective_value if sol_b.status == 'optimal' else 0.0
                    f_b = sol_b.fluxes if sol_b.status == 'optimal' else {}
                except: g_b=0; f_b={}
                
                try:
                    sol_p = model_p.optimize()
                    g_p = sol_p.objective_value if sol_p.status == 'optimal' else 0.0
                    f_p = sol_p.fluxes if sol_p.status == 'optimal' else {}
                except: g_p=0; f_p={}
                
                bio_b += (g_b - DEATH_RATE) * bio_b * DT
                bio_p += (g_p - DEATH_RATE) * bio_p * DT
                
                curr_med = update_medium(curr_med, f_b, bio_b, DT)
                curr_med = update_medium(curr_med, f_p, bio_p, DT)
            
            coculture_biomasses.append(bio_p)
            del model_b, model_p, curr_med
            gc.collect()
            
        except Exception as e:
            print(f"Error in {b_name} iter {n}: {e}")

    # Store Data
    biomass_distributions[b_name] = coculture_biomasses
    
    # Statistical Test (Welch's T-test)
    t_stat, p_val = stats.ttest_ind(control_biomasses, coculture_biomasses, equal_var=False)
    
    # Effect Size
    mean_control = np.mean(control_biomasses)
    mean_coculture = np.mean(coculture_biomasses)
    
    # Check for empty results (crashes)
    if mean_control > 0:
        interaction_score = (mean_coculture - mean_control) / mean_control
    else:
        interaction_score = 0.0
    
    stats_results.append({
        'Bacteria': b_name,
        'Mean_Penicillium_Biomass': mean_coculture,
        'Std_Penicillium_Biomass': np.std(coculture_biomasses),
        'Interaction_Score': interaction_score,
        'P_Value': p_val,
        'Significant': 'Yes' if p_val < 0.05 else 'No'
    })
    
    print(f"   Done. Score: {interaction_score:.2f}, P-value: {p_val:.4e}")

# --- 5. SAVE ---

df_stats = pd.DataFrame(stats_results)
df_stats.to_csv(OUTPUT_STATS_FILE, index=False)
print(f"\nStats saved to {OUTPUT_STATS_FILE}")
print(df_stats)

# --- 6. PLOT ---

# Load from file if available, else use a placeholder dict as this is a plotting request.
# But I must use the uploaded file if I can access it.
csv_file = 'dFBA_Sensitivity_Analysis_All.csv'
df_stats = pd.read_csv(csv_file)


df_stats['Bacteria_Display'] = df_stats['Bacteria'].map(rename_map).fillna(df_stats['Bacteria'])


# Sort by Mean Biomass
df_stats.sort_values('Mean_Penicillium_Biomass', inplace=True)

# --- 4. Plot Horizontal ---
plt.figure(figsize=(10, 12))

names = df_stats['Bacteria_Display'].tolist()
means = df_stats['Mean_Penicillium_Biomass'].tolist()
stds = df_stats['Std_Penicillium_Biomass'].tolist()
bar_colors = [final_color_map.get(n, '#808080') for n in names]

# Control Line
# Assuming control_mean = 29.10 based on previous files
control_mean = 29.10 
# std? Assuming small or 0 if unknown. Let's use 0.5 for visual shading if exact not known
control_std = 0.5 

plt.axvline(x=control_mean, color='black', linestyle='--', label='Control Mean')
plt.axvspan(control_mean - control_std, control_mean + control_std, color='black', alpha=0.1, label='Control StdDev')

y_pos = np.arange(len(names))
plt.barh(y_pos, means, xerr=stds, align='center', alpha=0.9, color=bar_colors, ecolor='black', capsize=3)

plt.yticks(y_pos, names, fontsize=10)
plt.xlabel('Penicillium Biomass (gDW/L)', fontsize=12)
plt.title(f'Sensitivity Analysis: Bacterial Impact on Fungal Growth', fontsize=14)
plt.legend(loc='lower right')
plt.tight_layout()

output_file = 'Sensitivity_Analysis_All_Models_Horizontal.tiff'
plt.savefig(output_file, dpi=300)
print(f"Plot saved as {output_file}")





# --- 1. Load Data ---
csv_file = 'dFBA_Sensitivity_Analysis_All.csv'
try:
    df_stats = pd.read_csv(csv_file)
except:
    print("Error: CSV not found. Using dummy data for demonstration.")
    # Placeholder for structure if file is missing
    df_stats = pd.DataFrame(columns=['Bacteria', 'Mean_Penicillium_Biomass', 'Std_Penicillium_Biomass'])

# --- 2. Rename Bacteria (Scientific Names) ---
rename_map = {
    'GSMM_11_gapfilled_final': 'A. pittii',
    'GSMM_11': 'A. pittii',
    'GSMM_9_gapfilled': 'P. juntendii',
    'GSMM_9': 'P. juntendii',
    'GSMM_23_gapfilled': 'P. endophytica',
    'GSMM_23': 'P. endophytica',
    'GSMM_8_gapfilled': 'K. cowanii',
    'GSMM_8': 'K. cowanii',
    'GSMM_19_gapfilled': 'C. terrigena',
    'GSMM_19': 'C. terrigena',
    'GSMM_15_gapfilled_final': 'B. nasdae',
    'GSMM_15': 'B. nasdae',
    'GSMM_20_gapfilled': 'P. vagans',
    'GSMM_20': 'P. vagans',
    'GSMM_5_gapfilled_final': 'S. yunnanensis',
    'GSMM_5': 'S. yunnanensis',
    'GSMM_MR_gapfilled': 'M. radiotolerans',
    'GSMM_MR': 'M. radiotolerans'
}

df_stats['Bacteria_Display'] = df_stats['Bacteria'].map(rename_map).fillna(df_stats['Bacteria'])

# --- 3. Split Data ---
df_stats['Is_Gapfilled'] = df_stats['Bacteria'].str.contains('gapfilled', case=False)
df_gapfilled = df_stats[df_stats['Is_Gapfilled']].copy()
df_original = df_stats[~df_stats['Is_Gapfilled']].copy()

# Sort by Mean Biomass
df_gapfilled.sort_values('Mean_Penicillium_Biomass', inplace=True)
df_original.sort_values('Mean_Penicillium_Biomass', inplace=True)

# --- 4. Define Colors ---
palette_hex = [
    '#E6194B', # Red (A. pittii)
    '#3CB44B', # Green (M. radiotolerans)
    '#FFE119', # Yellow (P. juntendii)
    '#4363D8', # Blue (P. endophytica)
    '#F58231', # Orange (K. cowanii)
    '#911EB4', # Purple (C. terrigena)
    '#42D4F4', # Cyan (B. nasdae)
    '#F032E6', # Magenta (P. vagans)
    '#BFEF45'  # Lime (S. yunnanensis)
]
# Map pairs to colors using the unique display name content
pair_codes = ['A. pittii', 'M. radiotolerans', 'P. juntendii', 'P. endophytica', 
              'K. cowanii', 'C. terrigena', 'B. nasdae', 'P. vagans', 'S. yunnanensis']
pair_color_dict = {code: color for code, color in zip(pair_codes, palette_hex)}

def lighten_color(color, amount=0.6):
    try: c = mcolors.cnames[color]
    except: c = color
    c = mcolors.to_rgb(c)
    white = np.array([1, 1, 1])
    c_arr = np.array(c)
    new_c = c_arr + (white - c_arr) * amount
    return mcolors.to_hex(new_c)

def get_color(name, is_gapfilled):
    # Match specific names from the pair_codes list
    matched_code = None
    for code in pair_codes:
        if code in name:
            matched_code = code
            break
            
    if matched_code:
        base = pair_color_dict[matched_code]
        return base if is_gapfilled else lighten_color(base, 0.6)
    return '#808080'

# Assign colors
df_gapfilled['Color'] = df_gapfilled['Bacteria_Display'].apply(lambda x: get_color(x, True))
df_original['Color'] = df_original['Bacteria_Display'].apply(lambda x: get_color(x, False))

# --- 5. Plotting Function ---
def plot_sensitivity(df, title, filename):
    plt.figure(figsize=(10, 8))
    
    names = df['Bacteria_Display'].tolist()
    means = df['Mean_Penicillium_Biomass'].tolist()
    stds = df['Std_Penicillium_Biomass'].tolist()
    colors = df['Color'].tolist()
    
    # Control Line
    control_mean = 29.10 
    control_std = 0.5 

    plt.axvline(x=control_mean, color='black', linestyle='--', label='Control Mean')
    plt.axvspan(control_mean - control_std, control_mean + control_std, color='black', alpha=0.1, label='Control StdDev')

    y_pos = np.arange(len(names))
    plt.barh(y_pos, means, xerr=stds, align='center', alpha=0.9, color=colors, ecolor='black', capsize=3)

    plt.yticks(y_pos, names, fontsize=10, style='italic') # Italic for scientific names
    plt.xlabel('Penicillium Biomass (gDW/L)', fontsize=12)
    plt.title(title, fontsize=14)
    
    # Extend X Axis
    plt.xlim(0, 33) 
    
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Plot saved as {filename}")

# --- 6. Generate Plots ---
if not df_gapfilled.empty:
    plot_sensitivity(df_gapfilled, 'Sensitivity Analysis: Gapfilled Models', 'Sensitivity_Gapfilled.tiff')
if not df_original.empty:
    plot_sensitivity(df_original, 'Sensitivity Analysis: Original Models', 'Sensitivity_Original.tiff')







# === Corelation seperate for gapfilled and original models ===



# --- 1. CONFIGURATION ---
SIM_FILE = 'dFBA_Checkpoint_Results.csv'
SENSITIVITY_FILE = 'dFBA_Sensitivity_Analysis_All.csv'
EXP_FILES = {
    'Diffusable': 'diffusable_ex_results.csv',
    'Fruit': 'fruits_ex_results.csv'
}

# --- 2. LOAD DATA ---
try: df_sim = pd.read_csv(SIM_FILE)
except: df_sim = pd.DataFrame()

try: df_sens = pd.read_csv(SENSITIVITY_FILE)
except: df_sens = pd.DataFrame()

def load_exp(filename, label):
    if not os.path.exists(filename): return pd.DataFrame()
    try:
        df = pd.read_csv(filename)
        df.columns = [c.strip() for c in df.columns]
        gsmm_col = next((c for c in df.columns if 'GSMM' in c), None)
        val_col = [c for c in df.columns if 'Day' in c][-1]
        
        if not gsmm_col or not val_col: return pd.DataFrame()
        
        expanded = []
        for _, row in df.iterrows():
            raw = str(row[gsmm_col])
            val = row[val_col]
            for s in raw.split(','):
                s_clean = s.strip().replace('"', '').replace("'", "").replace('gapfillied', 'gapfilled')
                if s_clean == 'MR_gapfilled': s_clean = 'GSMM_MR_gapfilled'
                expanded.append({'Bacteria': s_clean, 'Value': val})
        
        df_agg = pd.DataFrame(expanded).groupby('Bacteria')['Value'].mean().reset_index()
        return df_agg.rename(columns={'Value': f'Exp_{label}'})
    except: return pd.DataFrame()

df_diff = load_exp(EXP_FILES['Diffusable'], 'Diffusable')
df_fruit = load_exp(EXP_FILES['Fruit'], 'Fruit')

# --- 3. MERGE & PROCESS ---
if not df_sim.empty:
    df_corr = pd.merge(df_sim, df_diff, on='Bacteria', how='inner')
    if not df_fruit.empty: df_corr = pd.merge(df_corr, df_fruit, on='Bacteria', how='left')
    
    # Normalize Z-Scores
    def zscore(col):
        if col.std() == 0: return col - col.mean() # Zero variance fallback
        return (col - col.mean()) / col.std()

    df_corr['Sim_Z'] = zscore(df_corr['Interaction_Score'])
    df_corr['Exp_Diff_Z'] = zscore(df_corr['Exp_Diffusable'])
    df_corr['Exp_Fruit_Z'] = zscore(df_corr['Exp_Fruit'])
else:
    df_corr = pd.DataFrame() # Fallback

# --- 4. COLORS & RENAMING (Updated Names) ---
rename_map = {
    'GSMM_11_gapfilled_final': 'A. pittii', 'GSMM_11': 'A. pittii',
    'GSMM_9_gapfilled': 'P. juntendii', 'GSMM_9': 'P. juntendii',
    'GSMM_23_gapfilled': 'P. endophytica', 'GSMM_23': 'P. endophytica',
    'GSMM_8_gapfilled': 'K. cowanii', 'GSMM_8': 'K. cowanii',
    'GSMM_19_gapfilled': 'C. terrigena', 'GSMM_19': 'C. terrigena',
    'GSMM_15_gapfilled_final': 'B. nasdae', 'GSMM_15': 'B. nasdae',
    'GSMM_20_gapfilled': 'P. vagans', 'GSMM_20': 'P. vagans',
    'GSMM_5_gapfilled_final': 'S. yunnanensis', 'GSMM_5': 'S. yunnanensis',
    'GSMM_MR_gapfilled': 'M. radiotolerans', 'GSMM_MR': 'M. radiotolerans'
}

if not df_corr.empty:
    df_corr['Bacteria_Display'] = df_corr['Bacteria'].map(rename_map).fillna(df_corr['Bacteria'])
if not df_sens.empty:
    df_sens['Bacteria_Display'] = df_sens['Bacteria'].map(rename_map).fillna(df_sens['Bacteria'])

# Colors
palette_hex = ['#E6194B', '#3CB44B', '#FFE119', '#4363D8', '#F58231', '#911EB4', '#42D4F4', '#F032E6', '#BFEF45']
pair_codes = ['A. pittii', 'M. radiotolerans', 'P. juntendii', 'P. endophytica', 'K. cowanii', 'C. terrigena', 'B. nasdae', 'P. vagans', 'S. yunnanensis']
pair_dict = {c: h for c, h in zip(pair_codes, palette_hex)}

def get_color(name, is_gapfilled):
    # Match updated names
    code = next((p for p in pair_codes if p in name), None)
    if not code: return '#808080'
    base = pair_dict[code]
    if is_gapfilled: return base
    try: c = mcolors.to_rgb(base)
    except: return base
    white = np.array([1,1,1])
    return mcolors.to_hex(np.array(c) + (white - np.array(c))*0.6)

if not df_corr.empty:
    df_corr['Is_Gapfilled'] = df_corr['Bacteria'].str.contains('gapfilled', case=False)
    df_corr['Color'] = df_corr.apply(lambda x: get_color(x['Bacteria_Display'], x['Is_Gapfilled']), axis=1)

if not df_sens.empty:
    df_sens['Is_Gapfilled'] = df_sens['Bacteria'].str.contains('gapfilled', case=False)
    df_sens['Color'] = df_sens.apply(lambda x: get_color(x['Bacteria_Display'], x['Is_Gapfilled']), axis=1)

# --- 5. PLOTTING ---
def generate_figure(is_gapfilled, title_main, filename):
    subset_corr = pd.DataFrame()
    if not df_corr.empty:
        subset_corr = df_corr[df_corr['Is_Gapfilled'] == is_gapfilled].copy()
    
    subset_sens = pd.DataFrame()
    if not df_sens.empty:
        subset_sens = df_sens[df_sens['Is_Gapfilled'] == is_gapfilled].copy()
        subset_sens.sort_values('Mean_Penicillium_Biomass', inplace=True)

    fig, axes = plt.subplots(1, 3, figsize=(20, 10))
    fig.suptitle(title_main, fontsize=16)

    # Common Aesthetics
    for ax in axes:
        ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7, zorder=0)
        ax.set_axisbelow(True) # Grid behind bars/points

    # Plot 1: Diffusable
    if not subset_corr.empty:
        if subset_corr['Sim_Z'].std() > 0:
            r, p = stats.spearmanr(subset_corr['Sim_Z'], subset_corr['Exp_Diff_Z']) 
            sns.regplot(x='Sim_Z', y='Exp_Diff_Z', data=subset_corr, ax=axes[0], scatter=False, color='blue', ci=95)
            stats_text = f"Spearman rho={r:.2f}, p={p:.3f}"
        else:
            stats_text = "No Variance (N/A)"
            
        axes[0].plot([-2, 1], [-2, 1], 'k--', alpha=0.3, label='Perfect Match')
        axes[0].scatter(subset_corr['Sim_Z'], subset_corr['Exp_Diff_Z'], c=subset_corr['Color'], s=150, edgecolors='black', zorder=10)
        axes[0].set_title(f"Diffusable Experiment Correlation\n{stats_text}", fontsize=12)
        axes[0].set_xlabel('Simulation Z-Score')
        axes[0].set_ylabel('Experiment Z-Score')
        axes[0].autoscale(enable=True, axis='both', tight=True)
        axes[0].margins(x=0.05, y=0.05) 

    # Plot 2: Fruit
    if not subset_corr.empty and not subset_corr['Exp_Fruit_Z'].isnull().all():
        sub_fruit = subset_corr.dropna(subset=['Exp_Fruit_Z'])
        if not sub_fruit.empty:
            if sub_fruit['Sim_Z'].std() > 0:
                r, p = stats.spearmanr(sub_fruit['Sim_Z'], sub_fruit['Exp_Fruit_Z'])
                sns.regplot(x='Sim_Z', y='Exp_Fruit_Z', data=sub_fruit, ax=axes[1], scatter=False, color='green', ci=95)
                stats_text = f"Spearman rho={r:.2f}, p={p:.3f}"
            else:
                stats_text = "No Variance (N/A)"

            axes[1].plot([-2, 1], [-2, 1], 'k--', alpha=0.3)
            axes[1].scatter(sub_fruit['Sim_Z'], sub_fruit['Exp_Fruit_Z'], c=sub_fruit['Color'], s=150, edgecolors='black', zorder=10)
            axes[1].set_title(f"Fruit Experiment Correlation\n{stats_text}", fontsize=12)
            axes[1].set_xlabel('Simulation Z-Score')
            axes[1].set_ylabel('Experiment Z-Score')
            axes[1].autoscale(enable=True, axis='both', tight=True)
            axes[1].margins(x=0.05, y=0.05)

    # Plot 3: Sensitivity Bar Plot
    if not subset_sens.empty:
        names = subset_sens['Bacteria_Display'].tolist()
        means = subset_sens['Mean_Penicillium_Biomass'].tolist()
        stds = subset_sens['Std_Penicillium_Biomass'].tolist()
        colors = subset_sens['Color'].tolist()
        y_pos = np.arange(len(names))
        
        control_mean = 29.10; control_std = 0.5
        axes[2].axvline(x=control_mean, color='black', linestyle='--', label='Control Mean')
        axes[2].axvspan(control_mean - control_std, control_mean + control_std, color='black', alpha=0.1)
        axes[2].barh(y_pos, means, xerr=stds, align='center', alpha=0.9, color=colors, ecolor='black', capsize=3)
        axes[2].set_yticks(y_pos)
        axes[2].set_yticklabels(names, fontsize=9, style='italic') # Italic names
        axes[2].set_xlabel('Penicillium Biomass (gDW/L)')
        axes[2].set_title('Sensitivity Analysis (Predicted Biomass)', fontsize=12)
        axes[2].legend(loc='lower right') # Legend placement
        
        # Extended X-Axis
        axes[2].set_xlim(0, 33) 
        axes[2].margins(y=0.02) 

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

# EXECUTE
generate_figure(True, "Gapfilled Models: Correlations & Sensitivity", "Combined_Analysis_Gapfilled_Scientific.tiff")
generate_figure(False, "Original Models: Correlations & Sensitivity", "Combined_Analysis_Original_Scientific.tiff")








'''
barplots for the lab experiments
'''


from statsmodels.stats.multicomp import pairwise_tukeyhsd
import networkx as nx

# Define codes to scientific names mapping
code_to_name = {
    'AP': 'A. pittii',
    'MR': 'M. radiotolerans',
    'PJ': 'P. juntendii',
    'PE': 'P. endophytica',
    'KC': 'K. cowanii',
    'CT': 'C. terrigena',
    'BN': 'B. nasdae',
    'PV': 'P. vagans',
    'SY': 'S. yunnanensis',
    'CONTROL': 'Control',
    'Control': 'Control',
    'Water': 'Control',
    'C': 'Control'
}

# Define colors (same hex codes, mapped to new names)
colors_hex = [
    '#E6194B', # Red (A. pittii)
    '#3CB44B', # Green (M. radiotolerans)
    '#FFE119', # Yellow (P. juntendii)
    '#4363D8', # Blue (P. endophytica)
    '#F58231', # Orange (K. cowanii)
    '#911EB4', # Purple (C. terrigena)
    '#42D4F4', # Cyan (B. nasdae)
    '#F032E6', # Magenta (P. vagans)
    '#BFEF45', # Lime (S. yunnanensis)
    '#808080'  # Gray (Control)
]

ordered_names = [
    'A. pittii', 'M. radiotolerans', 'P. juntendii', 'P. endophytica',
    'K. cowanii', 'C. terrigena', 'B. nasdae', 'P. vagans', 'S. yunnanensis', 'Control'
]

color_map = dict(zip(ordered_names, colors_hex))

# Load data
try:
    fruits_df = pd.read_csv('fruits_ex_results.csv')
    diffusable_df = pd.read_csv('diffusable_ex_results.csv')
except FileNotFoundError:
    print("Error: Files not found.")
    fruits_df = pd.DataFrame()
    diffusable_df = pd.DataFrame()

def clean_and_map(df):
    if df.empty: return df
    df.columns = [c.strip() for c in df.columns]
    # Find Isolates column
    iso_cols = [c for c in df.columns if 'Isolate' in c]
    if iso_cols:
        iso_col = iso_cols[0]
        df['Isolates_Raw'] = df[iso_col].astype(str).str.strip()
        df['Isolates_Display'] = df['Isolates_Raw'].map(code_to_name).fillna(df['Isolates_Raw'])
        # Identify last day column
        day_cols = [c for c in df.columns if 'Day' in c]
        if day_cols:
            df['Value'] = df[day_cols[-1]]
        return df
    return df

fruits_df = clean_and_map(fruits_df)
diffusable_df = clean_and_map(diffusable_df)

def get_cld(df, value_col, group_col):
    if df.empty or df[group_col].nunique() < 2:
        return {g: '' for g in df[group_col].unique()}
        
    try:
        tukey = pairwise_tukeyhsd(endog=df[value_col], groups=df[group_col], alpha=0.05)
        tukey_results = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])
        
        groups = df[group_col].unique()
        G = nx.Graph()
        G.add_nodes_from(groups)
        
        for _, row in tukey_results.iterrows():
            if not row['reject']:
                G.add_edge(row['group1'], row['group2'])
        
        cliques = list(nx.find_cliques(G))
        means = df.groupby(group_col)[value_col].mean()
        
        def get_clique_mean(clique):
            return means[list(clique)].mean()
        
        cliques.sort(key=get_clique_mean, reverse=True)
        
        import string
        letters = string.ascii_lowercase
        group_letters = {g: '' for g in groups}
        
        for i, clique in enumerate(cliques):
            letter = letters[i] if i < len(letters) else str(i)
            for group in clique:
                if letter not in group_letters[group]:
                    group_letters[group] += letter
        
        for g in group_letters:
            group_letters[g] = ''.join(sorted(group_letters[g]))
            
        return group_letters
    except Exception as e:
        print(f"Stats Error: {e}")
        return {g: '' for g in df[group_col].unique()}

def process_dataset(df):
    if df.empty or 'Value' not in df.columns: return pd.DataFrame()
    df_clean = df.dropna(subset=['Value', 'Isolates_Display'])
    
    summary = df_clean.groupby('Isolates_Display')['Value'].agg(['mean', 'std']).reset_index()
    summary = summary.sort_values('mean', ascending=False)
    
    cld = get_cld(df_clean, 'Value', 'Isolates_Display')
    summary['cld'] = summary['Isolates_Display'].map(cld)
    
    return summary

fruits_summary = process_dataset(fruits_df)
diffusable_summary = process_dataset(diffusable_df)

# Plotting
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

def plot_bar(ax, summary, title):
    if summary.empty:
        ax.text(0.5, 0.5, 'No Data', ha='center')
        ax.set_title(title)
        return

    bar_colors = [color_map.get(x, '#333333') for x in summary['Isolates_Display']]
    
    ax.bar(summary['Isolates_Display'], summary['mean'], yerr=summary['std'], 
           capsize=5, color=bar_colors, edgecolor='black')
    
    ax.set_title(title, fontsize=14)
    ax.set_ylabel('Day 7 Measurement', fontsize=12)
    # Set italic labels for scientific names
    ax.set_xticklabels(summary['Isolates_Display'], rotation=45, ha='right', style='italic')
    
    # Fix Y-Axis and Labels
    y_max = summary['mean'].max() + (summary['std'].max() if pd.notna(summary['std'].max()) else 0)
    ax.set_ylim(0, y_max * 1.15)
        
    for i, (idx, row) in enumerate(summary.iterrows()):
        h = row['mean'] + (row['std'] if pd.notna(row['std']) else 0)
        ax.text(i, h + (y_max*0.02), row['cld'], 
                ha='center', va='bottom', fontsize=10, fontweight='bold')

plot_bar(axes[0], fruits_summary, 'Fruits Experiment')
plot_bar(axes[1], diffusable_summary, 'Plates')

plt.tight_layout()
plt.savefig('experiment_barplots.tiff', dpi=300)
print("Plot saved as experiment_barplots.tiff")




# Define colors
colors_hex = [
    '#E6194B', # Red (AP)
    '#3CB44B', # Green (MR)
    '#FFE119', # Yellow (PJ)
    '#4363D8', # Blue (PE)
    '#F58231', # Orange (KC)
    '#911EB4', # Purple (CT)
    '#42D4F4', # Cyan (BN)
    '#F032E6', # Magenta (PV)
    '#BFEF45'  # Lime (SY)
]
isolate_keys = ['AP', 'MR', 'PJ', 'PE', 'KC', 'CT', 'BN', 'PV', 'SY']
color_map = dict(zip(isolate_keys, colors_hex))
color_map['CONTROL'] = '#808080' # Gray for Control

# Load data
volatiles_df = pd.read_csv('volatiles_ex_results.csv')

def get_cld(df, value_col, group_col):
    # Get mean for sorting
    means = df.groupby(group_col)[value_col].mean()
    
    # Run Tukey HSD
    tukey = pairwise_tukeyhsd(endog=df[value_col], groups=df[group_col], alpha=0.05)
    tukey_results = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])
    
    # Create graph for non-significant differences
    groups = df[group_col].unique()
    G = nx.Graph()
    G.add_nodes_from(groups)
    
    # Add edges for non-significant pairs
    for _, row in tukey_results.iterrows():
        if not row['reject']:
            G.add_edge(row['group1'], row['group2'])
            
    # Find cliques
    cliques = list(nx.find_cliques(G))
    
    # Sort cliques by their average mean value
    def get_clique_mean(clique):
        return means[list(clique)].mean()
    
    cliques.sort(key=get_clique_mean, reverse=True)
    
    # Assign letters
    import string
    letters = string.ascii_lowercase
    group_letters = {g: '' for g in groups}
    
    for i, clique in enumerate(cliques):
        letter = letters[i] if i < len(letters) else str(i)
        for group in clique:
            if letter not in group_letters[group]:
                group_letters[group] += letter
            
    # Sort letters
    for g in group_letters:
        group_letters[g] = ''.join(sorted(group_letters[g]))
        
    return group_letters

def process_dataset(df, name):
    # Aggregating
    summary = df.groupby('Isolates')['Day7'].agg(['mean', 'std']).reset_index()
    summary = summary.sort_values('mean', ascending=False)
    
    # Get CLD
    cld = get_cld(df, 'Day7', 'Isolates')
    summary['cld'] = summary['Isolates'].map(cld)
    
    return summary

volatiles_summary = process_dataset(volatiles_df, 'Volatiles')

# Plotting
fig, ax = plt.subplots(figsize=(8, 6))

def plot_bar(ax, summary, title):
    bars = ax.bar(summary['Isolates'], summary['mean'], yerr=summary['std'], capsize=5, 
                  color=[color_map.get(x, '#333333') for x in summary['Isolates']], edgecolor='black')
    
    ax.set_title(title, fontsize=14)
    ax.set_ylabel('Day 7 Measurement', fontsize=12)
    ax.set_xticklabels(summary['Isolates'], rotation=45, ha='right')
    
    # Add letters
    y_max = summary['mean'].max() + summary['std'].max()
    ylim = ax.get_ylim()
    # Ensure there is enough space for letters
    ax.set_ylim(0, max(ylim[1], y_max * 1.15)) 
    
    ylim_updated = ax.get_ylim()
    
    for i, (idx, row) in enumerate(summary.iterrows()):
        ax.text(i, row['mean'] + row['std'] + (ylim_updated[1]*0.02), row['cld'], 
                ha='center', va='bottom', fontsize=10, fontweight='bold')

plot_bar(ax, volatiles_summary, 'Volatiles Experiment')

plt.tight_layout()
plt.savefig('volatiles_barplot.png', dpi=300)

print(volatiles_summary)









'''
This script executes a dynamic Flux Balance Analysis (dFBA)
to simulate the metabolic competition between Penicillium fungus and various
bacterial strains over a 72-hour period.
It iteratively simulates the co-culture growth, tracking the specific consumption
of key nutrients (glucose, fructose, sucrose, and ammonium) and dynamically updating
the shared medium. Finally, it classifies the interaction as inhibition or support
based on the fungal biomass and saves the detailed flux results to a CSV file.
'''



# --- 1. CONFIGURATION ---
MODELS_DIR = r"C:\Users\rotemb\OneDrive\Ph.D\binning_90_5\dFBA_no_bz_rib\final_GSMMs"
MEDIUM_FILE = r"C:\Users\rotemb\OneDrive\Ph.D\binning_90_5\dFBA_no_bz_rib\Apple_Medium.csv"
PENICILLIUM_FILE = "PE_gapfilled.xml"
OUTPUT_FILE = "dFBA_flux_Results_new.csv" 

# Simulation Parameters
TOTAL_TIME = 72.0       # Hours
DT = 1.0                # Time step
INITIAL_BIOMASS = 0.01  # gDW/L
DEATH_RATE = 0.005      # 1/h
DEFAULT_VMAX = 10.0     # mmol/gDW/h

# Additives & Volatiles
#ADDITIVES = {'EX_bz_e': 10.0, 'EX_rib__D_e': 10.0}
VOLATILES = ['EX_etoh_e', 'EX_ac_e', 'EX_acald_e'] 
EVAPORATION_RATE = 0.01

warnings.filterwarnings('ignore')

# --- 2. HELPER FUNCTIONS ---
def get_medium_dict(csv_path):
    df = pd.read_csv(csv_path)
    medium = pd.Series(df.Flux.values, index=df.Exchange).to_dict()

    return medium

def safe_apply_bounds(model, medium, biomass, dt):
    for reaction in model.exchanges:
        rid = reaction.id
        reaction.upper_bound = 1000.0 # Reset
        
        if rid in medium:
            available = medium[rid]
            avail_limit = available / (max(biomass, 1e-9) * dt)
            limit = min(DEFAULT_VMAX, avail_limit)
            new_lb = -limit
            
            if new_lb > reaction.upper_bound:
                reaction.upper_bound = new_lb + 10.0
            reaction.lower_bound = new_lb
        else:
            reaction.lower_bound = 0.0

def update_medium(medium, fluxes, biomass, dt):
    for rid, val in fluxes.items():
        if not rid.startswith('EX_'): continue
        change = val * biomass * dt
        medium[rid] = max(0.0, medium.get(rid, 0.0) + change)
    
    for v in VOLATILES:
        if v in medium: medium[v] *= (1 - EVAPORATION_RATE * dt)
    return medium

# --- 3. CHECKPOINT INIT ---
processed_bacteria = []
if os.path.exists(OUTPUT_FILE):
    try:
        existing_df = pd.read_csv(OUTPUT_FILE)
        # Check if the file has the new columns. If not, we should probably stop or warn.
        if 'Glc_Consumed' not in existing_df.columns:
            print("WARNING: Existing file does not have Nutrient columns. Please delete it and restart.")
            sys.exit(1)
            
        processed_bacteria = existing_df['Bacteria'].unique().tolist()
        print(f"Found existing results. Skipping {len(processed_bacteria)} bacteria.")
    except:
        print("Could not read existing file. Will overwrite.")
else:
    # --- FIX IS HERE: Added the new headers to match the data ---
    with open(OUTPUT_FILE, 'w') as f:
        f.write("Bacteria,Penicillium_Final,Bacteria_Final,Control_Ref,Interaction_Score,Type,Glc_Consumed,Fru_Consumed,Sucr_Consumed,Nh4_Consumed\n")    
        print("Created new results file.")

# --- 4. RUN CONTROL (PENICILLIUM ONLY) ---
print("\n--- Step A: Calculating Control (Penicillium Monoculture) ---")
pen_path = os.path.join(MODELS_DIR, PENICILLIUM_FILE)
initial_medium = get_medium_dict(MEDIUM_FILE)

p_model = cobra.io.read_sbml_model(pen_path)
p_med = initial_medium.copy()
p_bio = INITIAL_BIOMASS

for t in np.arange(0, TOTAL_TIME, DT):
    safe_apply_bounds(p_model, p_med, p_bio, DT)
    try:
        sol = p_model.optimize()
        g = sol.objective_value if sol.status == 'optimal' else 0.0
        flx = sol.fluxes if sol.status == 'optimal' else {}
    except: g = 0.0; flx = {}
        
    p_bio += (g - DEATH_RATE) * p_bio * DT
    p_med = update_medium(p_med, flx, p_bio, DT)

final_control_biomass = p_bio
print(f"Control Final Biomass: {final_control_biomass:.4f}")
del p_model, p_med, p_bio
gc.collect()

# --- 5. RUN COMPETITION LOOP ---
bact_files = glob.glob(os.path.join(MODELS_DIR, "GSMM_*.xml"))
print(f"\n--- Step B: Starting Competition for {len(bact_files)} Bacteria ---")

for i, b_file in enumerate(bact_files):
    b_name = os.path.basename(b_file).replace('.xml', '')
    
    if b_name in processed_bacteria: continue
        
    print(f"Processing [{i+1}/{len(bact_files)}]: {b_name} ...", end=" ")
    
    try:
        model_b = cobra.io.read_sbml_model(b_file)
        model_p = cobra.io.read_sbml_model(pen_path)
        
        curr_med = initial_medium.copy()
        bio_b = INITIAL_BIOMASS
        bio_p = INITIAL_BIOMASS
        
        # --- 1. Reset Counters for this bacteria ---
        total_glc = 0.0
        total_fru = 0.0
        total_sucr = 0.0
        total_nh4 = 0.0
        
        for t in np.arange(0, TOTAL_TIME, DT):
            safe_apply_bounds(model_b, curr_med, bio_b, DT)
            safe_apply_bounds(model_p, curr_med, bio_p, DT)
            
            # Solve Bacteria
            try:
                sol_b = model_b.optimize()
                g_b = sol_b.objective_value if sol_b.status == 'optimal' else 0.0
                f_b = sol_b.fluxes.to_dict() if sol_b.status == 'optimal' else {} 
            except: g_b=0.0; f_b={}
            
            # Solve Penicillium
            try:
                sol_p = model_p.optimize()
                g_p = sol_p.objective_value if sol_p.status == 'optimal' else 0.0
                f_p = sol_p.fluxes if sol_p.status == 'optimal' else {}
            except: g_p=0.0; f_p={}
            
            # --- 2. Accumulate Consumption (Flux * Biomass * DT) ---
            # Note: Flux is typically negative for uptake, so we add the raw negative value.
            
            # Glucose
            total_glc += f_b.get('EX_glc__D_e', f_b.get('EX_glc_e', 0.0)) * bio_b * DT
            
            # Fructose (Try standard IDs)
            total_fru += f_b.get('EX_fru_e', f_b.get('EX_fru__D_e', 0.0)) * bio_b * DT
            
            # Sucrose
            total_sucr += f_b.get('EX_sucr_e', 0.0) * bio_b * DT
            
            # Ammonium
            total_nh4 += f_b.get('EX_nh4_e', 0.0) * bio_b * DT
            
            # Update State
            bio_b += (g_b - DEATH_RATE) * bio_b * DT
            bio_p += (g_p - DEATH_RATE) * bio_p * DT
            bio_b = max(1e-9, bio_b)
            bio_p = max(1e-9, bio_p)
            
            curr_med = update_medium(curr_med, f_b, bio_b, DT)
            curr_med = update_medium(curr_med, f_p, bio_p, DT)
            
        # Score
        score = (bio_p - final_control_biomass) / final_control_biomass
        if score < -0.1: i_type = "Inhibition"
        elif score > 0.1: i_type = "Support"
        else: i_type = "Neutral"
        
        # --- 3. Save with New Columns ---
        with open(OUTPUT_FILE, 'a') as f:
            line = f"{b_name},{bio_p:.6f},{bio_b:.6f},{final_control_biomass:.6f},{score:.6f},{i_type},{total_glc:.4f},{total_fru:.4f},{total_sucr:.4f},{total_nh4:.4f}\n"
            f.write(line)
            
        print(f"Saved. ({i_type})")
        
        del model_b, model_p, sol_b, sol_p, curr_med
        gc.collect()
        
    except Exception as e:
        print(f"Error with {b_name}: {e}")
        gc.collect()




'''
This script performs a targeted secretion analysis on a specific list of bacterial
strains (e.g., Strain 15, 19, MR) to identify potential volatile inhibitors.
It runs a 72-hour dFBA simulation for each target model, accumulating the total amount
of every metabolite secreted (positive flux) into the medium.
The results are aggregated into a pivot table and saved as a CSV file,
enabling direct comparison of toxin production profiles across the different
inhibitory and non-inhibitory strains.
'''


# --- CONFIGURATION ---
MODELS_DIR = r"C:\Users\rotemb\OneDrive\Ph.D\binning_90_5\dFBA_no_bz_rib\final_GSMMs"
MEDIUM_FILE = r"C:\Users\rotemb\OneDrive\Ph.D\binning_90_5\dFBA_no_bz_rib\Apple_Medium.csv"
OUTPUT_FILE = "Full_Secretion_Profile.csv"

# LIST OF TARGETS (Based on your Volatiles Experiment)
# We map the Isolate to the Model Name pattern
TARGET_PATTERNS = [
    'GSMM_15', 'GSMM_20', 'GSMM_MR', 'GSMM_19', 
    'GSMM_9', 'GSMM_8', 'GSMM_23', 'GSMM_11', 'GSMM_5'
]

# Simulation Params
TOTAL_TIME = 72.0
DT = 1.0
INITIAL_BIOMASS = 0.01
DEFAULT_VMAX = 10.0
warnings.filterwarnings('ignore')

# --- HELPERS ---
def get_medium_dict(csv_path):
    df = pd.read_csv(csv_path)
    return pd.Series(df.Flux.values, index=df.Exchange).to_dict()

def safe_apply_bounds(model, medium, biomass, dt):
    for reaction in model.exchanges:
        rid = reaction.id
        reaction.upper_bound = 1000.0
        if rid in medium:
            limit = min(DEFAULT_VMAX, medium[rid] / (max(biomass, 1e-9) * dt))
            reaction.lower_bound = -limit
        else:
            reaction.lower_bound = 0.0

def update_medium(medium, fluxes, biomass, dt):
    for rid, val in fluxes.items():
        if not rid.startswith('EX_'): continue
        change = val * biomass * dt
        medium[rid] = max(0.0, medium.get(rid, 0.0) + change)
    return medium

# --- MAIN LOOP ---
print(f"Scanning secretions for {len(TARGET_PATTERNS)} strains...")
all_secretions = []
initial_medium = get_medium_dict(MEDIUM_FILE)

# Find files
files = glob.glob(os.path.join(MODELS_DIR, "*.xml"))
target_files = []
for f in files:
    fname = os.path.basename(f)
    # Only process if it matches our list and is a 'gapfilled' model (active)
    if 'gapfilled' in fname and any(t in fname for t in TARGET_PATTERNS):
        target_files.append(f)

print(f"Found {len(target_files)} model files.")

for b_file in target_files:
    name = os.path.basename(b_file).replace('.xml', '')
    print(f"Running {name}...", end=" ")
    
    try:
        model = cobra.io.read_sbml_model(b_file)
        medium = initial_medium.copy()
        biomass = INITIAL_BIOMASS
        secretion_totals = {} 
        
        for t in np.arange(0, TOTAL_TIME, DT):
            safe_apply_bounds(model, medium, biomass, DT)
            try:
                sol = model.optimize()
                fluxes = sol.fluxes.to_dict() if sol.status == 'optimal' else {}
                g = sol.objective_value if sol.status == 'optimal' else 0.0
            except: fluxes={}; g=0.0
            
            # Record Positive Fluxes (Secretion)
            for rxn_id, flux_val in fluxes.items():
                if rxn_id.startswith('EX_') and flux_val > 0.001: # Filter noise
                    amount = flux_val * biomass * DT
                    secretion_totals[rxn_id] = secretion_totals.get(rxn_id, 0.0) + amount
            
            biomass += (g - 0.005) * biomass * DT
            medium = update_medium(medium, fluxes, biomass, DT)
            
        for met, amount in secretion_totals.items():
            all_secretions.append({'Bacteria': name, 'Metabolite': met, 'Total_Secretion_mmol': amount})
        
        print("Done.")
        del model
        gc.collect()
        
    except Exception as e:
        print(f"Error: {e}")

# Save
df_sec = pd.DataFrame(all_secretions)
if not df_sec.empty:
    # Pivot: Index=Metabolite, Columns=Bacteria
    df_pivot = df_sec.pivot_table(index='Metabolite', columns='Bacteria', values='Total_Secretion_mmol').fillna(0)
    df_pivot.to_csv(OUTPUT_FILE)
    print(f"\nSaved full secretion profile to {OUTPUT_FILE}")
else:
    print("No secretions found.")







'''
This script generates a targeted heatmap of bacterial volatile secretions
(excluding the pathogen) and comparative barplots for Diffusable, Fruit, and Volatiles
inhibition assays. It applies a Tukey HSD post-hoc test to statistically group
the experimental results, annotating the sorted barplots with significance
letters to distinguish effective inhibitors from non-inhibitory strains.
'''




# --- 1. LOAD DATA ---
try:
    df_sec = pd.read_csv('Full_Secretion_Profile.csv', index_col=0) # Index is Metabolite
except Exception as e:
    print(f"Error loading file: {e}")
    df_sec = pd.DataFrame()

# --- 2. CONFIGURATION ---
# New scientific names mapping
rename_map_strains = {
    'GSMM_11_gapfilled_final': 'A. pittii', 'GSMM_11': 'A. pittii',
    'GSMM_9_gapfilled': 'P. juntendii', 'GSMM_9': 'P. juntendii',
    'GSMM_23_gapfilled': 'P. endophytica', 'GSMM_23': 'P. endophytica',
    'GSMM_8_gapfilled': 'K. cowanii', 'GSMM_8': 'K. cowanii',
    'GSMM_19_gapfilled': 'C. terrigena', 'GSMM_19': 'C. terrigena',
    'GSMM_15_gapfilled_final': 'B. nasdae', 'GSMM_15': 'B. nasdae',
    'GSMM_20_gapfilled': 'P. vagans', 'GSMM_20': 'P. vagans',
    'GSMM_5_gapfilled_final': 'S. yunnanensis', 'GSMM_5': 'S. yunnanensis',
    'GSMM_MR_gapfilled': 'M. radiotolerans', 'GSMM_MR': 'M. radiotolerans'
}

# Target Metabolites (Original + New)
target_mets = [
    'EX_etoh_e', 'EX_acald_e', 'EX_ac_e', 'EX_lac__D_e', 
    'EX_lac__L_e', 'EX_succ_e', 'EX_h2_e', 'EX_for_e',
    'EX_nh4_e', 'EX_no2_e', 'EX_leu__L_e', 'EX_slnt_e'
]

# Metabolite Name Mapping
met_rename_map = {
    'EX_etoh_e': 'Ethanol',
    'EX_acald_e': 'Acetaldehyde',
    'EX_ac_e': 'Acetate',
    'EX_lac__D_e': 'D-Lactate',
    'EX_lac__L_e': 'L-Lactate',
    'EX_succ_e': 'Succinate',
    'EX_h2_e': 'Hydrogen',
    'EX_for_e': 'Formate',
    'EX_nh4_e': 'Ammonium (NH4)',
    'EX_no2_e': 'Nitrite (NO2)',
    'EX_leu__L_e': 'Leucine',
    'EX_slnt_e': 'Selenite'
}

# --- 3. HEATMAP GENERATION ---
if not df_sec.empty:
    df_heatmap = df_sec.copy()
    
    # 1. Rename Columns (Strains)
    df_heatmap.rename(columns=rename_map_strains, inplace=True)
    
    # 2. Exclude PE (Pathogen)
    # Filter out columns containing "PE" (unless it's P. endophytica which is mapped to 'P. endophytica')
    # BE CAREFUL: "PE" is in "P. endophytica". 
    # Original logic: 'PE' not in c. 
    # If the column name is NOW 'P. endophytica', 'PE' is NOT in it (case sensitive).
    # But we should filter BEFORE renaming to be safe, OR check for "Penicillium".
    # Assuming the input dataframe has original names like "GSMM_PE...", we can filter first.
    
    # Safer approach: Filter columns that are NOT in our target list of bacteria
    valid_bacteria = list(rename_map_strains.values())
    # Keep only columns that match one of our known bacterial names
    cols_to_keep = [c for c in df_heatmap.columns if c in valid_bacteria]
    df_heatmap = df_heatmap[cols_to_keep]
    
    # 3. Filter Rows (Metabolites)
    df_heatmap = df_heatmap.loc[df_heatmap.index.intersection(target_mets)]
    
    # 4. Rename Rows (Metabolites)
    df_heatmap.rename(index=met_rename_map, inplace=True)
    
    # Sort for aesthetics
    df_heatmap.sort_index(axis=0, inplace=True) # Sort metabolites
    df_heatmap.sort_index(axis=1, inplace=True) # Sort strains
    
    # Plot
    plt.figure(figsize=(12, 8))
    sns.heatmap(df_heatmap, cmap="Reds", linewidths=.5, annot=True, fmt=".1f")
    plt.title("Bacterial Secretion Profile (mmol/72h)")
    
    # Rotate x-labels for better readability
    plt.xticks(rotation=45, ha='right', style='italic') # Italic for species names
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    plt.savefig("Volatile_Secretion_Heatmap_NoPE_Final.png", dpi=300)
    plt.show()
    print("Heatmap saved.")
else:
    print("Secretion data empty. Please ensure 'Full_Secretion_Profile.csv' is in the directory.")




