# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 08:47:58 2025

@author: rotemb
"""


import cobra
import os
import sys
import pandas as pd
import json


print("Starting Growth Test Workflow...")

# --- 1. SET FILE PATHS ---
MODELS_DIR = r"C:\Users\rotemb\OneDrive\Ph.D\binning_90_5\CarvMe_new_version\models"
CURATION_DIR = r"C:\Users\rotemb\OneDrive\Ph.D\binning_90_5\CarvMe_new_version\curation\curated_models"
UNIVERSAL_REACTIONS_DIR = r"C:\Users\rotemb\OneDrive\Ph.D\binning_90_5\CarvMe_new_version\curation"
MISSING_REACTIONS_DIR = r"C:\Users\rotemb\OneDrive\Ph.D\binning_90_5\CarvMe_new_version\curation\missing_reactions"


original_model_name = 'GSMM_8.xml'
missing_reactions_file = 'bin.8_(GSMM_8)_missing_reactions.tsv'
universal_reactions_file = 'universal_bigg_reactions.json'
output_model_name = 'GSMM_8_gapfilled.xml'

original_model_path = os.path.join(MODELS_DIR, original_model_name)
missing_reactions_path = os.path.join(MISSING_REACTIONS_DIR, missing_reactions_file)
universal_reactions_path = os.path.join(UNIVERSAL_REACTIONS_DIR, universal_reactions_file)
output_model_path = os.path.join(CURATION_DIR, output_model_name)


# --- 2. LOAD THE ORIGINAL MODEL ---
try:
    model = cobra.io.read_sbml_model(original_model_path)
    print(f"Successfully loaded model '{model.id}'")
except Exception as e:
    print(f"FATAL ERROR: Could not read the model file: {e}")
    sys.exit()

# --- 3. GET OBJECTIVE REACTION ---
try:
    growth_reaction = model.reactions.get_by_id('Growth')
    print(f"Found growth reaction: '{growth_reaction.id}'")
except KeyError:
    print(f"Warning: Could not find 'Growth'. Trying 'R_Growth'.")
    try:
        growth_reaction = model.reactions.get_by_id('R_Growth')
        print(f"Found growth reaction: '{growth_reaction.id}'")
    except KeyError:
        print(f"FATAL ERROR: Could not find 'Growth' or 'R_Growth' reaction.")
        sys.exit()

# Found growth reaction: 'Growth'


# --- 4. DEFINE FULL MINIMAL MEDIUM (Base + Trace) ---
base_lab_medium = {
    "EX_nh4_e": 10.0, "EX_so4_e": 10.0, "EX_k_e": 10.0, "EX_pi_e": 10.0,
    "EX_na1_e": 10.0, "EX_cit_e": 2.0, "EX_mg2_e": 10.0, "EX_h2o_e": 10.0,
    "EX_o2_e": 10.0, "EX_sucr_e": 10.0
}
trace_elements_exchanges = [
    'EX_fe2_e', 'EX_fe3_e', 'EX_zn2_e', 'EX_mn2_e', 'EX_cu2_e',
    'EX_cobalt2_e', 'EX_cl_e', 'EX_ca2_e', 'EX_mobd_e', 'EX_ni2_e',
    'EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_h_e'
]


full_minimal_medium = base_lab_medium.copy()
for ex_id in trace_elements_exchanges:
    if ex_id not in full_minimal_medium:
        full_minimal_medium[ex_id] = 10.0

# --- 5. HELPER FUNCTION TO SET MEDIUM ---
def set_minimal_medium(model, medium_dict):
    print(f"\nSetting full minimal medium ({len(medium_dict)} components)...")
    for ex in model.exchanges:
        ex.lower_bound = 0.0
    for ex_id, uptake in medium_dict.items():
        try:
            exchange_reaction = model.reactions.get_by_id(ex_id)
            exchange_reaction.lower_bound = -abs(uptake)
        except KeyError:
            pass # We know 'EX_tungs_e' is missing, that's fine.

set_minimal_medium(model, full_minimal_medium)


# --- 6. FIND BLOCKED PRECURSORS ---
print("Testing production of all biomass precursors...")
blocked_precursors = []
components = growth_reaction.reactants

print(f"Found {len(components)} biomass components in '{growth_reaction.id}'.")


# Testing production of all biomass precursors...
# Found 53 biomass components in 'Growth'.


for met in components:
    with model:
        # Create a temporary demand reaction for this metabolite
        # A 'demand' reaction allows the metabolite to leave the system
        demand_reaction = model.add_boundary(met, type='demand')
        
        # Set the objective to maximize the production of this one metabolite
        model.objective = demand_reaction
        
        solution = model.optimize()
        
        # If flux is zero (or very close), we can't make it
        if solution.objective_value < 0.000001:
            blocked_precursors.append(met)

# --- 7. PRINT THE LIST" ---
print("--- BLOCKED PRECURSOR LIST ---")
print("The model CANNOT produce the following biomass components on this medium:")
if blocked_precursors:
    for met in blocked_precursors:
        print(f"  - {met.id} ({met.name})")
else:
    print("  None!")
    print("  Growth is 0 but all precursors can be made,")
    print("  A different problem (e.g., a blocked 'atp_c' reaction,")
    print("  'Growth' reaction itself has bounds=[0, 0]).")



# ================================================================================
# --- BLOCKED PRECURSOR TO-DO LIST ---
# The model CANNOT produce the following biomass precursors on this medium:
    #   - 10fthf_c (10-Formyltetrahydrofolate)
    #   - mlthf_c (5,10-Methylenetetrahydrofolate)
    #   - mql8_c (Menaquinol 8)
    #   - phe__L_c (L-Phenylalanine)
    #   - thf_c (5,6,7,8-Tetrahydrofolate)
    #   - thmpp_c (Thiamine diphosphate)
    #   - trp__L_c (L-Tryptophan)
    #   - tyr__L_c (L-Tyrosine)
# ================================================================================



# --- 8. ITERATIVE PATHWAY TRACER ---
print("--- TRACING ALL BLOCKED PATHWAYS (1 by 1) ---")

# The list of 8 precursors you provided
blocked_list = [
    '10fthf_c', 'mlthf_c', 'mql8_c', 'phe__L_c', 
    'thf_c', 'thmpp_c', 'trp__L_c', 'tyr__L_c'
]

# This will store our findings for the DataFrame
all_trace_results = []

for start_met_id in blocked_list:
    print(f"\nTracing pathway for: {start_met_id}")
    met_stack = [start_met_id]
    checked_metabolites = set()
    found_root_cause = False

    while met_stack:
        current_met_id = met_stack.pop(0)
        
        if current_met_id in checked_metabolites:
            continue
        checked_metabolites.add(current_met_id)

        try:
            met = model.metabolites.get_by_id(current_met_id)
        except KeyError:
            print(f"  [WARNING] Precursor '{current_met_id}' not found in model. Skipping.")
            continue

        # Test if this metabolite can be produced
        with model:
            demand = model.add_boundary(met, type='demand')
            model.objective = demand
            flux = model.optimize().objective_value

        if flux > 1e-6:
            # This metabolite CAN be produced. It is not the problem.
            continue
        
        # This metabolite CANNOT be produced. This is a blocked point.
        
        producing_reactions = [rxn for rxn in met.reactions if rxn.get_coefficient(met.id) > 0]
        
        if not producing_reactions:
            print(f"  [ROOT CAUSE FOUND] Cannot produce '{current_met_id}'")
            result = {
                "Blocked_Biomass_Precursor": start_met_id,
                "Root_Cause_Gap_Metabolite": current_met_id,
                "Notes": "No producing reaction in model"
            }
            all_trace_results.append(result)
            found_root_cause = True
            break # Stop tracing this pathway, move to the next in the for-loop
        else:
            # Not a root cause, add its precursors to the stack
            for rxn in producing_reactions:
                for precursor in rxn.reactants:
                    if precursor.id not in checked_metabolites:
                        met_stack.append(precursor.id)
    
    if not found_root_cause:
        print(f"  [NO ROOT CAUSE FOUND] Could not find a root gap for {start_met_id}.")
        result = {
            "Blocked_Biomass_Precursor": start_met_id,
            "Root_Cause_Gap_Metabolite": "Unknown",
            "Notes": "Pathway is blocked, but no metabolite with 0 producers was found."
        }
        all_trace_results.append(result)


# Tracing pathway for: 10fthf_c
#   [ROOT CAUSE FOUND] Cannot produce 'skm_e'
# 
# Tracing pathway for: mlthf_c
#   [ROOT CAUSE FOUND] Cannot produce 'skm_e'
# 
# Tracing pathway for: mql8_c
#   [ROOT CAUSE FOUND] Cannot produce 'dmso_e'
# 
# Tracing pathway for: phe__L_c
#   [ROOT CAUSE FOUND] Cannot produce 'phe__L_e'
# 
# Tracing pathway for: thf_c
#   [ROOT CAUSE FOUND] Cannot produce 'skm_e'
# 
# Tracing pathway for: thmpp_c
#   [ROOT CAUSE FOUND] Cannot produce '2ahethmpp_c'
# 
# Tracing pathway for: trp__L_c
#   [ROOT CAUSE FOUND] Cannot produce 'trp__L_e'
# 
# Tracing pathway for: tyr__L_c
#   [ROOT CAUSE FOUND] Cannot produce 'tyr__L_e'



# --- 6. CREATE AND PRINT THE DATAFRAME ---
print("--- GSMM_8 GAP-FILLING TO-DO LIST ---")

# Convert the results into a DataFrame
results_df = pd.DataFrame(all_trace_results)

# Clean up and display the DataFrame
if not results_df.empty:
    # Remove duplicate root causes (e.g., chorismate)
    unique_root_causes_df = results_df.drop_duplicates(subset='Root_Cause_Gap_Metabolite')
    print("The following root-cause gaps were found:")
    print(unique_root_causes_df.to_string(index=False))
    
    print("\n--- Full Trace Report ---")
    print(results_df.to_string(index=False))
else:
    print("No blocked pathways were found. This is unexpected.")


# --- GSMM_8 GAP-FILLING TO-DO LIST ---
# The following root-cause gaps were found:
# Blocked_Biomass_Precursor Root_Cause_Gap_Metabolite                          Notes
#                  10fthf_c                     skm_e No producing reaction in model
#                    mql8_c                    dmso_e No producing reaction in model
#                  phe__L_c                  phe__L_e No producing reaction in model
#                   thmpp_c               2ahethmpp_c No producing reaction in model
#                  trp__L_c                  trp__L_e No producing reaction in model
#                  tyr__L_c                  tyr__L_e No producing reaction in model
# 
# --- Full Trace Report ---
# Blocked_Biomass_Precursor Root_Cause_Gap_Metabolite                          Notes
#                  10fthf_c                     skm_e No producing reaction in model
#                   mlthf_c                     skm_e No producing reaction in model
#                    mql8_c                    dmso_e No producing reaction in model
#                  phe__L_c                  phe__L_e No producing reaction in model
#                     thf_c                     skm_e No producing reaction in model
#                   thmpp_c               2ahethmpp_c No producing reaction in model
#                  trp__L_c                  trp__L_e No producing reaction in model
#                  tyr__L_c                  tyr__L_e No producing reaction in model
# 



# --- 9. LOAD UNIVERSAL FILE MISSING REACTIONS FILE THE ORIGINAL MODEL ---

print(f"Loading universal reaction database: {universal_reactions_path}")
try:
    with open(universal_reactions_path, 'r') as f:
        reaction_db = json.load(f)
    print(f"Loaded {len(reaction_db)} reactions from JSON.")
except Exception as e:
    print(f"FATAL ERROR: Could not read universal reactions JSON file: {e}")
    sys.exit()

print(f"Loading missing reactions from: {missing_reactions_path}")
try:
    missing_reactions_df = pd.read_csv(missing_reactions_path, sep='\t')
    missing_reactions_set = set(missing_reactions_df['Missing_BiGG_Reaction'])
    print(f"Loaded {len(missing_reactions_set)} missing reaction IDs from your genome file.")
except Exception as e:
    print(f"FATAL ERROR: Could not read missing reactions file: {e}")
    sys.exit()

print(f"Loading original model: {original_model_path}")
try:
    model = cobra.io.read_sbml_model(original_model_path)
    # Create a set of all metabolite IDs that already exist in the model
    model_metabolites = {met.id for met in model.metabolites}
    print(f"Successfully loaded model '{model.id}' with {len(model_metabolites)} metabolites.")
except Exception as e:
    print(f"FATAL ERROR: Could not read the model file: {e}")
    sys.exit()



# --- 10. DEFINE TARGETS & HELPER FUNCTION ---
blocked_metabolites = [
    'skm_c',        # Shikimate (Internal)
    'phe__L_c',     # L-Phenylalanine (Internal)
    'trp__L_c',     # L-Tryptophan (Internal)
    'tyr__L_c',     # L-Tyrosine (Internal)
    'mql8_c',       # Menaquinol-8
    '2ahethmpp_c'   # Thiamine precursor
]

blocked_metabolites.extend(['sucbz_c', '2shchc_c', 'sbzcoa_c', 'dhna_c'])

def check_reaction_cost(reaction_data, model_mets):
    """Calculates the 'cost' of adding a reaction."""
    missing_precursors = []
    # Check all reactants (stoichiometry < 0)
    for met_info in reaction_data.get('metabolites', []):
        if met_info.get('stoichiometry', 0) < 0:
            met_id = f"{met_info['bigg_id']}_{met_info['compartment_bigg_id']}"
            if met_id not in model_mets:
                missing_precursors.append(met_id)
    return missing_precursors

# --- 11. RUN THE PROSPECTOR ---
print("--- GSMM_8 PATHWAY PROSPECTOR REPORT ---")
print("Searching BiGG database and cross-referencing with your genome file...")

candidate_pathways = {}

for rxn_id, reaction_data in reaction_db.items():
    for met_info in reaction_data.get('metabolites', []):
        if met_info.get('stoichiometry', 0) > 0:
            met_id_in_db = f"{met_info['bigg_id']}_{met_info['compartment_bigg_id']}"
            
            if met_id_in_db in blocked_metabolites:
                missing_precursors = check_reaction_cost(reaction_data, model_metabolites)
                cost = len(missing_precursors)
                
                # *** NEW LOGIC: Check if it's in the genome file ***
                in_genome = rxn_id in missing_reactions_set
                
                pathway_info = {
                    'id': rxn_id,
                    'name': reaction_data['name'],
                    'equation': reaction_data['reaction_string'],
                    'cost': cost,
                    'missing': missing_precursors,
                    'produces': met_id_in_db,
                    'in_genome': in_genome  # <-- Store the result
                }
                
                if met_id_in_db not in candidate_pathways:
                    candidate_pathways[met_id_in_db] = []
                
                candidate_pathways[met_id_in_db].append(pathway_info)

# --- 12. PRINT THE FINAL REPORT ---
if candidate_pathways:
    print("\n--- Candidate Pathway Options ---")
    print("Found the following reactions in BiGG that could fill the gaps.")
    print("Pathways are sorted by [In Genome] (True first), then by 'cost' (lowest first).\n")
    
    for met_id in blocked_metabolites: # Print in logical order
        if met_id in candidate_pathways:
            print(f"\n--- Options to produce: {met_id} ---")
            
            # *** NEW SORTING LOGIC ***
            # Sorts by:
            # 1. 'in_genome' (True comes first)
            # 2. 'cost' (0 comes first)
            sorted_options = sorted(candidate_pathways[met_id], key=lambda x: (not x['in_genome'], x['cost']))
            
            for option in sorted_options:
                # *** NEW PRINTING LOGIC ***
                print(f"  Reaction: {option['id']} (Cost: {option['cost']}) [In Genome: {option['in_genome']}]")
                print(f"    Name: {option['name']}")
                print(f"    Equation: {option['equation']}")
                if option['cost'] > 0:
                    print(f"    MISSING PRECURSORS: {option['missing']}")
                elif not option['in_genome']:
                    pass # Don't print "all precursors present" if it's not in the genome anyway
                else:
                    print(f"    (This reaction's precursors are ALL present in your model)")
                print("  " + "-"*30)
else:
    print("\n--- No Pathway Found ---")
    print("The script could not find any reaction in the BiGG database")
    print("that produces any of the known blocked precursors.")



# --- Options to produce: skm_c ---
#   Reaction: SHK3Dr (Cost: 0) [In Genome: False]
#     Name: Shikimate dehydrogenase
#     Equation: 3dhsk_c + h_c + nadph_c &#8652; nadp_c + skm_c
#   ------------------------------
# 
# --- Options to produce: phe__L_c ---
#   Reaction: GLYPHEHYc (Cost: 1) [In Genome: True]
#     Name: Hydrolysis of glycylphenylalanine
#     Equation: h2o_c + glyphe_c &#8652; gly_c + phe__L_c
#     MISSING PRECURSORS: ['glyphe_c']
#   ------------------------------
#   Reaction: AROH (Cost: 1) [In Genome: False]
#     Name: L-Arogenate hydro-lyase
#     Equation: Largn_c &#8652; co2_c + h2o_c + phe__L_c
#     MISSING PRECURSORS: ['Largn_c']
#   ------------------------------
# 
# --- Options to produce: trp__L_c ---
#   Reaction: TRPS1 (Cost: 0) [In Genome: False]
#     Name: Tryptophan synthase (indoleglycerol phosphate)
#     Equation: 3ig3p_c + ser__L_c &#8652; g3p_c + h2o_c + trp__L_c
#   ------------------------------
#   Reaction: TRPS2 (Cost: 0) [In Genome: False]
#     Name: Tryptophan synthase (indole)
#     Equation: indole_c + ser__L_c &#8652; h2o_c + trp__L_c
#   ------------------------------
# 
# --- Options to produce: tyr__L_c ---
#   Reaction: TYRTAi (Cost: 0) [In Genome: True]
#     Name: Tyrosine transaminase  irreversible
#     Equation: 34hpp_c + glu__L_c &#8652; akg_c + tyr__L_c
#     (This reaction's precursors are ALL present in your model)
#   ------------------------------
#   Reaction: PHETHPTOX (Cost: 1) [In Genome: False]
#     Name: L-Phenylalanine,tetrahydrobiopterin:oxygen oxidoreductase (4-hydroxylating) Phenylalanine, tyrosine and tryptophan biosynthesis EC:1.14.16.1
#     Equation: o2_c + phe__L_c + thbpt_c &#8652; h2o_c + tyr__L_c + dhbpt_c
#     MISSING PRECURSORS: ['thbpt_c']
#   ------------------------------
# 
# --- Options to produce: mql8_c ---
#   Reaction: AMMQLT8 (Cost: 0) [In Genome: False]
#     Name: S-adenosylmethione:2-demthylmenaquinole methyltransferase (menaquinone 8)
#     Equation: 2dmmql8_c + amet_c &#8652; ahcys_c + h_c + mql8_c
#   ------------------------------
#   Reaction: L_LACD3 (Cost: 0) [In Genome: False]
#     Name: L-Lactate dehydrogenase (menaquinone)
#     Equation: lac__L_c + mqn8_c &#8652; mql8_c + pyr_c
#   ------------------------------
# 
# --- Options to produce: sucbz_c ---
#   Reaction: SUCBZS (Cost: 0) [In Genome: False]
#     Name: O-succinylbenzoate-CoA synthase
#     Equation: 2shchc_c &#8652; h2o_c + sucbz_c
#   ------------------------------
# 
# --- Options to produce: 2shchc_c ---
#   Reaction: 2S6HCCi (Cost: 0) [In Genome: True]
#     Name: 2 succinyl 6 hydroxy 2 4 cyclohexadiene 1 carboxylate synthase
#     Equation: akg_c + h_c + ichor_c &#8652; 2shchc_c + co2_c + pyr_c
#     (This reaction's precursors are ALL present in your model)
#   ------------------------------
#   Reaction: SHCHCS3 (Cost: 0) [In Genome: False]
#     Name: 2-succinyl-6-hydroxy-2,4-cyclohexadiene 1-carboxylate synthase
#     Equation: 2sephchc_c &#8652; 2shchc_c + pyr_c
#   ------------------------------
# 
# --- Options to produce: sbzcoa_c ---
#   Reaction: SUCBZL (Cost: 0) [In Genome: False]
#     Name: O-succinylbenzoate-CoA ligase
#     Equation: atp_c + coa_c + sucbz_c &#8652; amp_c + ppi_c + sbzcoa_c
#   ------------------------------
# 
# --- Options to produce: dhna_c ---
#   Reaction: NPHS (Cost: 0) [In Genome: False]
#     Name: Naphthoate synthase
#     Equation: sbzcoa_c &#8652; coa_c + dhna_c
#   ------------------------------



from cobra import Reaction, Metabolite

# --- 13. ADD CURATION REACTIONS (based on Prospector Report) ---

print("\n--- Applying Manual Gap-Filling ---")

# --- Get Existing Metabolites (or create placeholders) ---
# Helper function to keep the code clean
def get_met(met_id, formula, name):
    try:
        return model.metabolites.get_by_id(met_id)
    except KeyError:
        return Metabolite(met_id, formula=formula, name=name, compartment='C_c')

# Core Metabolites
akg_c   = get_met('akg_c', 'C5H4O5', '2-Oxoglutarate')
ichor_c = get_met('ichor_c', 'C10H8O6', 'Isochorismate')
h_c     = get_met('h_c', 'H', 'H+')
co2_c   = get_met('co2_c', 'CO2', 'CO2')
pyr_c   = get_met('pyr_c', 'C3H3O3', 'Pyruvate')
h2o_c   = get_met('h2o_c', 'H2O', 'Water')
atp_c   = get_met('atp_c', 'C10H12N5O13P3', 'ATP')
coa_c   = get_met('coa_c', 'C21H32N7O16P3S', 'CoA')
amp_c   = get_met('amp_c', 'C10H12N5O7P', 'AMP')
ppi_c   = get_met('ppi_c', 'HO7P2', 'Diphosphate')
adp_c   = get_met('adp_c', 'C10H12N5O10P2', 'ADP')
pi_c    = get_met('pi_c', 'HO4P', 'Phosphate')

# Shikimate specific
met_3dhsk_c = get_met('3dhsk_c', 'C7H7O5', '3-Dehydroshikimate')
nadph_c     = get_met('nadph_c', 'C21H26N7O17P3', 'NADPH')
nadp_c      = get_met('nadp_c', 'C21H25N7O17P3', 'NADP+')
skm_c       = get_met('skm_c', 'C7H9O5', 'Shikimate')

# Menaquinone Intermediates (Fixed formula for 2shchc_c)
met_2shchc_c = get_met('2shchc_c', 'C11H10O6', '2-Succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate')
met_sucbz_c  = get_met('sucbz_c', 'C11H8O5', 'O-Succinylbenzoate')
met_sbzcoa_c = get_met('sbzcoa_c', 'C32H39N7O20P3S', 'O-Succinylbenzoyl-CoA')
dhna_c       = get_met('dhna_c', 'C11H7O4', '1,4-dihydroxy-2-naphthoate')
_2dmmql8_c   = get_met('2dmmql8_c', 'C50H74O2', '2-Demethylmenaquinol 8')
amet_c       = get_met('amet_c', 'C15H23N6O5S', 'S-Adenosyl-L-methionine')
ahcys_c      = get_met('ahcys_c', 'C14H20N6O5S', 'S-Adenosyl-L-homocysteine')
mql8_c       = get_met('mql8_c', 'C51H76O2', 'Menaquinol 8')

# Aromatic Amino Acids Intermediates
glu__L_c = get_met('glu__L_c', 'C5H8NO4', 'L-Glutamate')
_34hpp_c = get_met('34hpp_c', 'C9H7O4', '3-(4-Hydroxyphenyl)pyruvate')
tyr__L_c = get_met('tyr__L_c', 'C9H11NO3', 'L-Tyrosine')
ser__L_c = get_met('ser__L_c', 'C3H7NO3', 'L-Serine')
_3ig3p_c = get_met('3ig3p_c', 'C11H12NO6P', 'Indole-3-glycerol phosphate')
g3p_c    = get_met('g3p_c', 'C3H5O6P', 'Glyceraldehyde 3-phosphate')
trp__L_c = get_met('trp__L_c', 'C11H12N2O2', 'L-Tryptophan')


# --- DEFINE REACTIONS ---

# 1. 2S6HCCi (Synthase)
rxn_2s6hcci = Reaction('2S6HCCi')
rxn_2s6hcci.name = '2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase'
rxn_2s6hcci.lower_bound = -1000.0  # Reversible per equation
rxn_2s6hcci.upper_bound = 1000.0
rxn_2s6hcci.add_metabolites({akg_c: -1.0, h_c: -1.0, ichor_c: -1.0, met_2shchc_c: 1.0, co2_c: 1.0, pyr_c: 1.0})

# 2. SUCBZS (Synthase)
rxn_sucbzs = Reaction('SUCBZS')
rxn_sucbzs.name = 'O-succinylbenzoate-CoA synthase'
rxn_sucbzs.lower_bound = -1000.0
rxn_sucbzs.upper_bound = 1000.0
rxn_sucbzs.add_metabolites({met_2shchc_c: -1.0, h2o_c: 1.0, met_sucbz_c: 1.0})

# 3. SUCBZL (Ligase)
rxn_sucbzl = Reaction('SUCBZL')
rxn_sucbzl.name = 'O-succinylbenzoate-CoA ligase'
rxn_sucbzl.lower_bound = -1000.0
rxn_sucbzl.upper_bound = 1000.0
rxn_sucbzl.add_metabolites({atp_c: -1.0, coa_c: -1.0, met_sucbz_c: -1.0, amp_c: 1.0, ppi_c: 1.0, met_sbzcoa_c: 1.0})

# 4. SHK3Dr (Shikimate Dehydrogenase)
rxn_shk3dr = Reaction('SHK3Dr')
rxn_shk3dr.name = 'Shikimate dehydrogenase'
rxn_shk3dr.lower_bound = -1000.0 
rxn_shk3dr.upper_bound = 1000.0
rxn_shk3dr.add_metabolites({met_3dhsk_c: -1.0, h_c: -1.0, nadph_c: -1.0, skm_c: 1.0, nadp_c: 1.0})

# 5. ADK1 (Adenylate Kinase)
rxn_adk1 = Reaction('ADK1')
rxn_adk1.name = 'Adenylate kinase'
rxn_adk1.lower_bound = -1000.0 
rxn_adk1.upper_bound = 1000.0
rxn_adk1.add_metabolites({amp_c: -1.0, atp_c: -1.0, adp_c: 2.0})

# 6. PPA (Inorganic Pyrophosphatase)
rxn_ppa = Reaction('PPA')
rxn_ppa.name = 'Inorganic pyrophosphatase'
rxn_ppa.lower_bound = 0.0 
rxn_ppa.upper_bound = 1000.0
rxn_ppa.add_metabolites({ppi_c: -1.0, h2o_c: -1.0, pi_c: 2.0, h_c: 1.0})

# 7. TYRTAi (Tyrosine transaminase - In Genome)
rxn_tyrtai = Reaction('TYRTAi')
rxn_tyrtai.name = 'Tyrosine transaminase irreversible'
rxn_tyrtai.lower_bound = -1000.0 
rxn_tyrtai.upper_bound = 1000.0
rxn_tyrtai.add_metabolites({_34hpp_c: -1.0, glu__L_c: -1.0, akg_c: 1.0, tyr__L_c: 1.0})

# 8. TRPS1 (Tryptophan synthase)
rxn_trps1 = Reaction('TRPS1')
rxn_trps1.name = 'Tryptophan synthase (indoleglycerol phosphate)'
rxn_trps1.lower_bound = -1000.0 
rxn_trps1.upper_bound = 1000.0
rxn_trps1.add_metabolites({_3ig3p_c: -1.0, ser__L_c: -1.0, g3p_c: 1.0, h2o_c: 1.0, trp__L_c: 1.0})

# 9. NPHS (Naphthoate synthase)
rxn_nphs = Reaction('NPHS')
rxn_nphs.name = 'Naphthoate synthase'
rxn_nphs.lower_bound = -1000.0 
rxn_nphs.upper_bound = 1000.0
rxn_nphs.add_metabolites({met_sbzcoa_c: -1.0, coa_c: 1.0, dhna_c: 1.0})

# 10. AMMQLT8 (Menaquinone 8 final step)
rxn_ammqlt8 = Reaction('AMMQLT8')
rxn_ammqlt8.name = 'S-adenosylmethione:2-demthylmenaquinole methyltransferase'
rxn_ammqlt8.lower_bound = -1000.0 
rxn_ammqlt8.upper_bound = 1000.0
rxn_ammqlt8.add_metabolites({_2dmmql8_c: -1.0, amet_c: -1.0, ahcys_c: 1.0, h_c: 1.0, mql8_c: 1.0})


# --- 14. ADD TO MODEL ---
reactions_to_add = [
    rxn_2s6hcci, rxn_sucbzs, rxn_sucbzl, rxn_shk3dr, rxn_adk1, rxn_ppa,
    rxn_tyrtai, rxn_trps1, rxn_nphs, rxn_ammqlt8
]

print(f"Adding {len(reactions_to_add)} reactions to model...")

for rxn in reactions_to_add:
    if rxn.id not in model.reactions:
        model.add_reactions([rxn])
        print(f" + Added {rxn.id}")
    else:
        print(f" ! {rxn.id} already exists, skipping.")

# Adding 10 reactions to model...
#  + Added 2S6HCCi
#  ! SUCBZS already exists, skipping.
#  ! SUCBZL already exists, skipping.
#  + Added SHK3Dr
#  ! ADK1 already exists, skipping.
#  ! PPA already exists, skipping.
#  + Added TYRTAi
#  ! TRPS1 already exists, skipping.
#  ! NPHS already exists, skipping.
#  ! AMMQLT8 already exists, skipping.
 

# --- 15. TEST GROWTH ---
set_minimal_medium(model, full_minimal_medium)

print("\n--- Testing growth on gap-filled model ---")
try:
    growth_reaction = model.reactions.get_by_id('Growth')
    model.objective = growth_reaction
    solution = model.optimize()
    if solution.objective_value > 0.01:
        print(f"\n--- SUCCESS! ---")
        print(f"Model now grows! Growth rate: {solution.objective_value:.4f}")
        print("Pathways restored successfully.")
        print(f"Saving functional model to: {output_model_path}")
        cobra.io.write_sbml_model(model, output_model_path)
    else:
        print(f"--- FAILED. Growth rate: {solution.objective_value} ---")
except Exception as e:
    print(f"Error during optimization: {e}")

# =============================================================================
# --- SUCCESS! ---
# Model now grows! Growth rate: 0.8931
# =============================================================================



# --- 16. VALIDATION SUITE (Looping through Carbon Sources) ---

print("--- RUNNING CARBON SOURCE VALIDATION SUITE ---")

try:
    model = cobra.io.read_sbml_model(output_model_path)
    print(f"Successfully loaded model '{model.id}'")
except Exception as e:
    print(f"FATAL ERROR: Could not read the model file: {e}")
    sys.exit()

try:
    growth_reaction = model.reactions.get_by_id('Growth')
    print(f"Found growth reaction: '{growth_reaction.id}'")
except KeyError:
    print(f"Warning: Could not find 'Growth'. Trying 'R_Growth'.")
    try:
        growth_reaction = model.reactions.get_by_id('R_Growth')
        print(f"Found growth reaction: '{growth_reaction.id}'")
    except KeyError:
        print(f"FATAL ERROR: Could not find 'Growth' or 'R_Growth' reaction.")
        sys.exit()

# 1. Define the BASE minimal medium (Everything EXCEPT the carbon source)
base_medium = {
    "EX_nh4_e": 10.0, "EX_so4_e": 10.0, "EX_k_e": 10.0, "EX_pi_e": 10.0,
    "EX_na1_e": 10.0, "EX_mg2_e": 10.0, "EX_h2o_e": 10.0, "EX_o2_e": 10.0
}

trace_elements_exchanges = [
    'EX_fe2_e', 'EX_fe3_e', 'EX_zn2_e', 'EX_mn2_e', 'EX_cu2_e',
    'EX_cobalt2_e', 'EX_cl_e', 'EX_ca2_e', 'EX_mobd_e', 'EX_ni2_e',
    'EX_sel_e', 'EX_slnt_e', 'EX_tungs_e', 'EX_h_e'
]

# Add trace elements to the base medium
for ex_id in trace_elements_exchanges:
    base_medium[ex_id] = 10.0

# 2. Define the conditions to test
# Format: "Condition Name": {"Exchange_Reaction": Uptake_Rate}
carbon_tests = {
    "No Carbon": {},  # Empty dict means no carbon added
    "Sucrose":   {"EX_sucr_e": 10.0},
    "Fructose":  {"EX_fru_e": 10.0},
    "Glucose":   {"EX_glc__D_e": 10.0},
    "Sorbitol":  {"EX_sbt__D_e": 10.0}
}

# Dictionary to store results for the summary
validation_results = {}

# 3. Run the loop
for condition_name, carbon_source in carbon_tests.items():
    print(f"\n--- Testing growth with: {condition_name} ---")
    
    # Create a fresh copy of the base medium for this test
    current_medium = base_medium.copy()
    
    # Add the specific carbon source for this loop iteration
    current_medium.update(carbon_source)
    
    # Set the medium (this zeroes out all other exchanges first)
    set_minimal_medium(model, current_medium)
    
    # Optimize
    try:
        model.objective = model.reactions.get_by_id('Growth')
        solution = model.optimize()
        rate = solution.objective_value
        
        # Store the result
        validation_results[condition_name] = rate
        
        # Print iteration outcome
        if rate > 0.01:
            print(f"  -> SUCCESS! Growth rate: {rate:.4f}")
        else:
            print(f"  -> FAILED! No growth (Rate: {rate:.4f})")
            
    except Exception as e:
        print(f"  -> ERROR during optimization: {e}")
        validation_results[condition_name] = "Error"

print("--- VALIDATION SUMMARY ---")
print(f"{'Carbon Source':<15} | {'Predicted Growth Rate'}")
print("-" * 40)
for condition, rate in validation_results.items():
    if isinstance(rate, float):
        print(f"{condition:<15} | {rate:.4f}")
    else:
        print(f"{condition:<15} | {rate}")


# =============================================================================
# --- VALIDATION SUMMARY ---
# Carbon Source   | Predicted Growth Rate
# ----------------------------------------
# No Carbon       | 0.0000
# Sucrose         | 0.8931
# Fructose        | 0.5863
# Glucose         | 0.6648
# Sorbitol        | 0.6805
# =============================================================================

