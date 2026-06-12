import numpy as np
import json
import warnings
from pycalphad import Database, equilibrium, variables as v

warnings.filterwarnings('ignore')

db = Database('DeKeyzer2009_FeNiTi.tdb')
comps = ['FE', 'NI', 'TI', 'VA']

X_NI = 0.172
X_TI = 0.008
composition = {v.X('NI'): X_NI, v.X('TI'): X_TI}

# Temperature grid: 50 to 650 C, steps of 5 C
T_C_grid = list(range(50, 651, 5))
data = {}

print("Starting precomputation...")

for T_C in T_C_grid:
    T_K = T_C + 273.15
    grid_point = {}
    
    # 1. BCC + FCC + ETA + LIQUID equilibrium (Fig 2 / bulk equilibrium)
    try:
        eq_bulk = equilibrium(db, comps, ['BCC_A2', 'FCC_A1', 'ETA', 'LIQUID'], {v.T: T_K, v.P: 101325, **composition})
        ph_names = np.array(eq_bulk.Phase.values.squeeze(), dtype=str)
        fracs = np.array(eq_bulk.NP.values.squeeze(), dtype=float)
        
        grid_point['bulk_f_BCC'] = 0.0
        grid_point['bulk_f_FCC'] = 0.0
        grid_point['bulk_f_ETA'] = 0.0
        
        for i in range(len(ph_names)):
            ph = ph_names[i].strip()
            fr = fracs[i]
            if not ph or np.isnan(fr) or fr < 1e-12:
                continue
            if ph == 'BCC_A2':
                grid_point['bulk_f_BCC'] = float(fr)
            elif ph == 'FCC_A1':
                grid_point['bulk_f_FCC'] = float(fr)
            elif ph == 'ETA':
                grid_point['bulk_f_ETA'] = float(fr)
    except Exception as e:
        print(f"Error in bulk eq at {T_C} C: {e}")
        
    # 2. BCC + ETA equilibrium (precipitate solubility limit)
    try:
        eq_precip = equilibrium(db, comps, ['BCC_A2', 'ETA'], {v.T: T_K, v.P: 101325, **composition})
        ph_names = np.array(eq_precip.Phase.values.squeeze(), dtype=str)
        fracs = np.array(eq_precip.NP.values.squeeze(), dtype=float)
        
        grid_point['precip_C_Ti_BCC'] = 1.8e-5 
        grid_point['precip_C_Ni_BCC'] = X_NI
        grid_point['precip_C_Ti_ETA'] = 0.25
        grid_point['precip_C_Ni_ETA'] = 0.75
        grid_point['precip_f_ETA'] = 0.0
        
        for i in range(len(ph_names)):
            ph = ph_names[i].strip()
            fr = fracs[i]
            if ph == 'BCC_A2':
                grid_point['precip_C_Ti_BCC'] = float(eq_precip.X.sel(component='TI').values.squeeze()[i])
                grid_point['precip_C_Ni_BCC'] = float(eq_precip.X.sel(component='NI').values.squeeze()[i])
            elif ph == 'ETA':
                grid_point['precip_C_Ti_ETA'] = float(eq_precip.X.sel(component='TI').values.squeeze()[i])
                grid_point['precip_C_Ni_ETA'] = float(eq_precip.X.sel(component='NI').values.squeeze()[i])
                grid_point['precip_f_ETA'] = float(fr)
    except Exception as e:
        print(f"Error in bcc-eta eq at {T_C} C: {e}")
        
    # 3. BCC + FCC equilibrium (austenite interface limits)
    try:
        eq_aust = equilibrium(db, comps, ['BCC_A2', 'FCC_A1'], {v.T: T_K, v.P: 101325, **composition})
        ph_names = np.array(eq_aust.Phase.values.squeeze(), dtype=str)
        fracs = np.array(eq_aust.NP.values.squeeze(), dtype=float)
        
        grid_point['aust_C_Ni_BCC'] = X_NI
        grid_point['aust_C_Ni_FCC'] = X_NI
        grid_point['aust_C_Ti_BCC'] = X_TI
        grid_point['aust_C_Ti_FCC'] = X_TI
        grid_point['aust_f_FCC'] = 0.0
        
        for i in range(len(ph_names)):
            ph = ph_names[i].strip()
            fr = fracs[i]
            if ph == 'BCC_A2':
                grid_point['aust_C_Ni_BCC'] = float(eq_aust.X.sel(component='NI').values.squeeze()[i])
                grid_point['aust_C_Ti_BCC'] = float(eq_aust.X.sel(component='TI').values.squeeze()[i])
            elif ph == 'FCC_A1':
                grid_point['aust_C_Ni_FCC'] = float(eq_aust.X.sel(component='NI').values.squeeze()[i])
                grid_point['aust_C_Ti_FCC'] = float(eq_aust.X.sel(component='TI').values.squeeze()[i])
                grid_point['aust_f_FCC'] = float(fr)
    except Exception as e:
        print(f"Error in bcc-fcc eq at {T_C} C: {e}")
        
    data[str(T_C)] = grid_point
    print(f"  Precomputed T = {T_C} C")

# ================================================================
#  FULL EQUILIBRIUM SWEEP FOR FIG 2 PHASE DIAGRAM (0 - 1750 C)
# ================================================================
print("\nComputing full equilibrium sweep for Fig 2 (0 - 1750 C)...")
T_sweep_C = np.linspace(0, 1750, 36)
T_sweep_K = T_sweep_C + 273.15

fig2_data = {'T_C': [], 'f_ETA': [], 'f_FCC': []}
for tC, tK in zip(T_sweep_C, T_sweep_K):
    try:
        eq_full = equilibrium(db, comps, ['BCC_A2', 'FCC_A1', 'ETA', 'LIQUID'], {v.T: tK, v.P: 101325, **composition})
        ph_names = np.array(eq_full.Phase.values.squeeze(), dtype=str)
        fracs = np.array(eq_full.NP.values.squeeze(), dtype=float)
        
        f_eta = 0.0
        f_fcc = 0.0
        for i, ph in enumerate(ph_names):
            ph = ph.strip()
            if ph == 'ETA':
                f_eta = float(fracs[i])
            elif ph == 'FCC_A1':
                f_fcc = float(fracs[i])
                
        fig2_data['T_C'].append(float(tC))
        fig2_data['f_ETA'].append(f_eta)
        fig2_data['f_FCC'].append(f_fcc)
    except Exception as e:
        print(f"Error computing full eq at {tC} C: {e}")

data['fig2_full_equilibrium'] = fig2_data

with open('precomputed_thermo.json', 'w') as f:
    json.dump(data, f, indent=2)
print("\nSaved precomputed_thermo.json successfully (including Fig 2 data).")
