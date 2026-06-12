"""
Full Maraging Steel Precipitation Hardening Pipeline — CORRECTED
Forces BCC (A2) matrix (metastable martensite), excludes FCC (A1).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import json, os, warnings

print("="*70)
print("  MARAGING STEEL M300 — FULL PRECIPITATION HARDENING PIPELINE")
print("  (Metastable martensite: A2 matrix, A1 excluded)")
print("="*70)

from pycalphad import Database, equilibrium, variables as v
print("[OK] pycalphad imported.")

warnings.filterwarnings('ignore')
plt.rcParams.update({
    'figure.figsize': (10, 6), 'figure.dpi': 150, 'font.size': 12,
    'axes.labelsize': 13, 'axes.titlesize': 14, 'legend.fontsize': 10,
    'lines.linewidth': 2.0, 'figure.constrained_layout.use': True
})
COLORS = {400: '#264653', 475: '#2A9D8F', 500: '#E9C46A', 600: '#E63946'}
R_gas = 8.314
P_atm = 101325

def extract_eq_data(eq_result):
    phases = np.array(eq_result.Phase.values.squeeze(), dtype=str)
    fracs  = np.array(eq_result.NP.values.squeeze(), dtype=float)
    data = {}
    for i in range(len(phases)):
        ph = phases[i].strip()
        fr = fracs[i]
        if not ph or np.isnan(fr) or fr < 1e-12:
            continue
        x = {}
        for c in ['FE', 'NI', 'TI']:
            x[c] = float(eq_result.X.sel(component=c).values.squeeze()[i])
        data[ph] = {'f': fr, 'X': x}
    return data

# ================================================================
# SECTION 1 — PyCalphad Thermodynamics
# ================================================================
print("\n[1/10] Loading TDB ...")
db = Database('DeKeyzer2009_FeNiTi.tdb')
comps = ['FE', 'NI', 'TI', 'VA']

all_phases = sorted(db.phases.keys())
print(f"  All phases in DB: {all_phases}")

# Auto-detect phase names (supports both old and new TDB naming)
# BCC martensite: 'BCC_A2' or 'A2'
if 'BCC_A2' in all_phases:
    PH_BCC = 'BCC_A2'
elif 'A2' in all_phases:
    PH_BCC = 'A2'
else:
    raise RuntimeError("No BCC phase found in TDB! Expected 'BCC_A2' or 'A2'.")

# Ni3Ti precipitate (eta, D024): 'ETA' or 'NI3TI'
if 'ETA' in all_phases:
    PH_ETA = 'ETA'
elif 'NI3TI' in all_phases:
    PH_ETA = 'NI3TI'
else:
    raise RuntimeError("No Ni3Ti/eta phase found in TDB! Expected 'ETA' or 'NI3TI'.")

# FCC austenite: 'FCC_A1' or 'A1'
if 'FCC_A1' in all_phases:
    PH_FCC = 'FCC_A1'
elif 'A1' in all_phases:
    PH_FCC = 'A1'
else:
    PH_FCC = None
    print("  WARNING: No FCC phase found in TDB. Fig 2B (austenite) will be skipped.")

# Metastable phases: BCC + Ni3Ti only (excludes FCC for kinetics)
# Martensite is metastable -- excluding FCC (austenite) forces the
# thermodynamic calculation to use BCC as the matrix.
my_phases = [PH_BCC, PH_ETA, 'LIQUID']

# Full equilibrium phases: BCC + FCC + Ni3Ti + LIQUID (for Fig 2)
if PH_FCC:
    full_eq_phases = [PH_BCC, PH_FCC, PH_ETA, 'LIQUID']
else:
    full_eq_phases = my_phases

X_NI = 0.172
X_TI = 0.008
composition = {v.X('NI'): X_NI, v.X('TI'): X_TI}

print(f"  Detected BCC phase: {PH_BCC}")
print(f"  Detected Ni3Ti phase: {PH_ETA}")
print(f"  Detected FCC phase: {PH_FCC}")
print(f"  Metastable phases (kinetics): {my_phases}")
print(f"  Full equil. phases (Fig 2):   {full_eq_phases}")
print(f"  NOTE: FCC excluded from kinetics -- martensite is metastable BCC.")
print(f"  Composition: X(Ni)={X_NI}, X(Ti)={X_TI}, X(Fe)={1-X_NI-X_TI:.3f}")

# ---- Single-point equilibrium at 755 K ----
print("\n[2/10] Equilibrium at 755 K (482 C) ...")
T_age = 755
eq_755 = equilibrium(db, comps, my_phases,
                     {v.T: T_age, v.P: P_atm, **composition})
data_755 = extract_eq_data(eq_755)

print("=" * 70)
for ph, info in sorted(data_755.items()):
    print(f"  {ph:8s}  f = {info['f']:.6f}  "
          f"X(Fe)={info['X']['FE']:.5f}  "
          f"X(Ni)={info['X']['NI']:.5f}  "
          f"X(Ti)={info['X']['TI']:.5f}")
print("=" * 70)

f_eq_755 = data_755.get(PH_ETA, {}).get('f', 0.0)
_matrix_ph = PH_BCC
if PH_BCC not in data_755:
    _matrix_ph = max(data_755.keys(), key=lambda p: data_755[p]['f'])
    print(f"  WARNING: {PH_BCC} not found, using {_matrix_ph}")

x_ni_matrix = data_755[_matrix_ph]['X']['NI']
x_ti_matrix = data_755[_matrix_ph]['X']['TI']
x_ni_eta    = data_755.get(PH_ETA, {}).get('X', {}).get('NI', float('nan'))
x_ti_eta    = data_755.get(PH_ETA, {}).get('X', {}).get('TI', float('nan'))

mu_vals = eq_755.MU.values.squeeze()
comp_labels = [str(c) for c in eq_755.coords['component'].values if str(c) != 'VA']
print(f"\nMatrix phase: {_matrix_ph}")
print("Chemical potentials:")
for c, mu in zip(comp_labels, mu_vals):
    print(f"  mu({c}) = {mu:>12.2f} J/mol")

print(f"\n--- KEY RESULTS ---")
print(f"  f_eq({PH_ETA}) = {f_eq_755:.4f}  ({f_eq_755*100:.2f} %)")
if 0.05 < f_eq_755 < 0.20:
    print("  [OK] Within expected range 5-20%.")
elif 0.02 < f_eq_755 < 0.05:
    print("  [NOTE] Below expected 10-15%; ternary approx may underestimate.")
else:
    print(f"  [WARN] Outside expected range.")
print(f"  X(Ni)_matrix = {x_ni_matrix:.6f}")
print(f"  X(Ti)_matrix = {x_ti_matrix:.6f}")
print(f"  X(Ni)_eta    = {x_ni_eta:.6f}")
print(f"  X(Ti)_eta    = {x_ti_eta:.6f}")

# ---- Temperature sweep ----
print("\n[3/10] Temperature sweep (400-1800 K, 70 pts) ...")
T_range = np.linspace(400, 1800, 70)
sweep_T, sweep_fNI3TI, sweep_fA2 = [], [], []

for i, T_val in enumerate(T_range):
    try:
        eq = equilibrium(db, comps, my_phases,
                         {v.T: float(T_val), v.P: P_atm, **composition})
        d = extract_eq_data(eq)
    except Exception:
        d = {}
    sweep_T.append(float(T_val))
    sweep_fNI3TI.append(d.get(PH_ETA, {}).get('f', 0.0))
    sweep_fA2.append(d.get(PH_BCC, {}).get('f', 0.0))
    if (i+1) % 10 == 0:
        print(f"  ... {i+1}/{len(T_range)} done")

sweep_T      = np.array(sweep_T)
sweep_fNI3TI = np.array(sweep_fNI3TI)
sweep_fA2    = np.array(sweep_fA2)

df_sweep = pd.DataFrame({
    'T_K': sweep_T, 'T_C': sweep_T - 273.15,
    'f_NI3TI': sweep_fNI3TI, 'f_A2': sweep_fA2
})

# Solvus detection
threshold = 1e-4
above = sweep_fNI3TI > threshold
if np.any(above):
    idx_last = np.where(above)[0][-1]
    if idx_last + 1 < len(sweep_T):
        T1, T2 = sweep_T[idx_last], sweep_T[idx_last + 1]
        f1, f2 = sweep_fNI3TI[idx_last], sweep_fNI3TI[idx_last + 1]
        solvus_T_K = T1 + (threshold - f1) * (T2 - T1) / (f2 - f1) if f1 != f2 else T1
    else:
        solvus_T_K = float(sweep_T[idx_last])
else:
    solvus_T_K = float('nan')
solvus_T_C = solvus_T_K - 273.15

print(f"\nNi3Ti solvus temperature: {solvus_T_K:.1f} K ({solvus_T_C:.1f} C)")
if 1073 < solvus_T_K < 1323:
    print("[OK] Within expected 800-1050 C range.")
else:
    print(f"[NOTE] Solvus = {solvus_T_C:.0f} C (target 900-950 C for full alloy)")

# ---- PLOT 1 ----
print("\n  Generating Plot 1 ...")
fig, ax1 = plt.subplots(figsize=(11, 6))
ax1.plot(sweep_T, sweep_fNI3TI, 'o-', color='#E63946', ms=4, label=f'{PH_ETA} (eta, D024)')
ax1.plot(sweep_T, sweep_fA2,    's-', color='#457B9D', ms=3, alpha=0.8, label=f'{PH_BCC} (BCC martensite)')
if not np.isnan(solvus_T_K):
    ax1.axvline(solvus_T_K, ls='--', color='#E63946', alpha=0.6,
                label=f'Solvus = {solvus_T_K:.0f} K ({solvus_T_C:.0f} C)')
ax1.axvline(755, ls=':', color='#E9C46A', alpha=0.8, label='482 C aging temp')
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Equilibrium phase fraction')
ax1.set_title(f'Plot 1 - Phase fractions vs Temperature (cf. Fig 2A, Ahluwalia 2024)\n'
              f'Fe-{X_NI*100:.1f}Ni-{X_TI*100:.1f}Ti (at%), De Keyzer (2009) TDB, '
              f'metastable BCC matrix')
ax1.legend(loc='center right', framealpha=0.9)
ax1.set_xlim(400, 1800); ax1.set_ylim(-0.02, 1.05)
ax1.grid(True, alpha=0.25)
ax1b = ax1.twiny()
ax1b.set_xlim(400 - 273.15, 1800 - 273.15)
ax1b.set_xlabel('Temperature (C)')
plt.savefig('plot1_phase_fraction_vs_T.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: plot1_phase_fraction_vs_T.png")

# ============================================================
# FIG 2 — Full equilibrium phase fractions vs Temperature
# (BCC + FCC + Ni3Ti + LIQUID — true equilibrium, not metastable)
# ============================================================
print("\n[3b/10] Full equilibrium sweep for Fig 2 (0-1750 C, 100 pts) ...")
T_fig2_C = np.linspace(0, 1750, 100)
T_fig2_K = T_fig2_C + 273.15
fig2_fETA, fig2_fFCC, fig2_fBCC = [], [], []

for i, T_val in enumerate(T_fig2_K):
    try:
        eq = equilibrium(db, comps, full_eq_phases,
                         {v.T: float(T_val), v.P: P_atm, **composition})
        d = extract_eq_data(eq)
    except Exception:
        d = {}
    fig2_fETA.append(d.get(PH_ETA, {}).get('f', 0.0))
    fig2_fFCC.append(d.get(PH_FCC, {}).get('f', 0.0) if PH_FCC else 0.0)
    fig2_fBCC.append(d.get(PH_BCC, {}).get('f', 0.0))
    if (i+1) % 20 == 0:
        print(f"  ... {i+1}/{len(T_fig2_K)} done")

fig2_fETA = np.array(fig2_fETA)
fig2_fFCC = np.array(fig2_fFCC)
fig2_fBCC = np.array(fig2_fBCC)


# ---- Aging temperature table ----
print("\n[4/10] Aging temperature equilibrium table ...")
aging_T_C_list = [400, 450, 475, 480, 500, 510, 540, 600]
records = []
for T_C in aging_T_C_list:
    T_K = T_C + 273.15
    try:
        eq = equilibrium(db, comps, my_phases,
                         {v.T: T_K, v.P: P_atm, **composition})
        d = extract_eq_data(eq)
    except Exception:
        d = {}
    f_ni3ti = d.get(PH_ETA, {}).get('f', 0.0)
    mat = PH_BCC if PH_BCC in d else None
    x_ni_m = d[mat]['X']['NI'] if mat else float('nan')
    x_ti_m = d[mat]['X']['TI'] if mat else float('nan')
    records.append({'T_C': T_C, 'T_K': T_K, 'f_eq_Ni3Ti': f_ni3ti,
                    'matrix': mat or 'N/A',
                    'X_Ni_matrix': x_ni_m, 'X_Ti_matrix': x_ti_m})

df_aging = pd.DataFrame(records)
df_aging.to_csv('aging_equilibrium.csv', index=False, float_format='%.8f')
print(df_aging.to_string(index=False, float_format='{:.6f}'.format))
print("  Saved: aging_equilibrium.csv")

# ---- Full equilibrium (BCC+FCC+ETA) at aging temperatures ----
# This replaces the hardcoded f_gamma_sat in update_kinetics.py.
# Including FCC in equilibrium gives the true reverted austenite fraction.
print("\n[4b/10] Full equilibrium (BCC+FCC+ETA) at aging temperatures ...")
f_gamma_eq_at_aging_T = {}
f_eta_full_eq_at_aging_T = {}  # Ni3Ti fraction from full equil (BCC+FCC+ETA)
for T_C in [400, 475, 500, 600]:
    T_K = T_C + 273.15
    try:
        eq_full = equilibrium(db, comps, full_eq_phases,
                              {v.T: T_K, v.P: P_atm, **composition})
        d_full = extract_eq_data(eq_full)
        f_fcc  = d_full.get(PH_FCC, {}).get('f', 0.0) if PH_FCC else 0.0
        f_eta  = d_full.get(PH_ETA, {}).get('f', 0.0)
    except Exception as ex:
        print(f"    WARNING at T={T_C} C: {ex}")
        f_fcc, f_eta = 0.0, 0.0
    f_gamma_eq_at_aging_T[T_C]    = float(f_fcc)
    f_eta_full_eq_at_aging_T[T_C] = float(f_eta)
    print(f"  T={T_C} C: f_gamma(FCC)={f_fcc:.5f}  f_eta(ETA)={f_eta:.5f}")
print("  Stored in JSON as 'f_gamma_eq_at_aging_T' (eliminates hardcoded f_gamma_sat)")

# ---- PLOT 2 ----
print("\n  Generating Plot 2 ...")
fig2, ax2 = plt.subplots(figsize=(10, 5))
bars = ax2.bar([str(t) for t in df_aging['T_C']],
               df_aging['f_eq_Ni3Ti'] * 100,
               color='#457B9D', edgecolor='#264653', linewidth=1.2)
for bar, val in zip(bars, df_aging['f_eq_Ni3Ti']):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
             f'{val*100:.2f}%', ha='center', va='bottom', fontsize=9)
ax2.set_xlabel('Aging temperature (C)')
ax2.set_ylabel('Equilibrium Ni3Ti fraction (%)')
ax2.set_title('Plot 2 - Equilibrium Ni3Ti fraction at aging temps\n'
              '(data behind Fig 8A reference lines, Ahluwalia 2024)')
ax2.grid(axis='y', alpha=0.3)
plt.savefig('plot2_feq_barchart.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: plot2_feq_barchart.png")

# ---- Driving force ----
print("\n[5/10] Driving force for nucleation at 755 K ...")
Vm = 6.7e-6
alloy_X = {'NI': X_NI, 'TI': X_TI, 'FE': 1 - X_NI - X_TI}
eq_X    = {'NI': x_ni_matrix, 'TI': x_ti_matrix,
           'FE': 1 - x_ni_matrix - x_ti_matrix}
eta_X   = {'NI': x_ni_eta,    'TI': x_ti_eta,
           'FE': 1 - x_ni_eta - x_ti_eta}

sum_term = 0.0
for c in ['NI', 'TI', 'FE']:
    ratio = alloy_X[c] / eq_X[c] if eq_X[c] > 1e-12 else 1.0
    ln_ratio = np.log(max(ratio, 1e-30))
    contrib = eta_X[c] * ln_ratio
    sum_term += contrib
    print(f"  {c}: X_eta={eta_X[c]:.5f}  X_alloy={alloy_X[c]:.5f}  "
          f"X_eq={eq_X[c]:.5f}  ln(ratio)={ln_ratio:+.4f}  contrib={contrib:+.6f}")

dG_mol = R_gas * T_age * sum_term
dG_vol = dG_mol / Vm
print(f"\n  Driving force dG_nuc at {T_age} K:")
print(f"    Molar:      {dG_mol:>10.1f} J/mol  ({dG_mol/1000:.2f} kJ/mol)")
print(f"    Volumetric: {dG_vol:>10.2e} J/m3")

# ---- Save pipeline_inputs ----
pipeline_inputs = {
    'solvus_T_K':                round(float(solvus_T_K), 2),
    'solvus_T_C':                round(float(solvus_T_C), 2),
    'f_eq_755K':                 round(float(f_eq_755), 6),
    'X_Ni_matrix_755K':          round(float(x_ni_matrix), 6),
    'X_Ti_matrix_755K':          round(float(x_ti_matrix), 6),
    'X_Ni_eta':                  round(float(x_ni_eta), 6),
    'X_Ti_eta':                  round(float(x_ti_eta), 6),
    'driving_force_755K_J_per_mol': round(float(dG_mol), 2)
}
with open('pipeline_inputs.json', 'w') as f:
    json.dump(pipeline_inputs, f, indent=2)
df_sweep.to_csv('pycalphad_outputs.csv', index=False, float_format='%.8f')
print("\n  pipeline_inputs:")
for k, val in pipeline_inputs.items():
    print(f"    {k:35s} : {val}")
print("  Saved: pipeline_inputs.json + pycalphad_outputs.csv")

# ================================================================
# SECTION 2 — JMAK Precipitation Kinetics
# ================================================================
print("\n" + "="*70)
print("  SECTION 2 — JMAK KINETICS")
print("="*70)

pi = pipeline_inputs
df_aging_loaded = pd.read_csv('aging_equilibrium.csv')

n_av   = 1.5
k0     = 1.0e6
Q_jmak = 150000
aging_times = np.logspace(1, 5, 200)
jmak_T_C = [400, 475, 500, 600]

f_eq_lookup = {}
for T_C in jmak_T_C:
    row = df_aging_loaded.loc[df_aging_loaded['T_C'] == T_C]
    if len(row):
        f_eq_lookup[T_C] = float(row['f_eq_Ni3Ti'].iloc[0])
    else:
        interp_f = interp1d(df_sweep['T_K'], df_sweep['f_NI3TI'],
                            kind='linear', fill_value=0, bounds_error=False)
        f_eq_lookup[T_C] = float(interp_f(T_C + 273.15))

print(f"\n[6/10] JMAK: n={n_av}, k0={k0:.0e}, Q={Q_jmak/1000:.0f} kJ/mol")
for T_C, feq in f_eq_lookup.items():
    T_K = T_C + 273.15
    k_T = k0 * np.exp(-Q_jmak / (R_gas * T_K))
    print(f"  {T_C} C: f_eq = {feq:.5f},  k = {k_T:.3e}")

jmak_results = {}
for T_C in jmak_T_C:
    T_K   = T_C + 273.15
    k_T   = k0 * np.exp(-Q_jmak / (R_gas * T_K))
    X_t   = 1.0 - np.exp(-k_T * aging_times**n_av)
    f_eq  = f_eq_lookup[T_C]
    f_precip = X_t * f_eq
    jmak_results[T_C] = {'times': aging_times.copy(), 'X_t': X_t,
                          'f_precip': f_precip, 'f_eq': f_eq}

# ---- PLOT 3 ----
print("\n  Generating Plot 3 ...")
fig3, ax3 = plt.subplots(figsize=(11, 6))
for T_C in jmak_T_C:
    col = COLORS.get(T_C, '#333')
    ax3.semilogx(aging_times, jmak_results[T_C]['f_precip'] * 100,
                 '-', color=col, lw=2.2, label=f'{T_C} C')
    ax3.axhline(jmak_results[T_C]['f_eq'] * 100, ls='--', color=col, alpha=0.4)
ax3.set_xlabel('Aging time (s)')
ax3.set_ylabel('Ni3Ti precipitate fraction f_p (%)')
ax3.set_title(f'Plot 3 - JMAK: Ni3Ti volume fraction vs aging time\n'
              f'(cf. Fig 8A, Ahluwalia 2024)  n={n_av}, Q={Q_jmak/1000:.0f} kJ/mol')
ax3.legend(title='Aging T', loc='lower right', framealpha=0.9)
ax3.set_ylim(bottom=-0.2)
ax3.grid(True, alpha=0.25)
ax3b = ax3.twiny()
ax3b.set_xscale('log')
ax3b.set_xlim(np.array(ax3.get_xlim()) / 3600)
ax3b.set_xlabel('Aging time (hours)')
plt.savefig('plot3_jmak_fraction.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: plot3_jmak_fraction.png")

# ================================================================
# SECTION 3 — LSW Coarsening
# ================================================================
print("\n" + "="*70)
print("  SECTION 3 — LSW COARSENING")
print("="*70)

sigma_int = 0.2
D0_lsw    = 1.0e-4
Q_D       = 150000
r0        = 2.0e-9
C_eq_Ti   = pi['X_Ti_matrix_755K']

print(f"\n[7/10] LSW: sigma={sigma_int}, D0={D0_lsw:.0e}, Q_D={Q_D/1000:.0f} kJ/mol")
print(f"  Vm={Vm:.1e}, r0={r0*1e9:.1f} nm, C_eq={C_eq_Ti:.6f}")

lsw_results = {}
for T_C in jmak_T_C:
    T_K = T_C + 273.15
    D_T   = D0_lsw * np.exp(-Q_D / (R_gas * T_K))
    K_LSW = (8.0 * sigma_int * D_T * Vm * C_eq_Ti) / (9.0 * R_gas * T_K)
    r_t   = (r0**3 + K_LSW * aging_times) ** (1.0/3.0)
    r_nm  = r_t * 1e9
    lsw_results[T_C] = {'times': aging_times.copy(), 'r_m': r_t,
                         'r_nm': r_nm, 'K_LSW': K_LSW, 'D_T': D_T}
    print(f"  {T_C} C:  D={D_T:.2e}  K_LSW={K_LSW:.2e}  "
          f"r(3h)={((r0**3+K_LSW*10800)**(1/3))*1e9:.1f} nm  "
          f"r(28h)={((r0**3+K_LSW*100000)**(1/3))*1e9:.1f} nm")

# ---- PLOT 4 ----
print("\n  Generating Plot 4 ...")
fig4, ax4 = plt.subplots(figsize=(11, 6))
for T_C in jmak_T_C:
    col = COLORS.get(T_C, '#333')
    ax4.semilogx(aging_times, lsw_results[T_C]['r_nm'], '-', color=col, lw=2.2,
                 label=f'{T_C} C  (K={lsw_results[T_C]["K_LSW"]:.2e})')
ax4.axhline(5,  ls=':', color='gray', alpha=0.5)
ax4.axhline(20, ls=':', color='gray', alpha=0.5)
ax4.fill_between(aging_times, 10, 40, color='gray', alpha=0.08,
                 label='Expt. range 10-40 nm (Jagle 2017)')
ax4.set_xlabel('Aging time (s)')
ax4.set_ylabel('Mean precipitate radius (nm)')
ax4.set_title('Plot 4 - LSW coarsening: Ni3Ti radius vs aging time\n'
              '(cf. Fig 9A, Ahluwalia 2024)')
ax4.legend(title='Aging T', loc='upper left', fontsize=9, framealpha=0.9)
ax4.grid(True, alpha=0.25)
ax4b = ax4.twiny()
ax4b.set_xscale('log')
ax4b.set_xlim(np.array(ax4.get_xlim()) / 3600)
ax4b.set_xlabel('Aging time (hours)')
plt.savefig('plot4_lsw_coarsening.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: plot4_lsw_coarsening.png")

# ================================================================
# SECTION 4 — Orowan-Ashby Yield Strength
# ================================================================
print("\n" + "="*70)
print("  SECTION 4 — OROWAN-ASHBY YIELD STRENGTH")
print("="*70)

Y_alpha = 1000e6
Y_gamma = 250e6
b       = 0.286e-9
C44     = 65.7e9
ar      = 1.0
f_gamma = 0.0

print(f"\n[8/10] Orowan-Ashby (Eq. 22):")
print(f"  Y_alpha={Y_alpha/1e6:.0f} MPa, C44={C44/1e9:.1f} GPa, b={b*1e9:.3f} nm")

yield_results = {}
for T_C in jmak_T_C:
    f_p = jmak_results[T_C]['f_precip']
    r   = lsw_results[T_C]['r_m']

    sigma_P = np.zeros_like(f_p)
    valid = (r > b) & (f_p > 1e-12)
    sigma_P[valid] = (np.sqrt(ar) * 0.269 * b * C44
                      * np.sqrt(f_p[valid]) / r[valid]
                      * np.log(r[valid] / b))
    f_alpha = 1.0 - f_p
    YT = (Y_alpha + sigma_P) * f_alpha + Y_gamma * f_gamma

    yield_results[T_C] = {'times': aging_times.copy(),
                          'sigma_P': sigma_P, 'YT': YT}

    idx_peak = np.argmax(YT)
    print(f"  {T_C} C:  Peak YT = {YT[idx_peak]/1e6:.0f} MPa  "
          f"at t = {aging_times[idx_peak]/3600:.2f} h  "
          f"(sigma_P = {sigma_P[idx_peak]/1e6:.0f} MPa,  "
          f"r = {lsw_results[T_C]['r_nm'][idx_peak]:.1f} nm,  "
          f"f_p = {f_p[idx_peak]*100:.2f} %)")

# ---- PLOT 5 + PLOT 6 ----
print("\n[9/10] Generating Plots 5 + 6 ...")
t_hours = aging_times / 3600.0

fig56, (ax5, ax6) = plt.subplots(1, 2, figsize=(16, 6))
for T_C in jmak_T_C:
    col = COLORS.get(T_C, '#333')
    ax5.semilogx(t_hours, yield_results[T_C]['YT'] / 1e6,
                 '-', color=col, lw=2.2, label=f'{T_C} C')
ax5.axhline(1000, ls='--', color='gray', alpha=0.5, label='Martensite baseline')
ax5.axhspan(1700, 1900, color='#2A9D8F', alpha=0.10, label='Expt. peak (M300)')
ax5.set_xlabel('Aging time (hours)')
ax5.set_ylabel('Yield strength Y_T (MPa)')
ax5.set_title('Plot 5 - Total yield strength\n(cf. Fig 10A, Ahluwalia 2024)')
ax5.legend(loc='best', fontsize=9, framealpha=0.9)
ax5.grid(True, alpha=0.25)
ax5.set_xlim(t_hours[0], t_hours[-1])

for T_C in jmak_T_C:
    col = COLORS.get(T_C, '#333')
    ax6.semilogx(t_hours, yield_results[T_C]['sigma_P'] / 1e6,
                 '-', color=col, lw=2.2, label=f'{T_C} C')
ax6.set_xlabel('Aging time (hours)')
ax6.set_ylabel('Precipitate strengthening sigma_P (MPa)')
ax6.set_title('Plot 6 - Orowan-Ashby strengthening\n(cf. Fig 10B, Ahluwalia 2024)')
ax6.legend(loc='best', fontsize=9, framealpha=0.9)
ax6.grid(True, alpha=0.25)
ax6.set_xlim(t_hours[0], t_hours[-1])

plt.savefig('plot5_plot6_yield_strength.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: plot5_plot6_yield_strength.png")

# ---- Summary ----
print("\n[10/10] Summary ...")
peak_strength, peak_time_h, r_at_peak, f_eq_used = {}, {}, {}, {}
for T_C in jmak_T_C:
    YT  = yield_results[T_C]['YT']
    idx = np.argmax(YT)
    peak_strength[str(T_C)] = round(float(YT[idx] / 1e6), 1)
    peak_time_h[str(T_C)]   = round(float(aging_times[idx] / 3600), 3)
    r_at_peak[str(T_C)]     = round(float(lsw_results[T_C]['r_nm'][idx]), 2)
    f_eq_used[str(T_C)]     = round(float(f_eq_lookup[T_C]), 6)

final_outputs = {
    'peak_strength_MPa': peak_strength,
    'peak_time_hours':   peak_time_h,
    'solvus_T_C':        pi['solvus_T_C'],
    'r_at_peak_nm':      r_at_peak,
    'f_eq_used':         f_eq_used
}
with open('pipeline_outputs.json', 'w') as f:
    json.dump(final_outputs, f, indent=2)

print("\n" + "="*78)
print("                 FULL PIPELINE SUMMARY")
print("="*78)
print(f"{'T (C)':>8s} | {'f_eq':>8s} | {'Peak YT (MPa)':>14s} | "
      f"{'Peak time (h)':>13s} | {'r at peak (nm)':>15s}")
print("-" * 78)
for T_C in jmak_T_C:
    k = str(T_C)
    print(f"{T_C:>8d} | {f_eq_used[k]:>8.4f} | {peak_strength[k]:>14.1f} | "
          f"{peak_time_h[k]:>13.3f} | {r_at_peak[k]:>15.2f}")
print("-" * 78)
print(f"Solvus: {pi['solvus_T_C']:.1f} C  ({pi['solvus_T_K']:.1f} K)")
print(f"f_eq at 755 K: {pi['f_eq_755K']:.4f} ({pi['f_eq_755K']*100:.2f} %)")
print(f"Driving force: {pi['driving_force_755K_J_per_mol']:.1f} J/mol")

print("\n" + "="*70)
print("  ALL OUTPUT FILES:")
for fn in ['pycalphad_outputs.csv', 'aging_equilibrium.csv',
           'pipeline_inputs.json', 'pipeline_outputs.json',
           'plot1_phase_fraction_vs_T.png', 'plot2_feq_barchart.png',
           'plot3_jmak_fraction.png', 'plot4_lsw_coarsening.png',
           'plot5_plot6_yield_strength.png']:
    exists = os.path.exists(fn)
    sz = os.path.getsize(fn) if exists else 0
    print(f"    {'[OK]' if exists else '[!!]'}  {fn:40s}  ({sz:>8,d} bytes)")
print("="*70)
print("  PIPELINE COMPLETE.")

# ================================================================
# SECTION 5 — COMPREHENSIVE DATA EXPORT
# ================================================================
# Save ALL time-series arrays and parameters so run_pipeline_final.py
# can import them directly without any hardcoded values.
print("\n" + "=" * 70)
print("  SECTION 5 — EXPORTING FULL DATA TO JSON")
print("=" * 70)

export = {
    # ---- Metadata ----
    "description": "Full pipeline output from run_pipeline.py (pycalphad + JMAK + LSW + Orowan)",

    # ---- Time array ----
    "time_seconds": aging_times.tolist(),
    "time_hours": (aging_times / 3600.0).tolist(),
    "n_points": len(aging_times),

    # ---- Temperatures ----
    "temperatures_C": jmak_T_C,

    # ---- Physical constants used ----
    "constants": {
        "Y_alpha_Pa": float(Y_alpha),
        "Y_gamma_Pa": float(Y_gamma),
        "b_m": float(b),
        "C44_Pa": float(C44),
        "ar_assumed": float(ar),
        "f_gamma_assumed": float(f_gamma),
        "Vm_m3_per_mol": float(Vm),
        "R_gas": float(R_gas),
    },

    # ---- JMAK kinetic parameters ----
    "jmak_params": {
        "n_avrami": float(n_av),
        "k0": float(k0),
        "Q_jmak_J_per_mol": float(Q_jmak),
    },

    # ---- LSW coarsening parameters ----
    "lsw_params": {
        "sigma_interface_J_per_m2": float(sigma_int),
        "D0_m2_per_s": float(D0_lsw),
        "Q_diffusion_J_per_mol": float(Q_D),
        "r0_m": float(r0),
        "C_eq_Ti": float(C_eq_Ti),
    },

    # ---- Thermodynamic inputs from pycalphad ----
    "thermodynamic": pipeline_inputs,

    # ---- Equilibrium fractions at aging temperatures ----
    "f_eq_at_aging_T": {str(T): float(f_eq_lookup[T]) for T in jmak_T_C},

    # ---- Full equilibrium FCC (reverted austenite) at aging temps ----
    # Queried with BCC+FCC+ETA — eliminates hardcoded f_gamma_sat
    "f_gamma_eq_at_aging_T": {str(T): float(f_gamma_eq_at_aging_T.get(T, 0.0)) for T in jmak_T_C},
    "f_eta_full_eq_at_aging_T": {str(T): float(f_eta_full_eq_at_aging_T.get(T, 0.0)) for T in jmak_T_C},

    # ---- Phase fraction sweep - metastable (for Plot 1) ----
    "sweep": {
        "T_K": sweep_T.tolist(),
        "T_C": (sweep_T - 273.15).tolist(),
        "f_NI3TI": sweep_fNI3TI.tolist(),
        "f_A2": sweep_fA2.tolist(),
        "solvus_T_K": float(solvus_T_K),
        "solvus_T_C": float(solvus_T_C),
    },

    # ---- Fig 2: Full equilibrium sweep (BCC+FCC+ETA+LIQUID) ----
    "fig2_full_equilibrium": {
        "T_C": T_fig2_C.tolist(),
        "T_K": T_fig2_K.tolist(),
        "f_ETA": fig2_fETA.tolist(),
        "f_FCC": fig2_fFCC.tolist(),
        "f_BCC": fig2_fBCC.tolist(),
    },

    # ---- Aging equilibrium table (for Plot 2) ----
    "aging_equilibrium": {
        "T_C": df_aging['T_C'].tolist(),
        "f_eq_Ni3Ti": df_aging['f_eq_Ni3Ti'].tolist(),
    },

    # ---- Per-temperature time-series data ----
    "kinetic_data": {},
}

for T_C in jmak_T_C:
    T_K = T_C + 273.15
    k_T = k0 * np.exp(-Q_jmak / (R_gas * T_K))
    D_T = D0_lsw * np.exp(-Q_D / (R_gas * T_K))

    export["kinetic_data"][str(T_C)] = {
        # JMAK precipitation
        "f_precip": jmak_results[T_C]['f_precip'].tolist(),
        "X_transformed": jmak_results[T_C]['X_t'].tolist(),
        "f_eq": float(jmak_results[T_C]['f_eq']),
        "k_T": float(k_T),

        # LSW coarsening
        "r_nm": lsw_results[T_C]['r_nm'].tolist(),
        "r_m": lsw_results[T_C]['r_m'].tolist(),
        "K_LSW": float(lsw_results[T_C]['K_LSW']),
        "D_T": float(lsw_results[T_C]['D_T']),

        # Orowan-Ashby yield strength
        "sigma_P_Pa": yield_results[T_C]['sigma_P'].tolist(),
        "YT_Pa": yield_results[T_C]['YT'].tolist(),
    }

export_path = 'pipeline_full_timeseries.json'
with open(export_path, 'w') as f:
    json.dump(export, f, indent=2)

sz = os.path.getsize(export_path)
print(f"\n  Saved: {export_path}  ({sz:,d} bytes)")
print(f"  Contains {len(jmak_T_C)} temperatures × {len(aging_times)} time points")
print(f"  Arrays: f_precip, r_nm, sigma_P_Pa, YT_Pa for each temperature")
print(f"  Plus: sweep data, equilibrium table, all parameters")
print("=" * 70)
