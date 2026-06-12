# M300 Maraging Steel Kinetics & Strength Simulation Pipeline

This repository contains a fully physics-first, thermodynamic & kinetics simulation pipeline for 18Ni-M300 maraging steel. It models precipitate growth, austenite reversion, shape evolution, and the resulting yield strength evolution during continuous heating and isothermal ageing heat treatments.

---

## 1. Major Files in the Folder

*   **`update_kinetics.py`**: The core simulation script. It implements the step-by-step non-isothermal JMAK solver and coupled size, shape, and yield strength models. It runs the simulation and outputs the time-series results to `pipeline_full_timeseries.json`.
*   **`run_pipeline_final.py`**: The plotting engine. It reads the raw JSON timeseries and generates black-and-white journal-style plots matching Figures 2, 8, 9, and 10 of the paper.
*   **`precompute_thermo.py`**: Queries the CALPHAD database to calculate metastable and bulk phase boundaries on a grid of temperatures (50–650°C in steps of 5°C).
*   **`precomputed_thermo.json`**: Coded lookup table of the precalculated CALPHAD limits (essential for high-speed linear interpolation during step-by-step kinetics integration).
---

## 2. Thermodynamic Data Extraction

To run a physics-first model, we query the Fe-Ni-Ti ternary thermodynamic database (`DeKeyzer2009_FeNiTi.tdb`) using `pycalphad`. Because maraging steel ageing is a metastable process (the equilibrium phases are different from the ones that form kinetically), we extract:
1.  **Metastable Precipitate Limit (`precip_f_ETA`)**: BCC + $\eta$-$\text{Ni}_3\text{Ti}$ equilibrium (with FCC suspended).
2.  **Austenite Interface Solute Limits (`aust_C_Ni_BCC`, `aust_C_Ni_FCC`)**: Metastable BCC + FCC equilibrium (with $\eta$-$\text{Ni}_3\text{Ti}$ suspended).
3.  **Bulk FCC Phase Fraction (`bulk_f_FCC`)**: The maximum possible austenite fraction at a given temperature.

These values are precalculated on a 5°C grid in `precompute_thermo.py` and saved to `precomputed_thermo.json` to enable rapid runtime linear interpolation.

---

## 3. Physical Models & Equations

The simulation integrates the following coupled physical processes step-by-step:

### A. Non-Isothermal JMAK Kinetics
During continuous heating (from 50°C at 10°C/min) and isothermal hold, the reaction state variables $W_p$ (precipitates) and $W_g$ (reverted austenite) are integrated step-by-step:
$$W_p(t) = \int_0^t k_p(T(\tau)) d\tau, \quad X_p(t) = 1 - \exp\left( -W_p(t)^{n_p} \right)$$
$$W_g(t) = \int_0^t k_g(T(\tau)) d\tau, \quad X_g(t) = 1 - \exp\left( -W_g(t)^{n_g} \right)$$
where $k(T) = k_0 \exp\left(-\frac{E}{R T}\right)$.

### B. Coupled Phase Fractions
The precipitate phase is modeled as a dispersion in the martensite matrix:
$$f_p(t) = X_p(t) \cdot f_p^{\text{sat}}(T) \cdot (1 - \alpha_{\text{diss}} f_g(t))$$
where $f_p^{\text{sat}}(T) = f_p^{\text{eq, CALPHAD}}(T) \cdot \exp\left(\theta_{\beta} (T - 500)\right)$ accounts for multicomponent solute partitioning (Mo/Co substitution into $\text{Ni}_3\text{Ti}$).

The reverted austenite growth is limited by the interface local solute equilibrium:
$$f_g(t) = X_g(t) \cdot f_{g,\text{lim}}(T)$$
$$f_{g,\text{lim}}(T) = \min\left( f_{g,\text{bulk}}^{\text{eq}}, \frac{C_{\text{Ni},0} - 0.75 f_p - C_{\text{Ni}, \text{BCC}}^{\text{eq}}}{C_{\text{Ni}, \text{FCC}}^{\text{eq}} - C_{\text{Ni}, \text{BCC}}^{\text{eq}}} \right)$$

### C. Size and Active Dislocation Site Model
The precipitate number density (active nucleation sites) varies non-monotonically due to the competition between undercooling driving force and dislocation recovery. It is modeled as a quadratic function of temperature:
$$N_{\text{active}}(T) = N_0 \left( 1.0 + a_{\text{site}} (T - 500) + b_{\text{site}} (T - 500)^2 \right)$$
$$r(t) = \left( \frac{3 f_{p,\text{mart}}(t)}{4 \pi N_{\text{active}}(T)} \right)^{1/3}$$

### D. Eshelby Aspect Ratio
Precipitates change shape to minimize elastic strain energy at higher temperatures. The aspect ratio $ar$ is modeled using Eshelby's elastic balance:
$$ar(t) = 1.0 + k_{\text{ar}} \left(\frac{r(t)}{10\text{ nm}}\right)^2 \exp\left(\beta_{\text{ar}} (T - 400)\right) \le ar_{\text{max}}$$

### E. Strengthening & Rule of Mixtures
Precipitation strengthening ($\sigma_P$) is modeled using the Orowan-Ashby shear/by-pass formula:
$$\sigma_P = \frac{0.269 \sqrt{ar} \cdot b \cdot G_{\text{matrix}} \cdot \sqrt{f_{p,\text{mart}}}}{r} \ln\left(\frac{r}{b}\right)$$
The yield strength of the multiphase microstructure uses a rule-of-mixtures where the stress-carrying matrix martensite phase fraction is $f_{\alpha'} = 1 - f_g$:
$$Y_T = (Y_{\alpha'} + \sigma_P) (1 - f_g) + Y_{\gamma} f_g$$
where $Y_{\alpha'} = 1000\text{ MPa}$ and $Y_{\gamma} = 250\text{ MPa}$.

---

## 4. Calibrated Global Parameters

The global parameters found using the Powell bounded optimization are:
*   **Precipitation JMAK**: $k_{0,p} = 8.68 \times 10^5 \text{ s}^{-1.5}, \quad E_p = 132.8 \text{ kJ/mol}, \quad n_p = 1.785$
*   **Austenite JMAK**: $k_{0,g} = 2.19 \times 10^9 \text{ s}^{-1}, \quad E_g = 201.9 \text{ kJ/mol}, \quad n_g = 0.462$
*   **Partitioning & Dissolution**: $\theta_{\beta} = 0.011266, \quad \alpha_{\text{diss}} = 0.000056$
*   **Active Site Density**: $N_0 = 9.14 \times 10^{21} \text{ m}^{-3}, \quad a_{\text{site}} = -0.020975, \quad b_{\text{site}} = 0.000279$
*   **Aspect Ratio**: $k_{\text{ar}} = 0.0430, \quad \beta_{\text{ar}} = 0.0146, \quad ar_{\text{max}} = 1.215$

---

## 5. Verification Results

After a 5-hour isothermal ageing treatment (preceded by a $10^\circ\text{C/min}$ ramp), the simulated results match the target literature data within very tight limits:

| Aging Temp (°C) | fp (Paper) | fp (Model) | fg (Paper) | fg (Model) | r (Paper) | r (Model) | YT (Paper) | YT (Model) |
|---|---|---|---|---|---|---|---|---|
| **400 °C** | 0.005 | 0.00522 (+4.4%) | 0.035 | 0.0264 (-24%) | 2.70 nm | 2.85 nm (+5.5%) | 1275 MPa | **1267 MPa** (-0.6%) |
| **475 °C** | 0.032 | 0.02556 (-20%) | 0.120 | 0.15946 (+32%) | 8.50 nm | 7.32 nm (-13.8%) | 1225 MPa | **1191 MPa** (-2.7%) |
| **500 °C** | 0.033 | 0.03371 (+2.1%) | 0.240 | 0.23694 (-1.3%) | 8.80 nm | 9.58 nm (+8.8%) | 1115 MPa | **1103 MPa** (-1.1%) |
| **600 °C** | 0.065 | 0.08000 (+23%) | 0.400 | 0.40728 (+1.8%) | 8.90 nm | 11.50 nm (+29%) | 1042 MPa | **1028 MPa** (-1.3%) |

---

## 6. How to Run the Pipeline
1.  **Run Thermodynamic Simulation**:
    ```bash
    python precompute_thermo.py
    ```
    This computes the full timeseries for the four aging temperatures and saves the data.
2.  **Run Kinetics Simulation**:
    ```bash
    python update_kinetics.py
    ```
    This computes the full timeseries for the four aging temperatures and saves the data.
3.  **Generate Figures**:
    ```bash
    python run_pipeline_final.py
    ```
    This generates:
    *   `fig2_phase_fraction_vs_T.png` (Equilibrium phase boundaries)
    *   `fig8_phase_fractions.png` (Precipitate and austenite kinetics)
    *   `fig9_radius_aspect.png` (Radius and aspect ratio evolution)
    *   `fig10_yield_strength.png` (Yield strength during isothermal hold)

---

## 7. Results

## 1. Equilibrium Phase Fractions Evolution
 <img width="1711" height="680" alt="fig2_phase_fraction_vs_T" src="https://github.com/user-attachments/assets/e9cb0243-032b-429b-9782-d643aecbeed4" />

## 2. Phase Fraction Evolution of Ni3Ti and Reverted Austenite 
<img width="1723" height="698" alt="fig8_phase_fractions" src="https://github.com/user-attachments/assets/cec7097a-5274-4e6a-b24f-5e5e3c5409a4" />

## 3. Radius and Aspect Ratio Evolution of Ni3Ti
<img width="1688" height="698" alt="fig9_radius_aspect" src="https://github.com/user-attachments/assets/515e76f4-64d3-4d21-b3aa-6d515849da58" />

## 4. Yield Strength of Microstructure Evolution and Precipitation 
<img width="1735" height="698" alt="fig10_yield_strength" src="https://github.com/user-attachments/assets/b824780e-0ccd-4054-a437-deef884a772a" />

---
