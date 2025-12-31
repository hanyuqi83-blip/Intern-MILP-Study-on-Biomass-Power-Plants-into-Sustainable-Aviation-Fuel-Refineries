# MTJ (Methanol-to-Jet) biomass-to-jet model

This document describes the techno-economic model implemented in `Biomass_to_Jet_MTJ_MILP.py`.

This model was developed as part of an internship project at EDF (Électricité de France).
Due to confidentiality constraints related to equipment vendors, technology package providers, and project-specific operating data, the model in this repository is not the exact internship deliverable and does not contain any proprietary information. Instead, it is a sanitized and adapted version built solely from publicly available / open-source data and literature assumptions.

The purpose of releasing this code is to share the modeling workflow, decision logic, and analytical approach used during the internship (process-to-economics mapping, uncertainty handling, and investment evaluation), rather than to reproduce any confidential industrial design or performance results.

It supports:

- **Two-stage stochastic MILP (PuLP + solver)** — optimizes first-stage investment decisions and second-stage scenario operations, maximizing 20-year NPV under an 8-year financing gate, with optional risk aversion (CVaR).
- **Syngas-limited methanol** — methanol production is physically capped by CO/H₂ availability from gasification (no RWGS; WGS optional).

> **Wet-basis convention:** Biomass procurement is **as-received (wet basis)**. Drying reduces mass before gasification.

---

## 1. How to run

### 1) MILP optimization (requires PuLP + a MILP solver backend)
```bash
python Biomass_to_Jet_MTJ_MILP.py --data Biomass_to_Jet_MTJ_Data.json --out results/mtj.csv
```

### 2) Enforce minimum loading (avoid degenerate “build=1 but tiny biomass”)
```bash
python Biomass_to_Jet_MTJ_MILP.py --data Biomass_to_Jet_MTJ_Data.json --out results/mtj.csv   --min-load-frac 0.8
```

### 3) Counterfactual reporting
```bash
python Biomass_to_Jet_MTJ_MILP.py --data Biomass_to_Jet_MTJ_Data.json --out results/mtj.csv   --report-counterfactual --report-both-mtj
```

### 4) Risk-averse optimization (CVaR)
```bash
python Biomass_to_Jet_MTJ_MILP.py --data Biomass_to_Jet_MTJ_Data.json --out results/mtj.csv   --risk-alpha 0.9 --risk-lambda 0.3
```

Outputs:
- `*_policy.csv`: scenario-by-scenario rows for the “policy” run.
- (optional) `*_cf_*.csv`: counterfactual runs (forced build / relaxed gate / MTJ forced on/off).
- The script prints a JSON bundle to stdout containing **NPV**, **expected profit**, and **scenario rows**.

---

## 2. Model structure (high level)

Process chain:

**Biomass (wet, as-received)**  
→ **Drying / size reduction**  
→ **Entrained-flow gasification (O₂ supplied by ASU)**  
→ **Syngas (H₂, CO)**  
→ **(optional) WGS module (linearized)**  
→ **Methanol synthesis** *(syngas-limited; no RWGS in this model)*  
→ **MTJ conversion** *(MeOH → Jet + gasoline + naphtha)*  
→ **Products** *(sell MeOH if MTJ not built)*

Uncertainty is represented by discrete **scenarios** `s ∈ S` with probability `π_s`. Scenario parameters include: char yield, syngas yields, ASU electricity, drying heat, methanol electricity intensity, etc.

---

## 3. Sets and indices

- `r ∈ R`: biomass supply rings (distance bins)
- `s ∈ S`: scenarios, with probability `π_s`
- `k ∈ K`: **discrete** methanol capacity options (t/y)

---

## 4. Data (JSON) and units

### Key JSON blocks
- `feedstock.*` : ultimate analysis (dry wt%), dry LHV (MJ/kg)
- `process.*` :
  - moisture in/out (wet basis)
  - entrained-flow O₂ ratio
  - MTJ yields + MTJ H₂ demand + MTJ electricity
  - optional WGS switch
  - syngas conversion constant `Nm3_per_kmol`
- `biomass_rings[]` : `distance_km`, `max_supply_tpy` (as-received)
- `economic.*` : prices (RMB), operating hours (h/y), discount/degradation rates, fixed O&M  
  - includes methanol + jet (+ gasoline/naphtha if used)
- `capex.*` : discrete methanol capacity list + CAPEX0 per option, optional MTJ CAPEX0
- `scenarios[]` : scenario-specific yields / utilities / probabilities

### Units
- Mass flow: **t/y**
- Syngas: **Nm³/y**, converted to **kmol/y** using `Nm3_per_kmol`
- Electricity: **kWh/y**
- Heat: **MJ/y**
- Prices: **RMB** (per relevant unit)

---

## 5. Decision variables

### 5.1 First-stage (investment / sizing; scenario-independent)
- `build ∈ {0,1}` : build plant or not
- `cap_choice_k ∈ {0,1}` : choose one methanol capacity option `k`
- `cap_meoh_tpy ≥ 0` : methanol nameplate capacity implied by `cap_choice_k`
- `mtj_build ∈ {0,1}` : build MTJ or not (`mtj_build ≤ build`)
- `bio_ring_r_tpy ≥ 0` : biomass procured from ring `r` (as-received)
- `P_peak_kW ≥ 0` : contracted peak power for demand charge (covers all scenarios)

Discrete CAPEX0:
- `capex0_RMB = Σ_k capex0_k * cap_choice_k + capex_mtj * mtj_build`

### 5.2 Second-stage (operations; per scenario `s`)
Core variables per scenario:
- `dry_solids_tpy_s`, `water_removed_tpy_s`
- `O2_required_tpy_s`
- `syngas_H2_Nm3_y_s`, `syngas_CO_Nm3_y_s`
- `H2_kmol_y_s`, `CO_kmol_y_s`
- (optional) `CO_shift_kmol_y_s` (WGS)
- `CO_used_kmol_y_s`
- Methanol:
  - `methanol_tpy_s`
  - `methanol_to_mtj_tpy_s`
  - `methanol_sold_tpy_s`
- Products:
  - `jet_tpy_s`, `gasoline_tpy_s`, `naphtha_tpy_s`
- Electricity:
  - `electricity_kWh_y_s`
- Hydrogen:
  - `h2_makeup_tpy_s` *(for methanol synthesis)*
  - `h2_mtj_tpy_s` *(for MTJ upgrading)*
  - `h2_total_tpy_s`

---

## 6. Core constraints (what the code enforces)

### 6.1 Build & discrete capacity
- Choose exactly one capacity option if built: `Σ_k cap_choice_k = build`
- MTJ only if built: `mtj_build ≤ build`
- Ring supplies: `bio_ring_r_tpy ≤ max_supply_r * build`
- Total biomass: `B = Σ_r bio_ring_r_tpy`
- (optional) minimum loading:
  - `B ≥ min_load_frac * Bmax * build`

### 6.2 Moisture balance (wet → dry)
Let `x_in` be inlet moisture (wet basis) and `x_out` outlet moisture.
- `dry_solids = B * (1 − x_in)`
- `after_dry_ar = dry_solids / (1 − x_out)`
- `water_removed = B − after_dry_ar`

### 6.3 Oxygen requirement (ASU feed)
- `O2_required = ρ_O2 * dry_solids`

### 6.4 Char, ash correction, syngas yields (syngas branch)
Let `y_char_s` be char yield on dry solids; `ash_frac` from ultimate analysis.
- `char_dry = y_char_s * dry_solids`
- `char_daf = char_dry − ash_frac * dry_solids`, with `char_daf ≥ 0`
- `H2_Nm3 = Y_H2_s * char_daf * 1000`
- `CO_Nm3 = Y_CO_s * char_daf * 1000`
- Convert to kmol:
  - `H2_kmol = H2_Nm3 / Nm3_per_kmol`
  - `CO_kmol = CO_Nm3 / Nm3_per_kmol`

### 6.5 (Optional) linear WGS
If enabled:
- `0 ≤ CO_shift ≤ CO_kmol`
- Effective syngas:
  - `H2_eff = H2_kmol + CO_shift`
  - `CO_eff = CO_kmol − CO_shift`

If not enabled:
- `H2_eff = H2_kmol`, `CO_eff = CO_kmol`

### 6.6 Methanol synthesis and H₂ make-up (no RWGS)
Stoichiometry: `CO + 2H₂ → CH₃OH` (`R_meoh = 2`).

- `CO_used ≤ CO_eff`
- `H2_makeup_kmol ≥ R_meoh * CO_used − H2_eff`, `H2_makeup_kmol ≥ 0`
- `H2_makeup_tpy = H2_makeup_kmol * 2 / 1000`

Methanol from CO usage (with optional CO-to-MeOH efficiency `η_CO_s`):
- `meoh_raw = η_CO_s * CO_used * (32 kg/kmol) / 1000`

Capacity and feasibility:
- `methanol ≤ meoh_raw`
- `methanol ≤ cap_meoh_tpy`

> **Key idea:** methanol is “**capped by syngas**” (CO/H₂), not assumed as an independent yield.

### 6.7 Methanol allocation and MTJ activation
- `methanol_to_mtj ≤ methanol`
- MTJ only if built: `methanol_to_mtj ≤ M * mtj_build`  
  (Big-M tightened to the maximum methanol capacity option)
- `methanol_sold = methanol − methanol_to_mtj`

### 6.8 MTJ product yields
From `process.mtj.yields_t_per_t_meoh`:
- `jet = y_jet * methanol_to_mtj`
- `gasoline = y_gas * methanol_to_mtj`
- `naphtha = y_nap * methanol_to_mtj`

### 6.9 Electricity and peak contract sizing
Electricity components (linear):
- MeOH block: `e_meoh_s * methanol`
- ASU: `e_asu_s * O2_required`
- MTJ: `e_mtj * methanol_to_mtj`

Total:
- `E_kWh_s = e_meoh_s*methanol + e_asu_s*O2_required + e_mtj*methanol_to_mtj`

Peak sizing (first-stage) uses a load factor approximation:
- `E_kWh_s ≤ (H_op * peak_load_factor) * P_peak_kW` for all scenarios `s`

### 6.10 Hydrogen usage in MTJ
- `h2_mtj_tpy_s = (H2_kg_per_t_meoh_to_mtj * methanol_to_mtj) / 1000`
- `h2_total_tpy_s = h2_makeup_tpy_s + h2_mtj_tpy_s`

---

## 7. Economics

### 7.1 Biomass cost (ring-based)
For each ring:
- `cost_r = (bailing + transport * distance_km_r) * bio_ring_r_tpy`
Total:
- `C_bio = Σ_r cost_r`

### 7.2 Electricity tariff (two-part)
- `C_elec_s = p_kWh * E_kWh_s + 12 * p_demand * P_peak_kW`

### 7.3 Other OPEX
- `C_H2_s = p_H2 * h2_total_tpy_s`
- Drying heat:
  - `drying_heat_MJ_s = q_dry_s * water_removed_s * 1000`
  - `C_heat_s = p_heat * drying_heat_MJ_s`
- Fixed O&M:
  - `C_FOM = fixed_OM * build + mtj_fixed_OM * mtj_build`
- (optional) WGS utility:
  - `C_WGS_s = steam_cost_RMB_per_kmol_CO_shift * CO_shift_kmol_y_s`

Total:
- `OPEX_s = C_bio + C_elec_s + C_H2_s + C_heat_s + C_FOM + C_WGS_s`

### 7.4 Revenue
- `Rev_s = P_meoh * methanol_sold + P_jet * jet + P_gasoline * gasoline + P_naphtha * naphtha`

Profit per scenario:
- `Profit_s = Rev_s − OPEX_s`

Expected annual profit:
- `E[Profit] = Σ_s π_s * Profit_s`

---

## 8. Finance: PVF and NPV

Discount rate `i`, degradation rate `δ`.

Present value factor:
```text
PVF(N) = Σ_{t=1..N} (1/(1+i)^t) * (1-δ)^(t-1)
```

NPV definitions:
- `NPV_gate = −Capex0 + PVF(gate_years) * E[Profit]`
- `NPV_life = −Capex0 + PVF(lifetime_years) * E[Profit]` (risk-adjusted if CVaR enabled)

Policy financing gate:
- `NPV_gate ≥ 0`
If no configuration satisfies the gate, the optimal solution becomes `build = 0` (no build).

---

## 9. Risk-averse option (CVaR of annual profit)

If enabled via `--risk-alpha a --risk-lambda λ`, the model adds a downside-risk penalty based on the **CVaR** of annual profit.

Let `z` be VaR and `η_s` be scenario shortfall variables:
- `η_s ≥ z − Profit_s`, `η_s ≥ 0`
- `CVaR = z − 1/(1−a) * Σ_s π_s * η_s`

Risk-adjusted objective:
- maximize `NPV_life` with a penalty proportional to `E[Profit] − CVaR`.

Setting `λ=0` recovers the risk-neutral objective.

---

## 10. Solver

The MILP is built with **PuLP** and solved by a backend solver:
- Default: `CBC` (via `PuLP_CBC_CMD`)

(If desired, the script can be extended to support `--solver gurobi/highs` similarly to the FT model.)

---

## 11. Notes on “same gasification scenarios as FT”

If you want **FT and MTJ to share the same gasification uncertainty**, keep the same scenario fields controlling gasification:
- `char_yield_dry_fraction`
- `H2_yield_Nm3_per_kg_char_daf`
- `CO_yield_Nm3_per_kg_char_daf`
- (optional) `ASU_kWh_per_t_O2`, drying heat, etc.

Then:
- FT: syngas limits FT liquids and product split.
- MTJ: syngas limits methanol via CO/H₂ balance, then MTJ converts methanol into jet/gasoline/naphtha.

This makes the *upstream* uncertainty consistent while allowing *downstream* conversion differences.
