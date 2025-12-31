# FT (Fischer–Tropsch) biomass-to-jet model

This document describes the techno-economic model implemented in Biomass_to_Jet_FT_MILP.py.

This model was developed as part of an internship project at EDF (Électricité de France).
Due to confidentiality constraints related to equipment vendors, technology package providers, and project-specific operating data, the model in this repository is not the exact internship deliverable and does not contain any proprietary information. Instead, it is a sanitized and adapted version built solely from publicly available / open-source data and literature assumptions.

The purpose of releasing this code is to share the modeling workflow, decision logic, and analytical approach used during the internship (process-to-economics mapping, uncertainty handling, and investment evaluation), rather than to reproduce any confidential industrial design or performance results.

It supports two run modes:

Deterministic evaluation (no solver) — evaluates a fixed design/operation defined in model.fixed (useful for debugging and reproducing numbers based on the public dataset).

Two-stage stochastic MILP (PuLP + solver) — optimizes first-stage investment/contract decisions and second-stage scenario operations, maximizing 20-year NPV under an 8-year financing gate, with optional risk aversion (CVaR).

> **Wet-basis convention:** Biomass procurement is **as-received (wet basis)**. Drying reduces mass before gasification.

---

## 1. How to run

### 1) Deterministic evaluation (no solver required)
```bash
python Biomass_to_Jet_FT_MILP.py --data Biomass_to_Jet_FT_Data.json --out results/ft.csv --no-solver
```

### 2) MILP optimization (requires PuLP + a MILP solver backend)
```bash
python Biomass_to_Jet_FT_MILP.py --data Biomass_to_Jet_FT_Data.json --out results/ft.csv
```

### 3) Choose solver / show solver log
```bash
python Biomass_to_Jet_FT_MILP.py --data ... --out results/ft.csv --solver cbc
python Biomass_to_Jet_FT_MILP.py --data ... --out results/ft.csv --solver gurobi --solver-msg
```

### 4) Counterfactual reporting (MILP only)
```bash
python Biomass_to_Jet_FT_MILP.py --data ... --out results/ft_v2.csv \
  --report-counterfactual --report-both-upgrade
```

### 5) Risk-averse optimization (MILP only; CVaR)
```bash
python Biomass_to_Jet_FT_MILP.py --data ... --out results/ft.csv \
  --risk-alpha 0.9 --risk-lambda 0.3
```

Outputs:
- `*_policy.csv`: scenario-by-scenario rows for the “policy” run.
- (optional) `*_cf_*.csv`: counterfactual runs (forced build / relaxed gate).
- The script prints a JSON bundle to stdout containing **NPV**, **expected profit**, and **scenario rows**.

---

## 2. Model structure (high level)

Process chain:

**Biomass (wet, as-received)**  
→ **Drying / size reduction**  
→ **Entrained-flow gasification (O₂ supplied by ASU)**  
→ **Syngas (H₂, CO)**  
→ **(optional) WGS module (linearized)**  
→ **FT synthesis (no RWGS in this model)**  
→ **Liquid split** (Jet / Diesel / Naphtha / Wax)  
→ **(optional) Wax upgrading** (Wax → Jet/Diesel/Naphtha increments)

Uncertainty is represented by discrete **scenarios** `s ∈ S` with probability `π_s`. Scenario parameters include: char yield, syngas yields, ASU electricity, FT yield, product split, drying heat, etc.

---

## 3. Sets and indices

- `r ∈ R`: biomass supply rings (distance bins)
- `s ∈ S`: scenarios, with probability `π_s`
- `k ∈ K`: **discrete** jet capacity options (t/y)

---

## 4. Data (JSON) and units

### Key JSON blocks
- `feedstock.*` : ultimate analysis (dry wt%), dry LHV (MJ/kg)
- `process.*` : moisture in/out (wet basis), entrained-flow O₂ ratio, FT target H₂/CO, electricity intensities  
  - `process.wgs.enable_WGS` (optional, default false)  
  - `process.wgs.steam_cost_RMB_per_kmol_CO_shift` (optional)
- `biomass_rings[]` : `distance_km`, `max_supply_tpy` (as-received)
- `economic.*` : prices (RMB), operating hours (h/y), discount/degradation rates, fixed O&M  
  - `economic.peak_load_factor` (optional, default 1.0)
- `capex.*` : discrete capacity list + CAPEX0 per option, optional wax-upgrading CAPEX0
- `scenarios[]` : scenario-specific yields / utilities / splits / probabilities
- `model.bigM_wax` : Big-M for wax-upgrading activation
- `model.fixed` : fixed choices for deterministic evaluation
- (optional) `model.min_biomass_if_built_tpy` : prevents degenerate “build=1 but no feed”
- (optional) `model.wax_upgrading_max_tpy` : upgrading capacity upper bound

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
- `cap_choice_k ∈ {0,1}` : choose one capacity option `k`
- `cap_jet_tpy ≥ 0` : nameplate jet capacity implied by `cap_choice_k`
- `wax_upgrade ∈ {0,1}` : build wax upgrading or not
- `bio_ring_r_tpy ≥ 0` : biomass procured from ring `r` (as-received)
- `P_peak_kW ≥ 0` : contracted peak power for demand charge (covers all scenarios)

Discrete CAPEX0:
- `capex0_RMB = Σ_k capex0_k * cap_choice_k + capex_up * wax_upgrade`

### 5.2 Second-stage (operations; per scenario `s`)
Core variables per scenario:
- `dry_solids_tpy_s`, `water_removed_tpy_s`
- `char_dry_tpy_s`, `char_daf_tpy_s`
- `H2_Nm3_y_s`, `CO_Nm3_y_s` and `H2_kmol_y_s`, `CO_kmol_y_s`
- (optional) `CO_shift_kmol_y_s` (WGS)
- `CO_used_kmol_y_s`
- `H2_makeup_tpy_s`
- `liquids_tpy_s`
- Base products: `jet_base_tpy_s`, `diesel_base_tpy_s`, `naphtha_base_tpy_s`, `wax_tpy_s`
- Wax upgrading: `wax_upgraded_tpy_s`, `wax_sold_tpy_s`
- Final products: `jet_total_tpy_s`, `diesel_total_tpy_s`, `naphtha_total_tpy_s`
- Electricity: `E_kWh_s`

---

## 6. Core constraints (what the code enforces)

### 6.1 Build & discrete capacity
- Choose exactly one capacity option if built: `Σ_k cap_choice_k = build`
- Upgrading only if built: `wax_upgrade ≤ build`
- Ring supplies: `bio_ring_r_tpy ≤ max_supply_r * build`
- Total biomass: `B = Σ_r bio_ring_r_tpy`
- (optional) minimum feed if built: `B ≥ min_biomass_if_built_tpy * build`

### 6.2 Moisture balance (wet → dry)
Let `x_in` be inlet moisture (wet basis) and `x_out` outlet moisture.
- `dry_solids = B * (1 − x_in)`
- `after_dry_ar = dry_solids / (1 − x_out)`
- `water_removed = B − after_dry_ar`

### 6.3 Char, ash correction, syngas yields
- `char_dry = y_char_s * dry_solids`
- `char_daf = char_dry − ash_frac * dry_solids`, with `char_daf ≥ 0`
- `H2_Nm3 = Y_H2_s * char_daf * 1000`
- `CO_Nm3 = Y_CO_s * char_daf * 1000`
- Convert to kmol:
  - `H2_kmol = H2_Nm3 / Nm3_per_kmol`
  - `CO_kmol = CO_Nm3 / Nm3_per_kmol`

### 6.4 (Optional) linear WGS
If enabled:
- `0 ≤ CO_shift ≤ CO_kmol`
- Effective syngas:
  - `H2_eff = H2_kmol + CO_shift`
  - `CO_eff = CO_kmol − CO_shift`

If not enabled:
- `H2_eff = H2_kmol`, `CO_eff = CO_kmol`

### 6.5 FT conversion and H₂ make-up (no RWGS)
- `CO_used ≤ CO_eff`
- `H2_makeup_kmol ≥ R_FT * CO_used − H2_eff`, `H2_makeup_kmol ≥ 0`
- `H2_makeup_tpy = H2_makeup_kmol * 2 / 1000`
- `liquids = η_liq_s * CO_used`

### 6.6 Product split + **final jet capacity**
Base split from FT liquids:
- `jet_base = f_ker_s * liquids`
- `diesel_base = f_dis_s * liquids`
- `naphtha_base = f_nap_s * liquids`
- `wax = f_wax_s * liquids`

Wax upgrading (Big-M activation):
- `wax_upgraded ≤ wax`
- `wax_upgraded ≤ M_wax * wax_upgrade`
- `wax_sold = wax − wax_upgraded`
- (optional) upgrading capacity: `wax_upgraded ≤ wax_upgrading_max_tpy * wax_upgrade`

Final products:
- `jet_total = jet_base + γ_jet_s * wax_upgraded`
- `diesel_total = diesel_base + γ_dis_s * wax_upgraded`
- `naphtha_total = naphtha_base + γ_nap_s * wax_upgraded`

**Jet capacity constraint (nameplate):**
- `jet_total ≤ cap_jet_tpy`

> In this model, the jet capacity is enforced on the **final jet stream** including wax-upgrading increments.

### 6.7 Oxygen requirement + electricity peak sizing
O₂ requirement:
- `O2_required = ρ_O2 * dry_solids`

Electricity components (linear):
- Grinding/feed press: proportional to thermal input from `LHV_dry`
- ASU: `ASU_kWh_per_tO2_s * O2_required`
- FT: `kWh_per_t_liquids * liquids`
- Upgrading: `kWh_per_t_wax * wax_upgraded`

Peak sizing (first-stage) uses a linear “load factor” approximation:
- `E_kWh_s ≤ (H_op * peak_load_factor) * P_peak_kW` for all scenarios `s`
  - `peak_load_factor = 1.0` reduces to the average-power bound.

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
- `C_H2_s = p_H2 * H2_makeup_tpy_s`
- `C_heat_s = p_heat * drying_heat_MJ_s`
- `C_FOM = fixed_OM * build`
- (optional) WGS utility cost:
  - `C_WGS_s = steam_cost_RMB_per_kmol_CO_shift * CO_shift_kmol_y_s`

### 7.4 Revenue
- `Rev_s = P_jet * jet_total + P_diesel * diesel_total + P_naphtha * naphtha_total + P_wax * wax_sold`

Profit per scenario:
- `Profit_s = Rev_s − (C_bio + C_elec_s + C_H2_s + C_heat_s + C_FOM + C_WGS_s)`

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
- `NPV_8  = −Capex0 + PVF(8)  * E[Profit]`
- `NPV_20 = −Capex0 + PVF(20) * E[Profit]`

Policy financing gate:
- `NPV_8 ≥ 0`
If no configuration satisfies the gate, the optimal solution becomes `build = 0` (no build).

---

## 9. Risk-averse option (CVaR of annual profit)

If enabled via `--risk-alpha a --risk-lambda λ`, the model adds a downside-risk penalty based on the **CVaR** of annual profit.

Let `z` be VaR and `η_s` be scenario shortfall variables:
- `η_s ≥ z − Profit_s`, `η_s ≥ 0`
- `CVaR = z − 1/(1−a) * Σ_s π_s * η_s`

Risk penalty:
- `Penalty = E[Profit] − CVaR` (higher means worse downside risk)

Risk-adjusted objective (conceptually):
- maximize `NPV_20 − λ * PVF(20) * Penalty`

Setting `λ=0` recovers the risk-neutral objective.

---

## 10. Solver

The MILP is built with **PuLP** and solved by a backend solver:
- Default: `CBC` (`--solver cbc`)
- Optional: `HiGHS` (`--solver highs`, if supported by your PuLP install)
- Optional: `Gurobi` (`--solver gurobi`, if installed and licensed)

---

## 11. Deterministic evaluation mode (debugging)

If you run with `--no-solver` (or PuLP is not installed), the script evaluates `model.fixed`:
- Biomass allocation: greedy ring allocation from nearest rings until `biomass_in_ar_tpy` is met.
- CO usage is set to the maximum allowed by **final jet capacity** and available syngas.
- If `wax_upgrade=1`, wax is upgraded according to the deterministic rule; otherwise wax is sold.
- The same techno-economic calculations are performed per scenario and aggregated into NPV.

This mode is intended for reproducible calculations and sanity checks; it does **not** optimize decisions.
