# FT (Fischer–Tropsch) MILP Formulation (Basic)

This document describes the **basic** MILP used in `src/ft/biomass_to_jet_ft_milp.py`.
It is a two-stage stochastic MILP with:
- **First-stage** investment decisions (build / capacity / wax upgrading / biomass sourcing + grid connection sizing),
- **Second-stage** scenario-dependent operations (syngas yields, FT conversion, utilities, revenue, OPEX).

> Notation follows the implementation: biomass is **as-received (wet basis)**, drying reduces mass.

---

## Sets and indices

- $r \in \mathcal{R}$ : biomass supply rings (distance bins)
- $s \in \mathcal{S}$ : scenarios, with probability $\pi_s$
- $k \in \mathcal{K}$ : discrete plant capacity options (Jet capacity)

---

## Parameters (examples; see JSON)

### Biomass logistics
- $\overline{B}_r$ : max biomass supply in ring $r$ (t/y, as-received)
- $d_r$ : distance of ring $r$ (km)
- $c^{\text{bale}}$ : base biomass cost (RMB/t)
- $c^{\text{tr}}$ : transport cost (RMB/(t·km))

### Moisture & gasification (wet → dry → syngas)
- $x_{\text{in}}$ : inlet moisture (wet basis)
- $x_{\text{out}}$ : outlet moisture after drying (wet basis)
- $a$ : ash fraction in dry biomass
- $y^{\text{char}}_s$ : char yield on dry basis (scenario)
- $Y^{H2}_s,\;Y^{CO}_s$ : syngas yields (Nm³ per kg char_daf)
- $\nu$ : Nm³ per kmol (constant)
- $\rho_{O2}$ : O2-to-dry-biomass mass ratio (t O2 per t dry biomass)

### FT conversion & products
- $R^{FT}$ : target $H_2/CO$ at FT inlet (no WGS)
- $\eta^{liq}_s$ : FT liquids yield (t liquids per kmol CO used)
- $(f^{jet}_s, f^{dies}_s, f^{nap}_s, f^{wax}_s)$ : liquid split fractions
- Wax upgrading yields (per t wax upgraded): $\gamma^{jet}_s,\gamma^{dies}_s,\gamma^{nap}_s$

### Utilities & economics
- Electricity energy price $p^{e}$ (RMB/kWh)
- Electricity demand charge $p^{d}$ (RMB/(kW·month))
- Operating hours $H$ (h/y)
- Heat price $p^{h}$ (RMB/MJ)
- Hydrogen price $p^{H2}$ (RMB/t)
- Fixed O&M $C^{FOM}$ (RMB/y if built)
- Product prices: $P^{jet},P^{dies},P^{nap},P^{wax}$ (RMB/t)

### Finance
- Discount rate $i$, degradation rate $\delta$
- Present value factor (used in code):
$$
PVF(N)=\sum_{t=1}^{N}\frac{(1-\delta)^{t-1}}{(1+i)^t}
$$
- Gate horizon $N_g=8$, lifetime $N_L=20$

### CAPEX
- Discrete CAPEX per option $Capex_k$
- Wax-upgrading CAPEX $Capex^{up}$

---

## Decision variables

### First-stage (investment)
- $y \in \{0,1\}$ : build plant or not
- $x_k \in \{0,1\}$ : choose capacity option $k$
- $u \in \{0,1\}$ : build wax upgrading or not
- $B_r \ge 0$ : biomass procured from ring $r$ (t/y, as-received)
- $P^{peak} \ge 0$ : contracted peak power (kW) for demand charge

Derived:
- Total biomass $B=\sum_r B_r$
- Jet capacity (t/y): $Cap^{jet}=\sum_k Cap^{jet}_k x_k$
- Total CAPEX:
$$
Capex0=\sum_k Capex_k x_k + Capex^{up}u
$$

### Second-stage (operations, per scenario $s$)
Key mass/energy and product variables:
- $B^{dry}_s,\;B^{after}_s,\;W^{rem}_s$ : dry solids, after-drying mass, removed water
- $Char^{dry}_s,\;Char^{daf}_s$
- $H2^{Nm3}_s,\;CO^{Nm3}_s,\;H2^{kmol}_s,\;CO^{kmol}_s$
- $CO^{use}_s$ : CO used in FT (kmol/y)
- $H2^{make}_s$ : H2 make-up (t/y)
- $L^{liq}_s$ : FT liquids (t/y)
- Products: $Jet_s,\;Dies_s,\;Nap_s,\;Wax_s$
- Wax upgrading: $Wax^{up}_s,\;Wax^{sell}_s$
- Electricity $E_s$ (kWh/y), drying heat $Q^{dry}_s$ (MJ/y)

---

## Constraints

### (1) Capacity selection and build logic
Choose exactly one capacity option if built:
$$
\sum_{k} x_k = y
$$
Wax upgrading only if built:
$$
u \le y
$$
Biomass only if built:
$$
0 \le B_r \le \overline{B}_r \, y \quad \forall r
$$

> Implementation note: **avoid** having a “0-capacity” option inside $\mathcal{K}$ (or enforce $Cap^{jet}\ge Cap^{min}y$),
> otherwise the solver can pick $y=1$ with zero capacity.

---

### (2) Moisture balance (as-received → dry solids)
$$
B^{dry}_s = (1-x_{in})B
$$
$$
B^{after}_s = \frac{B^{dry}_s}{1-x_{out}},\qquad
W^{rem}_s = B - B^{after}_s
$$

---

### (3) Char and syngas yields
$$
Char^{dry}_s = y^{char}_s B^{dry}_s
$$
$$
Char^{daf}_s = Char^{dry}_s - a B^{dry}_s,\qquad Char^{daf}_s \ge 0
$$
$$
H2^{Nm3}_s = Y^{H2}_s \cdot Char^{daf}_s \cdot 1000,\qquad
CO^{Nm3}_s = Y^{CO}_s \cdot Char^{daf}_s \cdot 1000
$$
$$
H2^{kmol}_s = \frac{H2^{Nm3}_s}{\nu},\qquad
CO^{kmol}_s = \frac{CO^{Nm3}_s}{\nu}
$$

---

### (4) FT conversion, no WGS, H2 make-up
$$
0 \le CO^{use}_s \le CO^{kmol}_s
$$
$$
H2^{make,kmol}_s \ge R^{FT} \cdot CO^{use}_s - H2^{kmol}_s,\qquad H2^{make,kmol}_s\ge 0
$$
$$
H2^{make}_s = \frac{2}{1000} H2^{make,kmol}_s
$$
$$
L^{liq}_s = \eta^{liq}_s \cdot CO^{use}_s
$$

---

### (5) Product split and jet capacity
Base products:
$$
Jet^{base}_s=f^{jet}_s L^{liq}_s,\;
Dies^{base}_s=f^{dies}_s L^{liq}_s,\;
Nap^{base}_s=f^{nap}_s L^{liq}_s,\;
Wax_s=f^{wax}_s L^{liq}_s
$$

Wax upgrading (Big-M logic):
$$
0\le Wax^{up}_s \le Wax_s,\qquad Wax^{up}_s \le M^{wax}u
$$
$$
Wax^{sell}_s = Wax_s - Wax^{up}_s
$$

Upgraded products:
$$
Jet_s = Jet^{base}_s + \gamma^{jet}_s Wax^{up}_s
$$
$$
Dies_s = Dies^{base}_s + \gamma^{dies}_s Wax^{up}_s,\qquad
Nap_s = Nap^{base}_s + \gamma^{nap}_s Wax^{up}_s
$$

Jet capacity:
$$
Jet_s \le Cap^{jet}\qquad \forall s
$$

---

### (6) Oxygen requirement and electricity peak sizing
Oxygen:
$$
O2_s = \rho_{O2} B^{dry}_s
$$

Electricity balance (implementation uses linear intensities):
- grinding/feedpress based on thermal input from $LHV_{dry}$
- ASU: $k^{asu}_s \cdot O2_s$
- FT: $k^{ft} \cdot L^{liq}_s$
- upgrading: $k^{up} \cdot Wax^{up}_s$

Peak demand must cover all scenarios:
$$
E_s \le H \cdot P^{peak}\qquad \forall s
$$

---

## Economics (per scenario)

### Biomass cost (scenario-independent)
$$
C^{bio}=\sum_r (c^{bale}+c^{tr}d_r)\,B_r
$$

### Electricity: two-part tariff
$$
C^{elec}_s = p^{e}E_s + 12\,p^{d}P^{peak}
$$

### Other OPEX and revenue
$$
C^{H2}_s = p^{H2}\,H2^{make}_s,\qquad
C^{heat}_s=p^{h}Q^{dry}_s,\qquad
C^{FOM}=C^{FOM}\,y
$$
$$
Rev_s = P^{jet}Jet_s + P^{dies}Dies_s + P^{nap}Nap_s + P^{wax}Wax^{sell}_s
$$
$$
Profit_s = Rev_s - \left(C^{bio}+C^{elec}_s+C^{H2}_s+C^{heat}_s+C^{FOM}\right)
$$

Expected annual profit:
$$
\mathbb{E}[Profit]=\sum_s \pi_s Profit_s
$$

---

## NPV and objective

Gate-horizon NPV:
$$
NPV_{8} = -Capex0 + PVF(8)\cdot \mathbb{E}[Profit]
$$

Lifetime NPV:
$$
NPV_{20} = -Capex0 + PVF(20)\cdot \mathbb{E}[Profit]
$$

Financing gate (policy run):
$$
NPV_{8} \ge 0
$$

Objective:
$$
\max \; NPV_{20}
$$

---

## Policy vs counterfactual reporting (important)

- **Policy solve** applies the gate $NPV_8\ge 0$ → if not financeable, optimal decision is typically **$y=0$** (no build, no production).
- For analysis/plots, the code can also run **counterfactual** solves (force build and/or relax gate) to report:
  - best operational plan even if not financeable,
  - both wax-upgrading options for comparison.

This is intentional so you can see “how much it would produce / lose money” even when the investment decision is “do not build”.
