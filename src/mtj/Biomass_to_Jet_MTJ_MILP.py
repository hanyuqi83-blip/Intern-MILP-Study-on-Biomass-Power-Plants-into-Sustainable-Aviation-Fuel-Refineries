from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import pulp  # type: ignore
    HAS_PULP = True
except Exception:
    HAS_PULP = False


# -------------------------
# 工具函数
# -------------------------

def resolve_relative_to_script(p: Path) -> Path:
    """把相对路径解析为相对于当前脚本所在目录的绝对路径。"""
    if p.is_absolute():
        return p
    return (Path(__file__).resolve().parent / p).resolve()


def load_json(path: Path) -> Dict[str, Any]:
    """读取 JSON 配置。"""
    return json.loads(path.read_text(encoding="utf-8"))


def pv_factor(i: float, delta: float, years: int) -> float:
    """计算带退化的年金现值系数。"""
    q = (1.0 - delta) / (1.0 + i)
    if abs(1.0 - q) < 1e-12:
        return years / (1.0 + i)
    return (1.0 / (1.0 + i)) * (1.0 - (q ** years)) / (1.0 - q)


def ring_supply_scaled(data: Dict[str, Any], bio_total_override: Optional[float]) -> List[float]:
    """按圈层最大供给得到供应上限；如指定 bio_total_override 则按比例缩放。"""
    rings = data["biomass_rings"]
    ring_supplies = [float(r["max_supply_tpy"]) for r in rings]
    if bio_total_override is None:
        return ring_supplies
    s0 = sum(ring_supplies) if sum(ring_supplies) > 0 else 1.0
    scale = float(bio_total_override) / s0
    return [x * scale for x in ring_supplies]


def _get_o2_ratio(proc: Dict[str, Any]) -> float:
    """从 process 中取气化耗氧比，兼容 gasification / entrained_flow 两种命名。"""
    if "gasification" in proc and "O2_to_dry_biomass_mass_ratio" in proc["gasification"]:
        return float(proc["gasification"]["O2_to_dry_biomass_mass_ratio"])
    if "entrained_flow" in proc and "O2_to_dry_biomass_mass_ratio" in proc["entrained_flow"]:
        return float(proc["entrained_flow"]["O2_to_dry_biomass_mass_ratio"])
    raise ValueError("Missing O2_to_dry_biomass_mass_ratio in process.gasification or process.entrained_flow")


def validate_data(data: Dict[str, Any]) -> None:
    """检查关键字段是否齐全，避免模型在缺字段时静默退化。"""
    need = ["process", "feedstock", "economic", "capex", "biomass_rings", "scenarios"]
    for k in need:
        if k not in data:
            raise ValueError(f"Missing top-level key: {k}")

    proc = data["process"]
    econ = data["economic"]
    capex = data["capex"]
    feed = data["feedstock"]

    # 干燥参数
    xin = float(proc["moisture_in_wb"])
    xout = float(proc["moisture_out_wb"])
    if not (0.0 <= xin < 1.0) or not (0.0 <= xout < 1.0):
        raise ValueError("moisture_in_wb/moisture_out_wb must be in [0,1)")
    if xout > xin:
        raise ValueError("moisture_out_wb > moisture_in_wb: drying increases moisture?")

    # 气化耗氧比
    _ = _get_o2_ratio(proc)

    # 合成气体积-物质的量换算常数
    if "constants" not in proc or "Nm3_per_kmol" not in proc["constants"]:
        raise ValueError("process.constants.Nm3_per_kmol not found")

    # 灰分用于 char_daf 计算
    if "ultimate_dry_wt_percent" not in feed or "Ash" not in feed["ultimate_dry_wt_percent"]:
        raise ValueError("feedstock.ultimate_dry_wt_percent.Ash not found")

    # 产能选项
    caps = list(map(float, capex["capacity_options_meoh_tpy"]))
    capex_vals = list(map(float, capex["capex0_RMB_for_option"]))
    if len(caps) != len(capex_vals):
        raise ValueError("capex option length mismatch")

    # 场景概率
    psum = sum(float(sc.get("probability", 0.0)) for sc in data["scenarios"])
    if abs(psum - 1.0) > 1e-6:
        raise ValueError(f"Scenario probabilities must sum to 1.0, got {psum}")

    # 汇率
    if float(econ["FX_RMB_per_USD"]) <= 0:
        raise ValueError("FX_RMB_per_USD must be > 0")

    # 合成气产率字段
    for sc in data["scenarios"]:
        for k in ("char_yield_dry_fraction", "H2_yield_Nm3_per_kg_char_daf", "CO_yield_Nm3_per_kg_char_daf"):
            if k not in sc:
                raise ValueError(f"Scenario '{sc.get('name','?')}' missing key: {k}")

        if "meoh_CO_to_meoh_eff" in sc:
            eta = float(sc["meoh_CO_to_meoh_eff"])
            if not (0.0 <= eta <= 1.0):
                raise ValueError(f"Scenario '{sc.get('name','?')}' meoh_CO_to_meoh_eff must be in [0,1]")


# -------------------------
# 输出结构
# -------------------------

@dataclass
class RowOut:
    scenario: str
    probability: float

    build: int
    cap_meoh_tpy: float
    mtj_build: int

    biomass_in_ar_tpy: float
    dry_solids_tpy: float
    o2_required_tpy: float

    syngas_CO_Nm3_y: float
    syngas_H2_Nm3_y: float
    CO_used_kmol_y: float

    methanol_tpy: float
    methanol_to_mtj_tpy: float
    methanol_sold_tpy: float

    jet_tpy: float
    gasoline_tpy: float
    naphtha_tpy: float

    electricity_kWh_y: float
    peak_kW: float

    h2_makeup_tpy: float
    h2_mtj_tpy: float
    h2_total_tpy: float

    revenue_RMB_y: float
    opex_RMB_y: float
    annual_profit_RMB_y: float


def write_csv(rows: List[RowOut], out_path: Path) -> None:
    """把逐情景结果写入 CSV。"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()) if rows else ["empty"])
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))


# -------------------------
# MILP
# -------------------------

def solve_milp_mtj(
    data: Dict[str, Any],
    bio_total_override: Optional[float] = None,
    force_build: bool = False,
    relax_gate: bool = False,
    fix_mtj: Optional[int] = None,
    min_load_frac: Optional[float] = None,
    risk_alpha: Optional[float] = None,
    risk_lambda: float = 0.0,
) -> Dict[str, Any]:
    if not HAS_PULP:
        raise RuntimeError("PuLP not installed. Please install pulp (and CBC).")

    validate_data(data)

    proc = data["process"]
    feed = data["feedstock"]
    econ = data["economic"]
    capex = data["capex"]

    # 产能选项仅保留正值；y=0 时由约束保证不选任何产能
    caps_all = list(map(float, capex["capacity_options_meoh_tpy"]))
    capex_vals_all = list(map(float, capex["capex0_RMB_for_option"]))
    caps: List[float] = []
    capex_vals: List[float] = []
    for c, v in zip(caps_all, capex_vals_all):
        if c > 0:
            caps.append(c)
            capex_vals.append(v)
    if not caps:
        raise ValueError("Need at least one positive capacity option in capacity_options_meoh_tpy")

    cap_min = min(caps)
    cap_max = max(caps)

    rings = data["biomass_rings"]
    ring_supplies = ring_supply_scaled(data, bio_total_override)
    Bmax = sum(ring_supplies)

    # 现值系数
    i = float(econ["discount_rate"])
    delta = float(econ["degradation_rate"])
    gate_years = int(econ.get("gate_years", 8))
    life_years = int(econ.get("lifetime_years", 20))
    PVF_gate = pv_factor(i, delta, gate_years)
    PVF_life = pv_factor(i, delta, life_years)

    prob = pulp.LpProblem("Biomass_to_MeOH_to_Jet_MTJ_syngas", pulp.LpMaximize)

    # 一阶段：建厂与离散产能
    y = pulp.LpVariable("build", lowBound=0, upBound=1, cat=pulp.LpBinary)
    if force_build:
        prob += y == 1, "force_build"

    x = [pulp.LpVariable(f"cap_choice_{int(c)}", lowBound=0, upBound=1, cat=pulp.LpBinary) for c in caps]
    prob += pulp.lpSum(x) == y, "choose_one_capacity_if_built"

    cap_meoh = pulp.LpVariable("cap_meoh_tpy", lowBound=0)
    prob += cap_meoh == pulp.lpSum(c * xj for c, xj in zip(caps, x)), "cap_meoh_def"
    prob += cap_meoh >= cap_min * y, "cap_meoh_min_if_built"
    prob += cap_meoh <= cap_max * y, "cap_meoh_max_if_built"

    # 一阶段：MTJ 装置开关（只有建甲醇厂才允许建 MTJ）
    u = pulp.LpVariable("mtj_build", lowBound=0, upBound=1, cat=pulp.LpBinary)
    prob += u <= y, "mtj_only_if_built"
    if fix_mtj is not None:
        prob += u == int(fix_mtj), "fix_mtj"

    # 一阶段：按圈层采购生物质
    B_r = [pulp.LpVariable(f"bio_ring_{k}_tpy", lowBound=0) for k in range(len(rings))]
    for k, sup in enumerate(ring_supplies):
        prob += B_r[k] <= sup * y, f"ring_supply_{k}"

    B = pulp.LpVariable("biomass_in_ar_tpy", lowBound=0)
    prob += B == pulp.lpSum(B_r), "biomass_total"

    # 一阶段：避免建厂但不进料（可选）
    if min_load_frac is not None:
        frac = float(min_load_frac)
        if not (0.0 <= frac <= 1.0):
            raise ValueError("--min-load-frac must be in [0,1]")
        prob += B >= frac * Bmax * y, "min_load_biomass"

    # 一阶段：需量电峰值按最坏情景定容
    P_peak = pulp.LpVariable("P_peak_kW", lowBound=0)

    # 一阶段：CAPEX（产能选项 + MTJ 增量）
    capex0 = pulp.LpVariable("capex0_RMB", lowBound=0)
    base_capex0 = pulp.lpSum(v * xj for v, xj in zip(capex_vals, x))
    mtj_capex0 = float(capex.get("mtj_capex0_RMB", 0.0)) * u
    prob += capex0 == base_capex0 + mtj_capex0, "capex0_def"

    # 常数：干燥、灰分、耗氧、合成气换算
    x_in = float(proc["moisture_in_wb"])
    x_out = float(proc["moisture_out_wb"])
    ash_pct = float(feed["ultimate_dry_wt_percent"]["Ash"]) / 100.0
    o2_ratio = _get_o2_ratio(proc)
    Nm3_per_kmol = float(proc["constants"]["Nm3_per_kmol"])

    # 常数：甲醇合成化学计量与分子量
    R_meoh = 2.0            # CO + 2H2 -> CH3OH
    M_meoh_kg_per_kmol = 32.0

    # 常数：经济参数
    bailing = float(econ["biomass_base_cost_RMB_per_t"])
    transport = float(econ["biomass_transport_RMB_per_t_km"])
    p_kWh = float(econ["electricity_energy_RMB_per_kWh"])
    p_dem = float(econ["electricity_demand_RMB_per_kW_month"])
    H_op = float(econ["electricity_operating_hours_per_year"])
    p_heat = float(econ.get("heat_RMB_per_MJ", 0.0))
    p_H2 = float(econ["H2_RMB_per_t"])

    load_factor = float(econ.get("peak_load_factor", 1.0))
    load_factor = max(min(load_factor, 1.0), 1e-3)

    fixed_om = float(econ.get("fixed_OM_RMB_per_y", 0.0)) * y
    mtj_fixed_om = float(econ.get("mtj_fixed_OM_RMB_per_y", 0.0)) * u

    fx = float(econ["FX_RMB_per_USD"])
    P_meoh = float(econ["product_prices"]["methanol_USD_per_t"]) * fx
    P_jet = float(econ["product_prices"]["jet_USD_per_t"]) * fx
    P_gas = float(econ["product_prices"].get("gasoline_USD_per_t", 0.0)) * fx
    P_nap = float(econ["product_prices"].get("naphtha_USD_per_t", 0.0)) * fx

    # 常数：生物质采购成本按圈层距离累加
    C_bio = 0
    for Br, ring in zip(B_r, rings):
        d = float(ring["distance_km"])
        C_bio += (bailing + transport * d) * Br

    # 常数：把 MTJ 入口甲醇流量 big-M 收紧到产能上限
    M_meoh_to_mtj = cap_max

    scenario_profit: Dict[str, Tuple[float, Any]] = {}
    scenario_rows: Dict[str, Dict[str, Any]] = {}

    # 二阶段：逐情景物料衡算、能耗、经济
    for sc in data["scenarios"]:
        name = str(sc["name"])
        ps = float(sc["probability"])

        # 干燥相关变量
        dry_solids = pulp.LpVariable(f"dry_solids_tpy__{name}", lowBound=0)
        after_dry = pulp.LpVariable(f"after_dry_ar_tpy__{name}", lowBound=0)
        water_removed = pulp.LpVariable(f"water_removed_tpy__{name}", lowBound=0)

        # 耗氧
        o2_req = pulp.LpVariable(f"o2_required_tpy__{name}", lowBound=0)

        # char 与合成气
        char_dry = pulp.LpVariable(f"char_dry_tpy__{name}", lowBound=0)
        char_daf = pulp.LpVariable(f"char_daf_tpy__{name}", lowBound=0)

        syngas_H2_Nm3 = pulp.LpVariable(f"syngas_H2_Nm3_y__{name}", lowBound=0)
        syngas_CO_Nm3 = pulp.LpVariable(f"syngas_CO_Nm3_y__{name}", lowBound=0)

        H2_kmol = pulp.LpVariable(f"H2_kmol_y__{name}", lowBound=0)
        CO_kmol = pulp.LpVariable(f"CO_kmol_y__{name}", lowBound=0)

        # 甲醇合成使用的 CO 与补氢
        CO_used = pulp.LpVariable(f"CO_used_kmol_y__{name}", lowBound=0)
        H2_makeup_kmol = pulp.LpVariable(f"H2_makeup_kmol_y__{name}", lowBound=0)
        H2_makeup_tpy = pulp.LpVariable(f"H2_makeup_tpy__{name}", lowBound=0)

        # 甲醇产量
        meoh_raw = pulp.LpVariable(f"methanol_raw_tpy__{name}", lowBound=0)
        meoh = pulp.LpVariable(f"methanol_tpy__{name}", lowBound=0)

        # 甲醇分配：卖甲醇 vs 进 MTJ
        meoh_to_mtj = pulp.LpVariable(f"methanol_to_mtj_tpy__{name}", lowBound=0)
        meoh_sold = pulp.LpVariable(f"methanol_sold_tpy__{name}", lowBound=0)

        # MTJ 产品
        jet = pulp.LpVariable(f"jet_tpy__{name}", lowBound=0)
        gasoline = pulp.LpVariable(f"gasoline_tpy__{name}", lowBound=0)
        naphtha = pulp.LpVariable(f"naphtha_tpy__{name}", lowBound=0)

        # 能耗与氢耗
        E_kWh = pulp.LpVariable(f"electricity_kWh_y__{name}", lowBound=0)
        H2_mtj_tpy = pulp.LpVariable(f"H2_mtj_tpy__{name}", lowBound=0)
        H2_total_tpy = pulp.LpVariable(f"H2_total_tpy__{name}", lowBound=0)

        # 干燥物料衡算
        prob += dry_solids == B * (1.0 - x_in), f"dry_solids_def__{name}"
        prob += after_dry == dry_solids / (1.0 - x_out), f"after_dry_def__{name}"
        prob += water_removed == B - after_dry, f"water_removed_def__{name}"
        prob += water_removed >= 0, f"water_removed_nonneg__{name}"

        # 氧耗按干基比例估算
        prob += o2_req == o2_ratio * dry_solids, f"o2_req_def__{name}"

        # char 生成与 daf 修正
        y_char = float(sc["char_yield_dry_fraction"])
        prob += char_dry == y_char * dry_solids, f"char_dry_def__{name}"
        prob += char_daf == char_dry - ash_pct * dry_solids, f"char_daf_def__{name}"
        prob += char_daf >= 0, f"char_daf_nonneg__{name}"

        # 合成气产量（由 daf char 的体积产率得到）
        Y_H2 = float(sc["H2_yield_Nm3_per_kg_char_daf"])
        Y_CO = float(sc["CO_yield_Nm3_per_kg_char_daf"])
        prob += syngas_H2_Nm3 == Y_H2 * char_daf * 1000.0, f"syngas_H2_def__{name}"
        prob += syngas_CO_Nm3 == Y_CO * char_daf * 1000.0, f"syngas_CO_def__{name}"

        # 合成气体积转为 kmol
        prob += H2_kmol == syngas_H2_Nm3 / Nm3_per_kmol, f"H2_kmol_def__{name}"
        prob += CO_kmol == syngas_CO_Nm3 / Nm3_per_kmol, f"CO_kmol_def__{name}"

        # 使用的 CO 不能超过可用 CO
        prob += CO_used <= CO_kmol, f"CO_used_le_CO__{name}"

        # 补氢满足甲醇合成化学计量
        prob += H2_makeup_kmol >= R_meoh * CO_used - H2_kmol, f"H2_makeup_stoich__{name}"
        prob += H2_makeup_kmol >= 0, f"H2_makeup_nonneg__{name}"
        prob += H2_makeup_tpy == H2_makeup_kmol * 2.0 / 1000.0, f"H2_makeup_mass__{name}"

        # 由 CO_used 生成甲醇（效率用于反映净化/合成段的 CO 利用率）
        eta_CO = float(sc.get("meoh_CO_to_meoh_eff", 1.0))
        prob += meoh_raw == eta_CO * CO_used * (M_meoh_kg_per_kmol / 1000.0), f"meoh_raw_from_CO__{name}"

        # 甲醇产量受“合成气上限”和“装置产能上限”共同约束
        prob += meoh <= meoh_raw, f"meoh_le_raw__{name}"
        prob += meoh <= cap_meoh, f"meoh_cap__{name}"

        # 甲醇分配：只有建 MTJ 才允许送入 MTJ
        prob += meoh_to_mtj <= meoh, f"meoh_to_mtj_le_meoh__{name}"
        prob += meoh_to_mtj <= M_meoh_to_mtj * u, f"meoh_to_mtj_only_if_mtj__{name}"
        prob += meoh_sold == meoh - meoh_to_mtj, f"meoh_sold_def__{name}"

        # MTJ 产品产率
        y_jet = float(proc["mtj"]["yields_t_per_t_meoh"]["jet"])
        y_gas = float(proc["mtj"]["yields_t_per_t_meoh"].get("gasoline", 0.0))
        y_nap = float(proc["mtj"]["yields_t_per_t_meoh"].get("naphtha", 0.0))
        prob += jet == y_jet * meoh_to_mtj, f"jet_def__{name}"
        prob += gasoline == y_gas * meoh_to_mtj, f"gasoline_def__{name}"
        prob += naphtha == y_nap * meoh_to_mtj, f"naphtha_def__{name}"

        # 年电量：甲醇段电耗 + ASU 电耗 + MTJ 电耗
        e_meoh = float(sc["electricity_kWh_per_t_meoh"])
        e_asu = float(sc.get("ASU_kWh_per_t_O2", 0.0))
        e_mtj = float(proc["mtj"]["electricity_kWh_per_t_meoh_to_mtj"])
        prob += E_kWh == (e_meoh * meoh) + (e_asu * o2_req) + (e_mtj * meoh_to_mtj), f"elec_def__{name}"

        # 峰值需量覆盖：用 load_factor 把年电量映射到峰值合同容量
        prob += E_kWh <= (H_op * load_factor) * P_peak, f"peak_cover__{name}"

        # 干燥热耗成本
        q_dry = float(sc["drying_heat_MJ_per_kg_water_removed"])
        drying_heat_MJ = q_dry * water_removed * 1000.0
        C_heat = p_heat * drying_heat_MJ

        # MTJ 加氢消耗
        h2_kg_per_t_meoh_to_mtj = float(proc["mtj"]["H2_kg_per_t_meoh_to_mtj"])
        prob += H2_mtj_tpy == (h2_kg_per_t_meoh_to_mtj * meoh_to_mtj) / 1000.0, f"H2_mtj_def__{name}"

        # 总外购氢 = 甲醇补氢 + MTJ 加氢
        prob += H2_total_tpy == H2_makeup_tpy + H2_mtj_tpy, f"H2_total_def__{name}"

        # 成本与收入
        C_elec = p_kWh * E_kWh + 12.0 * p_dem * P_peak
        C_H2 = p_H2 * H2_total_tpy

        OPEX = C_bio + C_elec + C_heat + C_H2 + fixed_om + mtj_fixed_om
        REV = P_meoh * meoh_sold + P_jet * jet + P_gas * gasoline + P_nap * naphtha
        PROFIT = REV - OPEX

        scenario_profit[name] = (ps, PROFIT)
        scenario_rows[name] = dict(
            dry_solids=dry_solids,
            o2_req=o2_req,
            syngas_CO_Nm3=syngas_CO_Nm3,
            syngas_H2_Nm3=syngas_H2_Nm3,
            CO_used=CO_used,
            meoh=meoh,
            meoh_to_mtj=meoh_to_mtj,
            meoh_sold=meoh_sold,
            jet=jet,
            gasoline=gasoline,
            naphtha=naphtha,
            E_kWh=E_kWh,
            H2_makeup_tpy=H2_makeup_tpy,
            H2_mtj_tpy=H2_mtj_tpy,
            H2_total_tpy=H2_total_tpy,
            REV=REV,
            OPEX=OPEX,
            PROFIT=PROFIT,
        )

    # 期望利润
    exp_profit = pulp.lpSum(ps * prof for ps, prof in scenario_profit.values())

    # 可选：对年度利润做 CVaR 下行尾部度量
    if risk_alpha is not None and risk_lambda > 0:
        a = float(risk_alpha)
        if not (0.0 < a < 1.0):
            raise ValueError("--risk-alpha must be in (0,1)")
        if risk_lambda < 0:
            raise ValueError("--risk-lambda must be >= 0")

        z = pulp.LpVariable("VaR_profit_z", lowBound=None)
        eta = {sn: pulp.LpVariable(f"shortfall_eta__{sn}", lowBound=0) for sn in scenario_profit.keys()}
        for sn, (pi, prof) in scenario_profit.items():
            prob += eta[sn] >= z - prof, f"cvar_shortfall__{sn}"

        CVaR_profit = pulp.LpVariable("CVaR_profit", lowBound=None)
        prob += CVaR_profit == z - (1.0 / (1.0 - a)) * pulp.lpSum(
            float(scenario_profit[sn][0]) * eta[sn] for sn in scenario_profit.keys()
        ), "CVaR_def"

        risk_adj_profit = exp_profit - float(risk_lambda) * (exp_profit - CVaR_profit)
    else:
        CVaR_profit = None
        risk_adj_profit = exp_profit

    # NPV 定义
    NPV_gate = pulp.LpVariable("NPV_gate_RMB", lowBound=None)
    NPV_life = pulp.LpVariable("NPV_life_RMB", lowBound=None)
    prob += NPV_gate == -capex0 + PVF_gate * exp_profit, "NPV_gate_def"
    prob += NPV_life == -capex0 + PVF_life * risk_adj_profit, "NPV_life_def"

    # 融资门槛
    if not relax_gate:
        prob += NPV_gate >= 0, "financing_gate"

    # 目标函数
    prob += NPV_life, "maximize_NPV_life"

    # 求解
    solver = pulp.PULP_CBC_CMD(msg=False)
    prob.solve(solver)
    status = pulp.LpStatus.get(prob.status, str(prob.status))

    def val(v) -> float:
        try:
            return float(pulp.value(v))
        except Exception:
            return float("nan")

    # 汇总一阶段解
    build = int(round(val(y)))
    cap_meoh_sol = val(cap_meoh)
    mtj_sol = int(round(val(u)))
    B_sol = val(B)
    peak_sol = val(P_peak)
    capex0_sol = val(capex0)
    exp_profit_sol = val(exp_profit)
    npv_gate_sol = val(NPV_gate)
    npv_life_sol = val(NPV_life)

    # 汇总逐情景解
    rows_out: List[RowOut] = []
    for sc in data["scenarios"]:
        name = str(sc["name"])
        vs = scenario_rows[name]
        rows_out.append(RowOut(
            scenario=name,
            probability=float(sc["probability"]),
            build=build,
            cap_meoh_tpy=cap_meoh_sol,
            mtj_build=mtj_sol,
            biomass_in_ar_tpy=B_sol,
            dry_solids_tpy=val(vs["dry_solids"]),
            o2_required_tpy=val(vs["o2_req"]),
            syngas_CO_Nm3_y=val(vs["syngas_CO_Nm3"]),
            syngas_H2_Nm3_y=val(vs["syngas_H2_Nm3"]),
            CO_used_kmol_y=val(vs["CO_used"]),
            methanol_tpy=val(vs["meoh"]),
            methanol_to_mtj_tpy=val(vs["meoh_to_mtj"]),
            methanol_sold_tpy=val(vs["meoh_sold"]),
            jet_tpy=val(vs["jet"]),
            gasoline_tpy=val(vs["gasoline"]),
            naphtha_tpy=val(vs["naphtha"]),
            electricity_kWh_y=val(vs["E_kWh"]),
            peak_kW=peak_sol,
            h2_makeup_tpy=val(vs["H2_makeup_tpy"]),
            h2_mtj_tpy=val(vs["H2_mtj_tpy"]),
            h2_total_tpy=val(vs["H2_total_tpy"]),
            revenue_RMB_y=val(vs["REV"]),
            opex_RMB_y=val(vs["OPEX"]),
            annual_profit_RMB_y=val(vs["PROFIT"]),
        ))

    res: Dict[str, Any] = dict(
        status=status,
        build=build,
        cap_meoh_tpy=cap_meoh_sol,
        mtj_build=mtj_sol,
        biomass_in_ar_tpy=B_sol,
        P_peak_kW=peak_sol,
        capex0_RMB=capex0_sol,
        expected_annual_profit_RMB_y=exp_profit_sol,
        PVF_gate=PVF_gate,
        PVF_life=PVF_life,
        NPV_gate_RMB=npv_gate_sol,
        NPV_life_RMB=npv_life_sol,
        rows=[asdict(r) for r in rows_out],
    )

    if CVaR_profit is not None:
        res["CVaR_profit_RMB_y"] = val(CVaR_profit)
        res["risk_alpha"] = float(risk_alpha)
        res["risk_lambda"] = float(risk_lambda)

    return res


# -------------------------
# CLI
# -------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=str, default="Biomass_to_Jet_MTJ_Data_syngas.json")
    ap.add_argument("--out", type=str, default="meoh_jet_results.csv")
    ap.add_argument("--bio-total", type=float, default=None)
    ap.add_argument("--min-load-frac", type=float, default=None)
    ap.add_argument("--report-counterfactual", action="store_true")
    ap.add_argument("--report-both-mtj", action="store_true")
    ap.add_argument("--risk-alpha", type=float, default=None)
    ap.add_argument("--risk-lambda", type=float, default=0.0)

    args = ap.parse_args(argv)

    data_path = resolve_relative_to_script(Path(args.data))
    if not data_path.exists():
        raise SystemExit(f"[ERROR] Data file not found: {data_path}")

    data = load_json(data_path)

    out_base = resolve_relative_to_script(Path(args.out))

    def _suffix_path(sfx: str) -> Path:
        return out_base.with_name(out_base.stem + sfx + out_base.suffix)

    def _write_one(res_obj: Dict[str, Any], out_path: Path) -> None:
        rows = [RowOut(**r) for r in res_obj["rows"]]  # type: ignore
        if rows:
            write_csv(rows, out_path)
            print(f"\nWrote results CSV to: {out_path}")

    bundle: Dict[str, Any] = {}

    # policy：融资门槛开启
    res_policy = solve_milp_mtj(
        data,
        bio_total_override=args.bio_total,
        force_build=False,
        relax_gate=False,
        fix_mtj=None,
        min_load_frac=args.min_load_frac,
        risk_alpha=args.risk_alpha,
        risk_lambda=args.risk_lambda,
    )
    bundle["policy"] = res_policy
    _write_one(res_policy, _suffix_path("_policy"))

    # counterfactual：强制建厂 + 关闭融资门槛
    if args.report_counterfactual:
        res_cf_best = solve_milp_mtj(
            data,
            bio_total_override=args.bio_total,
            force_build=True,
            relax_gate=True,
            fix_mtj=None,
            min_load_frac=args.min_load_frac,
            risk_alpha=args.risk_alpha,
            risk_lambda=args.risk_lambda,
        )
        bundle["counterfactual_best"] = res_cf_best
        _write_one(res_cf_best, _suffix_path("_cf_best"))

        if args.report_both_mtj:
            res_cf_mtj0 = solve_milp_mtj(
                data,
                bio_total_override=args.bio_total,
                force_build=True,
                relax_gate=True,
                fix_mtj=0,
                min_load_frac=args.min_load_frac,
                risk_alpha=args.risk_alpha,
                risk_lambda=args.risk_lambda,
            )
            res_cf_mtj1 = solve_milp_mtj(
                data,
                bio_total_override=args.bio_total,
                force_build=True,
                relax_gate=True,
                fix_mtj=1,
                min_load_frac=args.min_load_frac,
                risk_alpha=args.risk_alpha,
                risk_lambda=args.risk_lambda,
            )
            bundle["counterfactual_mtj_0"] = res_cf_mtj0
            bundle["counterfactual_mtj_1"] = res_cf_mtj1
            _write_one(res_cf_mtj0, _suffix_path("_cf_mtj0"))
            _write_one(res_cf_mtj1, _suffix_path("_cf_mtj1"))

    print(json.dumps(bundle, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
