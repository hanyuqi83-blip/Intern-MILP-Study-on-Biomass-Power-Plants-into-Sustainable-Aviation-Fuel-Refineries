"""
本代码实现一个“系统级过程—经济”一体化的 MILP/随机规划模型，用于评估并优化
生物质→合成气→费托（FT）路线生产可持续航空燃料（SAF）的投资与运营决策。

────────────────────────────────────────────────────────
一、建模方法（How）
────────────────────────────────────────────────────────
1) 两阶段随机规划（Two-stage stochastic programming）
   - 第一阶段（First-stage）：做跨情景共享的投资/合同/供应链决策（是否建厂、离散产能选项、
     是否配置蜡升级、各距离圈层生物质采购量、需量电合同峰值 P_peak 等）。
   - 第二阶段（Second-stage）：在每个不确定情景 s 下做运营决策（物料衡算、合成气产率、
     FT 转化、产品分配、用电/用热、公用工程消耗、收入与 OPEX），得到情景利润 Profit_s。
   - 通过情景概率 π_s 将各情景利润汇总为期望利润 E[Profit]，并据此计算 NPV。

2) 混合整数线性规划（MILP）
   - “整数/离散”来自：建厂与否、离散产能选项选择、是否建设蜡升级（0/1）。
   - “线性”来自：以线性收率/线性强度描述干燥、气化、ASU、FT 与升级过程，以及线性化的成本与收益。
   - 求解由 PuLP 建模并调用 MILP 求解器（默认 CBC）完成。

3) 金融评价（NPV + 融资门槛）
   - 年度期望利润 → 通过折现因子 PVF(N) 计算 N 年期 NPV：
       NPV_N = -Capex0 + PVF(N) * E[Profit]
   - 目标：最大化 NPV_20（长期价值最大化）。
   - 约束：NPV_8 ≥ 0（融资可行性门槛），若无法满足，模型可选择 build=0（不建设）作为最优。

4) 两部制电价（Two-part tariff）
   - 电费由“电量电费 + 最大需量电费”组成：
       C_elec_s = p_energy * E_s + 12 * p_demand * P_peak
   - P_peak 为第一阶段决策（合同需量/接入容量），并通过约束保证其对所有情景的用电需求足够。

────────────────────────────────────────────────────────
二、实现的对象（What）
────────────────────────────────────────────────────────
1) 工艺链路（Process）
   生物质（湿基到厂态）→ 干燥/粉碎 → ASU 供氧 → Entrained-flow 气化 → 合成气（H2/CO）
   → H2 补充（默认不启用 WGS/RWGS）→ FT 合成 → 产品分配（Jet/Diesel/Naphtha/Wax）
   →（可选）蜡升级：将 Wax 转化为更多 Jet/Diesel/Naphtha。

2) 决策问题（Decisions）
   - 是否建设装置（build or not）
   - 选择离散航煤产能规模（capacity option）
   - 是否配置蜡升级单元（wax upgrading）
   - 生物质供应链：从不同距离圈层采购多少生物质（ring-based procurement）
   - 运营层：每个情景下合成气进入 FT 的规模、产品产量、用电/用热与成本收益等。

────────────────────────────────────────────────────────
三、输出结果（Outputs）
────────────────────────────────────────────────────────
- 最优/固定方案下：各情景的产量（Jet/Diesel/Naphtha/Wax）、关键物耗（O2、H2、用电等）、
  年度收入、OPEX、利润 Profit_s；
- 期望年度利润 E[Profit]；
- 8 年期与 20 年期 NPV（NPV_8、NPV_20）；
- 若启用求解器：返回最优的投资与运营方案；若不启用求解器：返回给定固定方案的可复现计算结果。

────────────────────────────────────────────────────────
四、推荐运行方式（Run）
────────────────────────────────────────────────────────
- MILP 优化求解（需要安装 PuLP 与 CBC 或其他 MILP 求解器）：
  python Biomass_to_Jet_FT_MILP.py --data Biomass_to_Jet_FT_Data.json --out ft_results.csv

- 固定方案评估（不求解，仅用于复现/校验；使用 model.fixed 给定的决策）：
  python Biomass_to_Jet_FT_MILP.py --data Biomass_to_Jet_FT_Data.json --out ft_results.csv --no-solver
"""

from __future__ import annotations
import argparse
import csv
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import pulp  # type: ignore
    HAS_PULP = True
except Exception:
    HAS_PULP = False


# -----------------------------
# Helpers / Validation
# -----------------------------

def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))

def safe_div(a: float, b: float) -> float:
    return a / b if abs(b) > 1e-12 else 0.0

def pv_factor(i: float, delta: float, years: int) -> float:
    """
    PVF_N = sum_{t=1..N} (1/(1+i)^t) * (1-delta)^(t-1)
    Closed form with q=(1-delta)/(1+i).
    """
    q = (1.0 - delta) / (1.0 + i)
    if abs(1.0 - q) < 1e-12:
        # 极端情况：q≈1，退化为 years/(1+i)
        return years / (1.0 + i)
    return (1.0 / (1.0 + i)) * (1.0 - (q ** years)) / (1.0 - q)

def resolve_relative_to_script(p: Path) -> Path:
    if p.is_absolute():
        return p
    return (Path(__file__).resolve().parent / p).resolve()

def normalize_split(fracs: Dict[str, float], keys: Tuple[str, ...]) -> Dict[str, float]:
    """把 product split 归一化到和=1，防止数据小误差导致模型不稳。"""
    vals = [float(fracs.get(k, 0.0)) for k in keys]
    s = sum(vals)
    if s <= 0:
        # 极端：全 0，给一个默认（不建议，但保证代码不炸）
        n = len(keys)
        return {k: 1.0 / n for k in keys}
    return {k: v / s for k, v in zip(keys, vals)}

def validate_data(data: Dict[str, Any]) -> None:
    """
    对数据做基本一致性与物理可行性检查。
    这一步能显著提升“可复现性 + 可解释性”，也更像研究生作品。
    """
    required_top = ["process", "feedstock", "economic", "capex", "biomass_rings", "scenarios", "model"]
    for k in required_top:
        if k not in data:
            raise ValueError(f"Missing top-level key: {k}")

    caps = data["capex"]["capacity_options_jet_tpy"]
    capex_vals = data["capex"]["capex0_RMB_for_option"]
    if len(caps) != len(capex_vals):
        raise ValueError("capex.capacity_options_jet_tpy and capex.capex0_RMB_for_option length mismatch")
    if any(float(c) < 0 for c in caps):
        raise ValueError("capacity_options_jet_tpy contains negative values")
    if all(float(c) == 0 for c in caps):
        raise ValueError("All capacity options are 0; model becomes degenerate")

    # 概率检查
    psum = sum(float(sc.get("probability", 0.0)) for sc in data["scenarios"])
    if abs(psum - 1.0) > 1e-6:
        raise ValueError(f"Scenario probabilities must sum to 1. Got {psum}")

    # 生物质圈层检查
    if len(data["biomass_rings"]) == 0:
        raise ValueError("biomass_rings is empty")
    for r in data["biomass_rings"]:
        if float(r.get("max_supply_tpy", 0.0)) < 0:
            raise ValueError("Found negative ring max_supply_tpy")
        if float(r.get("distance_km", 0.0)) < 0:
            raise ValueError("Found negative ring distance_km")

    # process 参数范围检查（只做最关键的）
    proc = data["process"]
    xin = float(proc["moisture_in_wb"])
    xout = float(proc["moisture_out_wb"])
    if not (0.0 <= xin < 1.0) or not (0.0 <= xout < 1.0):
        raise ValueError("moisture_in_wb/moisture_out_wb must be in [0,1)")
    if xout > xin:
        # 可以允许（如果你定义的是“出料含水更高”不合理），这里强提醒
        raise ValueError("moisture_out_wb > moisture_in_wb: drying increases moisture? Check data.")

    econ = data["economic"]
    if float(econ["discount_rate"]) < 0:
        raise ValueError("discount_rate must be >= 0")
    if float(econ["degradation_rate"]) < 0:
        raise ValueError("degradation_rate must be >= 0")

    # 场景 split 归一化检查（只提醒：我们会在计算时强制归一化）
    for sc in data["scenarios"]:
        sp = sc.get("product_split", {})
        ssum = sum(float(sp.get(k, 0.0)) for k in ("kerosene_frac", "diesel_frac", "naphtha_frac", "wax_frac"))
        if abs(ssum - 1.0) > 1e-3:
            # 不直接报错，避免你历史数据跑不了；但会在计算中自动归一化
            pass


@dataclass
class RowOut:
    scenario: str
    probability: float
    build: int
    cap_jet_tpy: float
    biomass_in_ar_tpy: float

    # Key intermediate requirements
    dry_solids_tpy: float
    o2_required_tpy: float
    h2_makeup_tpy: float

    # Syngas + conversion
    syngas_CO_Nm3_y: float
    syngas_H2_Nm3_y: float
    CO_used_kmol_y: float

    # Products (FINAL products; jet includes wax-upgrading increment)
    jet_tpy: float
    diesel_tpy: float
    naphtha_tpy: float
    wax_sold_tpy: float
    wax_upgraded_tpy: float

    # Economics (annual)
    revenue_RMB_y: float
    opex_RMB_y: float
    annual_profit_RMB_y: float


def write_csv(rows: List[RowOut], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()) if rows else ["empty"])
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))


# -----------------------------
# Deterministic evaluation (no solver)
# -----------------------------

def _capex_from_choice(data: Dict[str, Any], cap_jet: float, build: int, wax_upgrade: int) -> float:
    """离散 CAPEX：严格匹配 option，否则取最近。"""
    if build <= 0:
        return 0.0
    caps = data["capex"]["capacity_options_jet_tpy"]
    capex_vals = data["capex"]["capex0_RMB_for_option"]
    try:
        idx = list(map(float, caps)).index(float(cap_jet))
    except ValueError:
        idx = min(range(len(caps)), key=lambda k: abs(float(caps[k]) - float(cap_jet)))
    base_capex = float(capex_vals[idx])
    extra_up = float(data["capex"].get("wax_upgrading_capex0_RMB", 0.0)) if wax_upgrade else 0.0
    return base_capex + extra_up


def eval_fixed(data: Dict[str, Any], bio_total_override: Optional[float] = None) -> Dict[str, Any]:
    """
    无求解器模式：用 model.fixed 定死 build/cap/wax_upgrade/biomass_in，然后把每个情景算一遍。

    """
    model = data["model"]
    fixed = model.get("fixed", {})
    build = int(fixed.get("build", 1))
    cap_jet = float(fixed.get("cap_jet_tpy", 30000))
    wax_upgrade = int(fixed.get("wax_upgrade", 1))

    # biomass rings
    rings = data["biomass_rings"]
    ring_supplies = [float(r["max_supply_tpy"]) for r in rings]
    if bio_total_override is not None:
        s0 = sum(ring_supplies) if sum(ring_supplies) > 0 else 1.0
        scale = float(bio_total_override) / s0
        ring_supplies = [x * scale for x in ring_supplies]

    biomass_in = float(fixed.get("biomass_in_ar_tpy", sum(ring_supplies)))

    # simplest allocation: nearest rings first
    B_r: List[float] = []
    remaining = biomass_in
    for sup in ring_supplies:
        take = min(sup, remaining)
        B_r.append(take)
        remaining -= take
    if remaining > 1e-9:
        raise RuntimeError("Fixed biomass_in exceeds ring supplies. Increase bio-total or ring max supplies.")

    out_rows: List[RowOut] = []
    for sc in data["scenarios"]:
        out_rows.append(_scenario_calcs_fixed_like_milp(data, sc, build, cap_jet, wax_upgrade, biomass_in, B_r))

    econ = data["economic"]
    i = float(econ["discount_rate"])
    delta = float(econ["degradation_rate"])
    pv8 = pv_factor(i, delta, 8)
    pv20 = pv_factor(i, delta, 20)

    capex0 = _capex_from_choice(data, cap_jet, build, wax_upgrade)
    exp_profit = sum(r.probability * r.annual_profit_RMB_y for r in out_rows)

    npv8 = -capex0 + pv8 * exp_profit
    npv20 = -capex0 + pv20 * exp_profit

    return {
        "status": "EVALUATED_WITHOUT_SOLVER",
        "build": build,
        "cap_jet_tpy": cap_jet,
        "wax_upgrade": wax_upgrade,
        "biomass_in_ar_tpy": biomass_in,
        "capex0_RMB": capex0,
        "expected_annual_profit_RMB_y": exp_profit,
        "PVF_8": pv8,
        "PVF_20": pv20,
        "NPV_8_RMB": npv8,
        "NPV_20_RMB": npv20,
        "rows": [asdict(r) for r in out_rows],
    }


def _scenario_calcs_fixed_like_milp(
    data: Dict[str, Any],
    sc: Dict[str, Any],
    build: int,
    cap_jet: float,
    wax_upgrade: int,
    biomass_in_ar_tpy: float,
    B_r_list: List[float],
) -> RowOut:
    """
    deterministic 的情景计算：尽量与 MILP 的“口径”一致（尤其是产能约束与需量电）。
    研究级改进：capacity 约束默认施加在“最终 jet（含升级）”上。
    """
    proc = data["process"]
    feed = data["feedstock"]
    econ = data["economic"]

    # ---- Biomass cost (ring-based)
    bailing = float(econ["biomass_base_cost_RMB_per_t"])
    transport = float(econ["biomass_transport_RMB_per_t_km"])
    rings = data["biomass_rings"]

    C_bio = 0.0
    for Br, ring in zip(B_r_list, rings):
        d = float(ring["distance_km"])
        C_bio += (bailing + transport * d) * Br

    # ---- Moisture balance
    x_in = float(proc["moisture_in_wb"])
    x_out = float(proc["moisture_out_wb"])
    dry_solids_tpy = biomass_in_ar_tpy * (1.0 - x_in)
    after_dry_ar_tpy = safe_div(dry_solids_tpy, (1.0 - x_out))
    water_removed_tpy = max(0.0, biomass_in_ar_tpy - after_dry_ar_tpy)

    # ---- Char (dry and daf)
    ash_pct = float(feed["ultimate_dry_wt_percent"].get("Ash", 0.0)) / 100.0
    y_char = float(sc["char_yield_dry_fraction"])
    char_dry_tpy = y_char * dry_solids_tpy
    char_daf_tpy = max(0.0, char_dry_tpy - ash_pct * dry_solids_tpy)

    # ---- Syngas yields (Nm3/y)
    Y_H2 = float(sc["H2_yield_Nm3_per_kg_char_daf"])
    Y_CO = float(sc["CO_yield_Nm3_per_kg_char_daf"])
    syngas_H2 = Y_H2 * char_daf_tpy * 1000.0
    syngas_CO = Y_CO * char_daf_tpy * 1000.0

    # ---- O2 requirement (t/y)
    o2_ratio = float(proc["entrained_flow"]["O2_to_dry_biomass_mass_ratio"])
    o2_required_tpy = o2_ratio * dry_solids_tpy

    # ---- Convert syngas to kmol
    Nm3_per_kmol = float(proc["constants"]["Nm3_per_kmol"])
    H2_kmol = syngas_H2 / Nm3_per_kmol
    CO_kmol = syngas_CO / Nm3_per_kmol

    # ---- Optional: linear WGS module（默认关闭；研究生升级：给出“可扩展模块”）
    enable_wgs = bool(proc.get("wgs", {}).get("enable_WGS", False))
    CO_shift = 0.0
    if enable_wgs and build:
        # 这里为了 deterministic 简单起见：shift 到满足 H2/CO 所需为止（不一定最优，只是演示）
        R_ft = float(proc["ft"]["target_H2_CO"])
        # 如果 H2 不够，可 shift 一部分 CO 生成 H2：每 1 kmol CO_shift -> H2 +1
        # shift 后：H2_eff = H2 + CO_shift; CO_eff = CO - CO_shift
        # 要满足 H2_eff >= R * CO_used
        # deterministic 我们先不深入；MILP 中会优化决定
        pass

    # ---- FT conversion and product split
    eta_liq = float(sc["ft_liquids_t_per_kmol_CO"])

    split = normalize_split(sc.get("product_split", {}), ("kerosene_frac", "diesel_frac", "naphtha_frac", "wax_frac"))
    f_ker = split["kerosene_frac"]
    f_dis = split["diesel_frac"]
    f_nap = split["naphtha_frac"]
    f_wax = split["wax_frac"]

    # ---- If build=0, no production
    if build <= 0:
        CO_used_kmol = 0.0
        liquids_tpy = 0.0
        wax = 0.0
        wax_upgraded = 0.0
        wax_sold = 0.0
        jet = diesel = naphtha = 0.0
        h2_makeup_tpy = 0.0
    else:
        # 研究级改进：capacity 默认约束 final jet_total <= cap_jet
        # jet_total = f_ker*liquids + gamma_jet*wax_upgraded
        # 但 wax_upgraded 本身取决于 wax（=f_wax*liquids）与是否升级
        gam = sc["wax_upgrading_yields"]
        gamma_jet = float(gam["jet_per_t_wax"])
        gamma_dis = float(gam["diesel_per_t_wax"])
        gamma_nap = float(gam["naphtha_per_t_wax"])

        # 如果升级开启，我们采用“升级全部蜡”（deterministic 规则），否则不升级
        # wax = f_wax*liquids
        # jet_total = f_ker*liquids + gamma_jet*(wax if upgrade else 0)
        # => jet_total = (f_ker + gamma_jet*f_wax)*liquids  (if upgrade)
        # => jet_total = f_ker*liquids (if not upgrade)
        jet_coeff = f_ker + (gamma_jet * f_wax if wax_upgrade else 0.0)

        # CO_used limited by available CO and by final jet capacity
        # liquids = eta*CO_used => jet_total = jet_coeff*eta*CO_used <= cap
        CO_used_kmol = min(CO_kmol, safe_div(cap_jet, (jet_coeff * eta_liq)))

        R_ft = float(proc["ft"]["target_H2_CO"])
        H2_makeup_kmol = max(0.0, R_ft * CO_used_kmol - H2_kmol)
        h2_makeup_tpy = H2_makeup_kmol * 2.0 / 1000.0

        liquids_tpy = eta_liq * CO_used_kmol

        jet_base = f_ker * liquids_tpy
        diesel_base = f_dis * liquids_tpy
        naphtha_base = f_nap * liquids_tpy
        wax = f_wax * liquids_tpy

        wax_upgraded = wax if wax_upgrade else 0.0
        wax_sold = wax - wax_upgraded

        jet = jet_base + gamma_jet * wax_upgraded
        diesel = diesel_base + gamma_dis * wax_upgraded
        naphtha = naphtha_base + gamma_nap * wax_upgraded

        # 最终 jet 再做一次保险裁剪（避免浮点误差 > cap）
        jet = min(jet, cap_jet)

    # ---- Utilities / OPEX
    q_dry = float(sc["drying_heat_MJ_per_kg_water_removed"])
    drying_heat_MJ_y = q_dry * water_removed_tpy * 1000.0

    LHV_dry = float(feed["LHV_dry_MJ_per_kg"])
    thermal_input_kWhth = dry_solids_tpy * 1000.0 * LHV_dry / 3.6

    e_grind = float(sc.get("grinding_kWe_per_kWth", 0.0))
    e_feedpress = float(sc.get("feed_press_kWe_per_kWth", 0.0))
    elec_grind = e_grind * thermal_input_kWhth
    elec_feed = e_feedpress * thermal_input_kWhth

    asu_kWh_tO2 = float(sc.get("ASU_kWh_per_t_O2", 0.0))
    elec_asu = o2_required_tpy * asu_kWh_tO2

    ft_kWh_per_t_liq = float(proc["ft"]["electricity_kWh_per_t_liquids"])
    up_kWh_per_t_wax = float(proc["ft"]["upgrading_electricity_kWh_per_t_wax"])
    elec_ft = ft_kWh_per_t_liq * liquids_tpy
    elec_up = up_kWh_per_t_wax * wax_upgraded

    E_kWh = elec_grind + elec_feed + elec_asu + elec_ft + elec_up

    # 两部制电价：研究级改进（峰值近似）
    p_kWh = float(econ["electricity_energy_RMB_per_kWh"])
    p_dem = float(econ["electricity_demand_RMB_per_kW_month"])
    H_op = float(econ["electricity_operating_hours_per_year"])

    # 新增：load_factor（平均功率/峰值功率），默认 1.0 退化回你原来的做法
    load_factor = float(econ.get("peak_load_factor", 1.0))
    load_factor = max(min(load_factor, 1.0), 1e-3)

    P_peak = E_kWh / (H_op * load_factor) if H_op > 0 else 0.0
    C_elec = p_kWh * E_kWh + 12.0 * p_dem * P_peak

    p_heat = float(econ.get("heat_RMB_per_MJ", 0.0))
    C_heat = p_heat * drying_heat_MJ_y

    p_H2 = float(econ["H2_RMB_per_t"])
    C_H2 = p_H2 * h2_makeup_tpy

    C_FOM = float(econ.get("fixed_OM_RMB_per_y", 0.0)) * (1 if build else 0)

    P_jet = float(econ["product_prices_RMB_per_t"]["jet"])
    P_dis = float(econ["product_prices_RMB_per_t"].get("diesel", 0.0))
    P_nap = float(econ["product_prices_RMB_per_t"].get("naphtha", 0.0))
    P_wax = float(econ["product_prices_RMB_per_t"].get("wax", 0.0))

    revenue = P_jet * jet + P_dis * diesel + P_nap * naphtha + P_wax * wax_sold
    opex = C_bio + C_elec + C_heat + C_H2 + C_FOM
    profit = revenue - opex

    return RowOut(
        scenario=str(sc["name"]),
        probability=float(sc["probability"]),
        build=int(build),
        cap_jet_tpy=float(cap_jet),
        biomass_in_ar_tpy=float(biomass_in_ar_tpy),
        dry_solids_tpy=float(dry_solids_tpy),
        o2_required_tpy=float(o2_required_tpy),
        h2_makeup_tpy=float(h2_makeup_tpy),
        syngas_CO_Nm3_y=float(syngas_CO),
        syngas_H2_Nm3_y=float(syngas_H2),
        CO_used_kmol_y=float(CO_used_kmol),
        jet_tpy=float(jet),
        diesel_tpy=float(diesel),
        naphtha_tpy=float(naphtha),
        wax_sold_tpy=float(wax_sold),
        wax_upgraded_tpy=float(wax_upgraded),
        revenue_RMB_y=float(revenue),
        opex_RMB_y=float(opex),
        annual_profit_RMB_y=float(profit),
    )


# -----------------------------
# MILP with PuLP (research-grade upgrades)
# -----------------------------

def choose_solver(name: str, msg: bool = False):
    """
    研究级：允许切换求解器。
    - CBC: 默认，开源
    - HiGHS: 若你的 pulp 版本支持（部分环境）
    - GUROBI: 需要安装与 license
    """
    name = (name or "cbc").strip().lower()
    if name == "cbc":
        return pulp.PULP_CBC_CMD(msg=msg)
    if name == "gurobi":
        # 你装了 gurobi 并配置好 license 才能用
        try:
            return pulp.GUROBI_CMD(msg=msg)
        except Exception as e:
            raise RuntimeError(f"GUROBI solver not available in this environment: {e}")
    if name == "highs":
        # 不同 pulp 版本接口不同，这里做 best-effort
        if hasattr(pulp, "HiGHS_CMD"):
            return pulp.HiGHS_CMD(msg=msg)  # type: ignore
        raise RuntimeError("HiGHS_CMD not found in this PuLP version")
    raise ValueError(f"Unknown solver: {name}")


def solve_milp(
    data: Dict[str, Any],
    bio_total_override: Optional[float] = None,
    force_build: bool = False,
    relax_gate: bool = False,
    fix_upgrade: Optional[int] = None,  # 0/1 to force, None = let MILP decide
    solver_name: str = "cbc",
    solver_msg: bool = False,
    # --- risk / research options ---
    risk_alpha: Optional[float] = None,   # e.g. 0.9
    risk_lambda: float = 0.0,             # 0 -> risk neutral (your original)
) -> Dict[str, Any]:

    if not HAS_PULP:
        raise RuntimeError("PuLP not installed. Run with --no-solver for deterministic evaluation.")

    proc = data["process"]
    feed = data["feedstock"]
    econ = data["economic"]
    caps = list(map(float, data["capex"]["capacity_options_jet_tpy"]))
    capex_vals = list(map(float, data["capex"]["capex0_RMB_for_option"]))

    # Ring supplies
    rings = data["biomass_rings"]
    ring_supplies = [float(r["max_supply_tpy"]) for r in rings]
    if bio_total_override is not None:
        s0 = sum(ring_supplies) if sum(ring_supplies) > 0 else 1.0
        scale = float(bio_total_override) / s0
        ring_supplies = [x * scale for x in ring_supplies]

    # PV factors
    i = float(econ["discount_rate"])
    delta = float(econ["degradation_rate"])
    PVF8 = pv_factor(i, delta, 8)
    PVF20 = pv_factor(i, delta, 20)

    prob = pulp.LpProblem("Biomass_to_Jet_FT_basic_v2", pulp.LpMaximize)

    # ---------------- First-stage ----------------
    y = pulp.LpVariable("build", lowBound=0, upBound=1, cat=pulp.LpBinary)

    # Discrete capacity option (exactly one if built)
    x = [pulp.LpVariable(f"cap_choice_{int(c)}", lowBound=0, upBound=1, cat=pulp.LpBinary) for c in caps]
    prob += pulp.lpSum(x) == y, "choose_one_capacity_if_built"

    cap_jet = pulp.LpVariable("cap_jet_tpy", lowBound=0)
    prob += cap_jet == pulp.lpSum(c * xj for c, xj in zip(caps, x)), "cap_jet_def"

    # Wax upgrading choice
    u = pulp.LpVariable("wax_upgrade", lowBound=0, upBound=1, cat=pulp.LpBinary)
    prob += u <= y, "upgrade_only_if_built"

    if force_build:
        prob += y == 1, "force_build"
    if fix_upgrade is not None:
        prob += u == int(fix_upgrade), "force_upgrade_choice"

    # Biomass procurement by ring
    B_r = [pulp.LpVariable(f"bio_ring_{k}_tpy", lowBound=0) for k in range(len(rings))]
    for k, sup in enumerate(ring_supplies):
        prob += B_r[k] <= sup * y, f"ring_supply_{k}"

    B = pulp.LpVariable("biomass_in_ar_tpy", lowBound=0)
    prob += B == pulp.lpSum(B_r), "biomass_total"

    # 研究级：避免 y=1 但 B=0 的“空建厂”退化解（可选）
    B_min = float(data.get("model", {}).get("min_biomass_if_built_tpy", 0.0))
    if B_min > 0:
        prob += B >= B_min * y, "min_biomass_if_built"

    # Electricity demand (contracted peak kW)
    P_peak = pulp.LpVariable("P_peak_kW", lowBound=0)

    # CAPEX0
    base_capex0 = pulp.lpSum(v * xj for v, xj in zip(capex_vals, x))
    up_capex0 = float(data["capex"].get("wax_upgrading_capex0_RMB", 0.0)) * u

    capex0 = pulp.LpVariable("capex0_RMB", lowBound=0)
    prob += capex0 == base_capex0 + up_capex0, "capex0_def"

    # ---------------- constants ----------------
    x_in = float(proc["moisture_in_wb"])
    x_out = float(proc["moisture_out_wb"])
    ash_pct = float(feed["ultimate_dry_wt_percent"].get("Ash", 0.0)) / 100.0
    LHV_dry = float(feed["LHV_dry_MJ_per_kg"])
    o2_ratio = float(proc["entrained_flow"]["O2_to_dry_biomass_mass_ratio"])
    Nm3_per_kmol = float(proc["constants"]["Nm3_per_kmol"])
    R_ft = float(proc["ft"]["target_H2_CO"])

    # economics constants
    bailing = float(econ["biomass_base_cost_RMB_per_t"])
    transport = float(econ["biomass_transport_RMB_per_t_km"])
    p_kWh = float(econ["electricity_energy_RMB_per_kWh"])
    p_dem = float(econ["electricity_demand_RMB_per_kW_month"])
    H_op = float(econ["electricity_operating_hours_per_year"])
    p_heat = float(econ.get("heat_RMB_per_MJ", 0.0))
    p_H2 = float(econ["H2_RMB_per_t"])
    C_FOM = float(econ.get("fixed_OM_RMB_per_y", 0.0)) * y

    # 研究级：需量电峰值近似参数（默认1.0等价原版）
    load_factor = float(econ.get("peak_load_factor", 1.0))
    load_factor = max(min(load_factor, 1.0), 1e-3)

    P_jet = float(econ["product_prices_RMB_per_t"]["jet"])
    P_dis = float(econ["product_prices_RMB_per_t"].get("diesel", 0.0))
    P_nap = float(econ["product_prices_RMB_per_t"].get("naphtha", 0.0))
    P_wax = float(econ["product_prices_RMB_per_t"].get("wax", 0.0))

    # Biomass cost (scenario-independent)
    C_bio = 0
    for Br, ring in zip(B_r, rings):
        d = float(ring["distance_km"])
        C_bio += (bailing + transport * d) * Br

    # ---------------- Second-stage per scenario ----------------
    scenario_profit: Dict[str, Tuple[float, Any]] = {}
    scenario_rows: Dict[str, Dict[str, Any]] = {}

    # 可选：线性 WGS（默认关闭）
    enable_wgs = bool(proc.get("wgs", {}).get("enable_WGS", False))
    wgs_cost_per_kmol = float(proc.get("wgs", {}).get("steam_cost_RMB_per_kmol_CO_shift", 0.0))

    for sc in data["scenarios"]:
        name = str(sc["name"])
        prob_s = float(sc["probability"])

        # ------- variables -------
        dry_solids = pulp.LpVariable(f"dry_solids_tpy__{name}", lowBound=0)
        after_dry = pulp.LpVariable(f"after_dry_ar_tpy__{name}", lowBound=0)
        water_removed = pulp.LpVariable(f"water_removed_tpy__{name}", lowBound=0)

        char_dry = pulp.LpVariable(f"char_dry_tpy__{name}", lowBound=0)
        char_daf = pulp.LpVariable(f"char_daf_tpy__{name}", lowBound=0)

        H2_Nm3 = pulp.LpVariable(f"H2_Nm3_y__{name}", lowBound=0)
        CO_Nm3 = pulp.LpVariable(f"CO_Nm3_y__{name}", lowBound=0)

        CO_kmol = pulp.LpVariable(f"CO_kmol_y__{name}", lowBound=0)
        H2_kmol = pulp.LpVariable(f"H2_kmol_y__{name}", lowBound=0)

        # Optional WGS shift
        CO_shift = pulp.LpVariable(f"CO_shift_kmol_y__{name}", lowBound=0) if enable_wgs else None

        CO_used = pulp.LpVariable(f"CO_used_kmol_y__{name}", lowBound=0)
        H2_makeup_kmol = pulp.LpVariable(f"H2_makeup_kmol_y__{name}", lowBound=0)
        H2_makeup_tpy = pulp.LpVariable(f"H2_makeup_tpy__{name}", lowBound=0)

        liquids = pulp.LpVariable(f"liquids_tpy__{name}", lowBound=0)

        # base products from liquids
        jet_base = pulp.LpVariable(f"jet_base_tpy__{name}", lowBound=0)
        diesel_base = pulp.LpVariable(f"diesel_base_tpy__{name}", lowBound=0)
        naphtha_base = pulp.LpVariable(f"naphtha_base_tpy__{name}", lowBound=0)
        wax = pulp.LpVariable(f"wax_tpy__{name}", lowBound=0)

        wax_up = pulp.LpVariable(f"wax_upgraded_tpy__{name}", lowBound=0)
        wax_sell = pulp.LpVariable(f"wax_sold_tpy__{name}", lowBound=0)

        # final products (include upgrade increment)
        jet_total = pulp.LpVariable(f"jet_total_tpy__{name}", lowBound=0)
        diesel_total = pulp.LpVariable(f"diesel_total_tpy__{name}", lowBound=0)
        naphtha_total = pulp.LpVariable(f"naphtha_total_tpy__{name}", lowBound=0)

        # ------- moisture balances -------
        prob += dry_solids == B * (1.0 - x_in), f"dry_solids_def__{name}"
        prob += after_dry == dry_solids / (1.0 - x_out), f"after_dry_def__{name}"
        prob += water_removed == B - after_dry, f"water_removed_def__{name}"

        # ------- char -------
        y_char = float(sc["char_yield_dry_fraction"])
        prob += char_dry == y_char * dry_solids, f"char_dry_def__{name}"
        prob += char_daf == char_dry - ash_pct * dry_solids, f"char_daf_def__{name}"
        prob += char_daf >= 0, f"char_daf_nonneg__{name}"

        # ------- syngas yields -------
        Y_H2 = float(sc["H2_yield_Nm3_per_kg_char_daf"])
        Y_CO = float(sc["CO_yield_Nm3_per_kg_char_daf"])
        prob += H2_Nm3 == Y_H2 * char_daf * 1000.0, f"H2_yield__{name}"
        prob += CO_Nm3 == Y_CO * char_daf * 1000.0, f"CO_yield__{name}"

        # convert to kmol
        prob += H2_kmol == H2_Nm3 / Nm3_per_kmol, f"H2_kmol_def__{name}"
        prob += CO_kmol == CO_Nm3 / Nm3_per_kmol, f"CO_kmol_def__{name}"

        # ------- optional WGS (linear) -------
        if enable_wgs:
            # 物理：0 <= CO_shift <= CO_kmol
            prob += CO_shift <= CO_kmol, f"CO_shift_le_CO__{name}"  # type: ignore
            # shift 后有效合成气：H2_eff = H2 + CO_shift, CO_eff = CO - CO_shift
            H2_eff = H2_kmol + CO_shift  # type: ignore
            CO_eff = CO_kmol - CO_shift  # type: ignore
        else:
            H2_eff = H2_kmol
            CO_eff = CO_kmol

        # CO used cannot exceed available effective CO
        prob += CO_used <= CO_eff, f"CO_used_le_COeff__{name}"

        # H2 make-up meets target ratio (no RWGS; WGS only affects availability)
        prob += H2_makeup_kmol >= R_ft * CO_used - H2_eff, f"H2_makeup_ratio__{name}"
        prob += H2_makeup_kmol >= 0, f"H2_makeup_nonneg__{name}"
        prob += H2_makeup_tpy == H2_makeup_kmol * 2.0 / 1000.0, f"H2_makeup_mass__{name}"

        # liquids from CO used
        eta_liq = float(sc["ft_liquids_t_per_kmol_CO"])
        prob += liquids == eta_liq * CO_used, f"liquids_def__{name}"

        # product split (归一化后再进模型，避免 split 和 !=1 的数值问题)
        split = normalize_split(sc.get("product_split", {}), ("kerosene_frac", "diesel_frac", "naphtha_frac", "wax_frac"))
        f_ker = float(split["kerosene_frac"])
        f_dis = float(split["diesel_frac"])
        f_nap = float(split["naphtha_frac"])
        f_wax = float(split["wax_frac"])

        prob += jet_base == f_ker * liquids, f"jet_split__{name}"
        prob += diesel_base == f_dis * liquids, f"diesel_split__{name}"
        prob += naphtha_base == f_nap * liquids, f"naphtha_split__{name}"
        prob += wax == f_wax * liquids, f"wax_split__{name}"

        # wax upgrading logic: upgraded wax <= wax and only if u=1
        M_wax = float(data["model"].get("bigM_wax", 1e7))
        prob += wax_up <= wax, f"wax_up_le_wax__{name}"
        prob += wax_up <= M_wax * u, f"wax_up_only_if_upgrade__{name}"
        prob += wax_sell == wax - wax_up, f"wax_sell_def__{name}"

        # 研究级：可选升级能力上限（更像真实装置）
        up_cap = data.get("model", {}).get("wax_upgrading_max_tpy", None)
        if up_cap is not None:
            prob += wax_up <= float(up_cap) * u, f"wax_up_capacity__{name}"

        # upgrading yields (incremental)
        gam = sc["wax_upgrading_yields"]
        gamma_jet = float(gam["jet_per_t_wax"])
        gamma_dis = float(gam["diesel_per_t_wax"])
        gamma_nap = float(gam["naphtha_per_t_wax"])

        prob += jet_total == jet_base + gamma_jet * wax_up, f"jet_total_def__{name}"
        prob += diesel_total == diesel_base + gamma_dis * wax_up, f"diesel_total_def__{name}"
        prob += naphtha_total == naphtha_base + gamma_nap * wax_up, f"naphtha_total_def__{name}"

        # ✅ 研究级关键修正：产能约束作用于“最终 jet_total”
        prob += jet_total <= cap_jet, f"jet_cap_final__{name}"

        # ------- utilities -------
        q_dry = float(sc["drying_heat_MJ_per_kg_water_removed"])
        drying_heat_MJ = q_dry * water_removed * 1000.0

        thermal_input_kWhth = dry_solids * 1000.0 * LHV_dry / 3.6
        e_grind = float(sc.get("grinding_kWe_per_kWth", 0.0))
        e_feed = float(sc.get("feed_press_kWe_per_kWth", 0.0))
        elec_grind = e_grind * thermal_input_kWhth
        elec_feed = e_feed * thermal_input_kWhth

        # ASU electricity
        o2_required = o2_ratio * dry_solids
        asu_kWh_tO2 = float(sc.get("ASU_kWh_per_t_O2", 0.0))
        elec_asu = o2_required * asu_kWh_tO2

        # FT and upgrading electricity
        ft_kWh_per_t_liq = float(proc["ft"]["electricity_kWh_per_t_liquids"])
        up_kWh_per_t_wax = float(proc["ft"]["upgrading_electricity_kWh_per_t_wax"])
        elec_ft = ft_kWh_per_t_liq * liquids
        elec_up = up_kWh_per_t_wax * wax_up

        E_kWh = elec_grind + elec_feed + elec_asu + elec_ft + elec_up

        # ✅ 研究级：峰值近似（load_factor）
        #   P_peak >= E / (H_op * load_factor)
        prob += E_kWh <= (H_op * load_factor) * P_peak, f"peak_demand_cover__{name}"

        # ------- opex & revenue -------
        C_elec = p_kWh * E_kWh + 12.0 * p_dem * P_peak
        C_heat = p_heat * drying_heat_MJ
        C_H2 = p_H2 * H2_makeup_tpy

        C_WGS = (wgs_cost_per_kmol * CO_shift) if enable_wgs else 0  # type: ignore

        OPEX = C_bio + C_elec + C_heat + C_H2 + C_FOM + C_WGS
        REV = P_jet * jet_total + P_dis * diesel_total + P_nap * naphtha_total + P_wax * wax_sell
        PROFIT = REV - OPEX

        scenario_profit[name] = (prob_s, PROFIT)
        scenario_rows[name] = {
            "jet_total": jet_total,
            "diesel_total": diesel_total,
            "naphtha_total": naphtha_total,
            "wax_sell": wax_sell,
            "wax_up": wax_up,
            "o2_required": o2_required,
            "H2_makeup_tpy": H2_makeup_tpy,
            "CO_used": CO_used,
            "H2_Nm3": H2_Nm3,
            "CO_Nm3": CO_Nm3,
            "E_kWh": E_kWh,
            "dry_solids": dry_solids,
            "REV": REV,
            "OPEX": OPEX,
            "PROFIT": PROFIT,
        }

    # Expected annual profit
    exp_profit = pulp.lpSum(prob_s * prof for (prob_s, prof) in scenario_profit.values())

    # NPV expressions (linear)
    NPV8 = pulp.LpVariable("NPV_8_RMB", lowBound=None)
    NPV20 = pulp.LpVariable("NPV_20_RMB", lowBound=None)
    prob += NPV8 == -capex0 + PVF8 * exp_profit, "NPV8_def"
    prob += NPV20 == -capex0 + PVF20 * exp_profit, "NPV20_def"

    # Financing gate (policy)
    if not relax_gate:
        prob += NPV8 >= 0, "financing_gate_NPV8_nonnegative"

    # ---------------- Risk (CVaR) optional ----------------
    # 研究级升级：对随机收益做下行风险度量
    # 定义 lower-tail CVaR of profit: CVaR = z - 1/(1-a) * Σ π_s * η_s,  η_s >= z - Profit_s
    # risk_lambda>0 时：最大化  NPV20 - risk_lambda * PVF20 * (E[Profit] - CVaR)
    # 解释：惩罚“期望收益 - 下行尾部收益”，让解更稳健
    CVaR = None
    if risk_alpha is not None and risk_lambda > 0:
        a = float(risk_alpha)
        if not (0.0 < a < 1.0):
            raise ValueError("risk_alpha must be in (0,1)")
        z = pulp.LpVariable("VaR_profit_z", lowBound=None)
        eta = {name: pulp.LpVariable(f"shortfall_eta__{name}", lowBound=0) for name in scenario_profit.keys()}

        for name, (pi, prof) in scenario_profit.items():
            prob += eta[name] >= z - prof, f"cvar_shortfall__{name}"

        CVaR = pulp.LpVariable("CVaR_profit", lowBound=None)
        prob += CVaR == z - (1.0 / (1.0 - a)) * pulp.lpSum(
            float(scenario_profit[name][0]) * eta[name] for name in scenario_profit.keys()
        ), "CVaR_def"

        # 风险惩罚项：E[Profit]-CVaR >=0
        risk_penalty = exp_profit - CVaR
        # 目标：最大化 NPV20 - λ * PVF20 * (E[Profit]-CVaR)
        prob += NPV20 - float(risk_lambda) * PVF20 * risk_penalty, "risk_adjusted_objective"
    else:
        # Risk-neutral objective (your original)
        prob += NPV20, "maximize_NPV20"

    # ---------------- Solve ----------------
    solver = choose_solver(solver_name, msg=solver_msg)
    prob.solve(solver)

    status = pulp.LpStatus.get(prob.status, str(prob.status))

    def val(v):
        try:
            return float(pulp.value(v))
        except Exception:
            return float("nan")

    build_sol = int(round(val(y)))
    cap_sol = val(cap_jet)
    up_sol = int(round(val(u)))
    B_sol = val(B)
    capex0_sol = val(capex0)
    exp_profit_sol = val(exp_profit)
    npv8_sol = val(NPV8)
    npv20_sol = val(NPV20)

    # Scenario rows
    rows_out: List[RowOut] = []
    for sc in data["scenarios"]:
        name = str(sc["name"])
        vs = scenario_rows[name]
        rows_out.append(RowOut(
            scenario=name,
            probability=float(sc["probability"]),
            build=build_sol,
            cap_jet_tpy=cap_sol,
            biomass_in_ar_tpy=B_sol,
            dry_solids_tpy=val(vs["dry_solids"]),
            o2_required_tpy=val(vs["o2_required"]),
            h2_makeup_tpy=val(vs["H2_makeup_tpy"]),
            syngas_CO_Nm3_y=val(vs["CO_Nm3"]),
            syngas_H2_Nm3_y=val(vs["H2_Nm3"]),
            CO_used_kmol_y=val(vs["CO_used"]),
            jet_tpy=val(vs["jet_total"]),
            diesel_tpy=val(vs["diesel_total"]),
            naphtha_tpy=val(vs["naphtha_total"]),
            wax_sold_tpy=val(vs["wax_sell"]),
            wax_upgraded_tpy=val(vs["wax_up"]),
            revenue_RMB_y=val(vs["REV"]),
            opex_RMB_y=val(vs["OPEX"]),
            annual_profit_RMB_y=val(vs["PROFIT"]),
        ))

    res = {
        "status": status,
        "build": build_sol,
        "cap_jet_tpy": cap_sol,
        "wax_upgrade": up_sol,
        "biomass_in_ar_tpy": B_sol,
        "capex0_RMB": capex0_sol,
        "expected_annual_profit_RMB_y": exp_profit_sol,
        "PVF_8": PVF8,
        "PVF_20": PVF20,
        "NPV_8_RMB": npv8_sol,
        "NPV_20_RMB": npv20_sol,
        "rows": [asdict(r) for r in rows_out],
    }
    if CVaR is not None:
        res["CVaR_profit_RMB_y"] = val(CVaR)

    return res


# -----------------------------
# CLI main
# -----------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=str, default="biomass_to_jet_ft_basic_data.json", help="Input JSON data file")
    ap.add_argument("--out", type=str, default="ft_basic_results.csv", help="Output results CSV")
    ap.add_argument("--no-solver", action="store_true", help="Force deterministic evaluation (ignore PuLP)")
    ap.add_argument("--bio-total", type=float, default=None, help="Override total biomass availability (t/y)")
    ap.add_argument("--solver", type=str, default="cbc", help="MILP solver: cbc | highs | gurobi")
    ap.add_argument("--solver-msg", action="store_true", help="Show solver log")
    ap.add_argument("--report-counterfactual", action="store_true",
                    help="Also solve counterfactual cases (force build=1, relax NPV8 gate).")
    ap.add_argument("--report-both-upgrade", action="store_true",
                    help="When reporting counterfactual, also solve both wax_upgrade=0 and wax_upgrade=1.")
    # 风险选项（研究级）
    ap.add_argument("--risk-alpha", type=float, default=None, help="Enable CVaR with alpha (e.g., 0.9)")
    ap.add_argument("--risk-lambda", type=float, default=0.0,
                    help="Risk penalty weight (>0 enables risk-adjusted objective).")

    args = ap.parse_args(argv)

    data_path = resolve_relative_to_script(Path(args.data))
    if not data_path.exists():
        raise SystemExit(f"[ERROR] Data file not found: {data_path}")

    data = load_json(data_path)
    validate_data(data)

    def _write_one(res_obj: Dict[str, Any], out_path: Path) -> None:
        rows = [RowOut(**r) for r in res_obj["rows"]]  # type: ignore
        if rows:
            write_csv(rows, out_path)
            print(f"\nWrote results CSV to: {out_path}")

    out_base = resolve_relative_to_script(Path(args.out))

    def _suffix_path(sfx: str) -> Path:
        return out_base.with_name(out_base.stem + sfx + out_base.suffix)

    # ---------- policy solve ----------
    if args.no_solver or (not HAS_PULP):
        res_policy = eval_fixed(data, bio_total_override=args.bio_total)
    else:
        res_policy = solve_milp(
            data,
            bio_total_override=args.bio_total,
            solver_name=args.solver,
            solver_msg=args.solver_msg,
            risk_alpha=args.risk_alpha,
            risk_lambda=args.risk_lambda,
        )

    bundle: Dict[str, Any] = {"policy": res_policy}
    _write_one(res_policy, _suffix_path("_policy"))

    # ---------- counterfactual solves ----------
    if (not args.no_solver) and HAS_PULP and args.report_counterfactual:
        # 1) counterfactual best
        res_cf_best = solve_milp(
            data,
            bio_total_override=args.bio_total,
            force_build=True,
            relax_gate=True,
            fix_upgrade=None,
            solver_name=args.solver,
            solver_msg=args.solver_msg,
            risk_alpha=args.risk_alpha,
            risk_lambda=args.risk_lambda,
        )
        bundle["counterfactual_best"] = res_cf_best
        _write_one(res_cf_best, _suffix_path("_cf_best"))

        # 2) both wax options
        if args.report_both_upgrade:
            res_cf_no = solve_milp(
                data,
                bio_total_override=args.bio_total,
                force_build=True,
                relax_gate=True,
                fix_upgrade=0,
                solver_name=args.solver,
                solver_msg=args.solver_msg,
                risk_alpha=args.risk_alpha,
                risk_lambda=args.risk_lambda,
            )
            res_cf_up = solve_milp(
                data,
                bio_total_override=args.bio_total,
                force_build=True,
                relax_gate=True,
                fix_upgrade=1,
                solver_name=args.solver,
                solver_msg=args.solver_msg,
                risk_alpha=args.risk_alpha,
                risk_lambda=args.risk_lambda,
            )
            bundle["counterfactual_wax_upgrade_0"] = res_cf_no
            bundle["counterfactual_wax_upgrade_1"] = res_cf_up
            _write_one(res_cf_no, _suffix_path("_cf_upgrade0"))
            _write_one(res_cf_up, _suffix_path("_cf_upgrade1"))

    print(json.dumps(bundle, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
