#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 12:32:04 2025

@author: hashima
"""

import sys
from juliacall import Main as jl

try:
    # Juliaモジュールをロード
    # ファイルパスを適切に指定してください
    jl.include("power_flow.jl")
    jl.seval("using .power_flow")
    VSE = jl.power_flow
    # メインモジュールから関数にアクセス
    solve_power_flow = VSE.solve_power_flow

except Exception as e:
    print(f"Failed to load Julia module: {e}")
    print("Please ensure 'PowerFlow.jl' is in the same directory.")
    print("You may need to install NLsolve in Julia:")
    print("from juliacall import Main as jl")
    print("jl.seval(\"using Pkg; Pkg.add('NLsolve')\")")
    sys.exit(1)


# 表6.8  の Case 1 のデータを設定
# P, Q は MW/Mvar 単位
case3_data = {
    "V1": 1.08,
    "P2": 200.0,
    "V2": 1.01,
    "P3": -1500.0,
    "Q3": -150.0,
    "Bs3": 6.822
}

print("Solving for Case 1...")
result = solve_power_flow(case3_data)

# 結果の表示
if result["status"] == "Converged":
    print("Power Flow Converged.")
    print("-" * 30)
    print(f"LTC Tap Ratio (r): {result['LTC_tap_r']:.6f}")
    print("-" * 30)
    print(f"Bus 1 (Slack): V = {result['V1']:.4f} pu, Angle(θth base) = {result['θ1_deg']:.4f} rad")
    print(f"           P = {result['P1_MW']:.2f} MW, Q = {result['Q1_Mvar']:.2f} Mvar")
    print(f"Bus 2 (PV):    V = {result['V2']:.4f} pu, Angle(θth base) = {result['θ2_deg']:.4f} rad")
    print(f"           P = {case3_data['P2']:.2f} MW, Q = {result['Q2_Mvar']:.2f} Mvar")
    print(f"Bus 3 (PQ/LTC):V = {result['V3']:.4f} pu, Angle(θth base) = {result['θ3_deg']:.4f} rad")
    print(f"           P = {case3_data['P3']:.2f} MW, Q = {case3_data['Q3']:.2f} Mvar")
    print(f"Bus 4 (PQ):    V = {result['V4']:.4f} pu, Angle(θth base) = {result['θ4_deg']:.4f} rad")
    print(f"           P = 0.00 MW, Q = 0.00 Mvar")
    print(f"           Eth = {result['Eth']:.4f}, θth = {result['θth']:.4}")
    print(f"           δ = {result['δ']:.4f}, V2o = {result['V2o']:.4}")
else:
    print(f"Solver failed: {result['error']}")
    print(f"Details: {result['details']}")

# # --- Case 2 や Case 3 も同様に実行可能 ---
# #  (Table 6.8)
# case2_data = {
#     "V1": 1.10,
#     "P2": 450.0,
#     "V2": 1.00,
#     "P3": 1500.0,
#     "Q3": 150.0,
#     "Bs3": 6.0
# }

# case3_data = {
#     "V1": 1.08,
#     "P2": 200.0,
#     "V2": 1.01,
#     "P3": 1500.0,
#     "Q3": 150.0,
#     "Bs3": 6.822 # 
# }

# print("\n\nSolving for Case 2...")
# result_c2 = solve_power_flow(case2_data)
# print(f"Bus 3 Voltage: {result_c2.get('V3', 0):.4f} (Target: 1.0)")
# print(f"LTC Tap: {result_c2.get('LTC_tap_r', 0):.4f}")


# print("\n\nSolving for Case 3...")
# result_c3 = solve_power_flow(case3_data)
# print(f"Bus 3 Voltage: {result_c3.get('V3', 0):.4f} (Target: 1.0)")
# print(f"LTC Tap: {result_c3.get('LTC_tap_r', 0):.4f}")