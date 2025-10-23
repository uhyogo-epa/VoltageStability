#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 12:38:12 2025

@author: hashima
"""

module power_flow

using NLsolve
using LinearAlgebra

export solve_power_flow

# システムの固定パラメータ
struct SystemParameters
    X14::Float64 # 母線1-4間リアクタンス [cite: 788]
    X24::Float64 # 母線2-4間変圧器リアクタンス [cite: 788]
    X34::Float64 # 母線3-4間変圧器リアクタンス [cite: 790]
    n::Float64   # 母線2-4間変圧器タップ比 (固定) [cite: 789]
    V3_target::Float64 # 母線3のLTC電圧目標値 
end

# ケースごとの入力データ
mutable struct CaseData
    V1::Float64 # 母線1 電圧 (pu)
    P2::Float64 # 母線2 有効電力 (pu)
    V2::Float64 # 母線2 電圧 (pu)
    P3::Float64 # 母線3 有効電力 (pu)
    Q3::Float64 # 母線3 無効電力 (pu)
    Bs3::Float64 # 母線3 並列コンデンサ (pu)
end

"""
struct GenParams
    # Parameters from Section 6.7 (500-MVA -> 100-MVA base)
    Xd::Float64     
    Xq::Float64      
    Xd_p::Float64    
    # AVR parameters
    G::Float64       
end
"""


"""
タップ比rに基づいてアドミタンス行列(B行列)を構築する。
G (コンダクタンス) は全て0と仮定。
 (Table 6.2) のπ型モデルに基づく。
"""
function build_B_matrix(r::Float64, params::SystemParameters, case_data::CaseData)
    
    n = params.n
    X14 = params.X14
    X24 = params.X24
    X34 = params.X34
    Bs3 = case_data.Bs3

    # 1. ブランチアドミタンス (y_ij = -Y_ij)
    # y = g + jb. ここでは g=0 なので、サセプタンス b のみ
    b_14 = 1.0 / X14
    b_24 = 1.0 / (n * X24) # [cite: 282]
    b_34 = 1.0 / (r * X34) # [cite: 283]

    # 2. シャントアドミタンス (y_ii_sh)
    b_sh_1 = 0.0
    b_sh_2 = (1.0 - n) / (n * X24) # [cite: 282]
    b_sh_3 = (1.0 - r) / (r * X34) #+ Bs3 # [cite: 283]
    b_sh_4 = (n - 1.0) / (n^2 * X24) + (r - 1.0) / (r^2 * X34) # [cite: 282, 283]

    # 3. B行列 (Ybus の虚部)
    # B_ij = -Im(Y_ij)
    # Y_ij = -y_ij (i != j)
    # Y_ii = sum(y_ij) + y_ii_sh
    
    B = zeros(4, 4)

    # 対角要素 Bii = Im(Yii)
    B[1, 1] = b_14
    B[2, 2] = b_24 + b_sh_2 # (1/(n*X24)) + ((1-n)/(n*X24)) = n/(n*X24) = 1/X24
    B[3, 3] = b_34 + b_sh_3 # (1/(r*X34)) + ((1-r)/(r*X34)) + Bs3 = r/(r*X34) + Bs3 = 1/X34 + Bs3
    B[4, 4] = b_14 + b_24 + b_34 + b_sh_4

    # 非対角要素 Bij = Im(Yij)
    B[1, 4] = -b_14
    B[4, 1] = -b_14

    B[2, 4] = -b_24
    B[4, 2] = -b_24

    B[3, 4] = -b_34
    B[4, 3] = -b_34

    return B
end

"""
電力系統のミスマッチ方程式 F(x) = 0
x = [theta2, theta3, theta4, V3, V4, r]
"""
function power_flow_mismatch!(F::Vector, x::Vector, params::SystemParameters, case_data::CaseData)
    
    θ2, θ3, θ4 = x[1], x[2], x[3]
    V3, V4 = x[4], x[5]
    r = x[6]

    V = [case_data.V1, case_data.V2, V3, V4]
    θ = [0.0, θ2, θ3, θ4] 

    B = build_B_matrix(r, params, case_data)
    
    P_calc = zeros(4)
    Q_calc = zeros(4)

    # 標準的な潮流計算式 (G=0 の場合)
    # P_i = sum_j V_i V_j * B_ij * sin(th_i - th_j)
    # Q_i = -sum_j V_i V_j * B_ij * cos(th_i - th_j)
    for i in 1:4
        for j in 1:4
            θ_ij = θ[i] - θ[j]
            P_calc[i] += V[i] * V[j] * B[i, j] * sin(θ_ij)
            Q_calc[i] += -V[i] * V[j] * B[i, j] * cos(θ_ij)
        end
    end

    # 5. ミスマッチベクトルFを計算 (インデックスを修正)
    F[1] = case_data.P2 + P_calc[2] # dP2
    F[2] = case_data.P3 + P_calc[3] # dP3
    F[3] = 0.0 + P_calc[4]         # dP4
    F[4] = case_data.Q3 + Q_calc[3] # dQ3 (★ 3を参照するよう修正)
    F[5] = 0.0 + Q_calc[4]         # dQ4
    F[6] = params.V3_target - V[3] # dV3
end


function solve_th!(res::Vector,x::Vector,V1::Float64,P1::Float64,Q1::Float64)
    
    Eth,θth = x
    
    Xth = 0.01

    res[1]  = P1 - Eth * V1 * sin(θth) /Xth
    res[2]  = Q1 - ( Eth * V1 * cos(θth) - V1^2)/Xth  

end

function solve_V2o!(res::Vector,x::Vector,V2::Float64,θ2::Float64,P2::Float64,Q2::Float64)
    δ,V2o = x
    
    Xq = 2.1 * 100 / 500
    Xd = 2.1 * 100 / 500
    G  = 50.0
    
    
    Xq_in = 1 / Xq
    Xd_in = 1 / Xd
    
    res[1] =  P2 / V2 - G * (V2o - V2) * sin(δ-θ2) / Xd - V2 *(Xq_in - Xd_in) * sin(2*(δ-θ2)) /2 
    res[2] =  Q2 / V2 - G * (V2o - V2) * cos(δ-θ2) / Xd + V2 *(sin(δ-θ2)^2 / Xq + cos(δ-θ2)^2 / Xd)
    
end

"""
Python (juliacall) から呼び出すメイン関数
"""
function solve_power_flow(case_data_dict)
    
    # 1. 固定パラメータを設定
    params = SystemParameters(
        0.0277, # X14 [cite: 788]
        0.016,  # X24 [cite: 788]
        0.004,  # X34 [cite: 790]
        1.04,   # n [cite: 789]
        1.0     # V3_target 
    )


  

    # 2. 入力辞書からCaseDataを作成 (pu変換)
    #  (Table 6.8) の MW/Mvar を 100MVAベース [cite: 786] の pu に変換
    baseMVA = 100.0
    case_data = CaseData(
        case_data_dict["V1"],
        case_data_dict["P2"] / baseMVA,
        case_data_dict["V2"],
        case_data_dict["P3"] / baseMVA,
        case_data_dict["Q3"] / baseMVA,
        case_data_dict["Bs3"]
    )

    # 3. NLsolve の設定
    # f! は F(x) を F に「インプレースで」変更する関数
    f!(F, x) = power_flow_mismatch!(F, x, params, case_data)
    
    # 初期推定値: [θ2, θ3, θ4, V3, V4, r]
    initial_x = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    # 4. 非線形方程式を解く
    try
        result = nlsolve(f!, initial_x, method = :trust_region)

        if !converged(result)
            return Dict("error" => "Solver did not converge.", "details" => result)
        end

        # 5. 収束した場合、結果を整形
        x_sol = result.zero
        θ2, θ3, θ4 = x_sol[1], x_sol[2], x_sol[3]
        V3, V4 = x_sol[4], x_sol[5]
        r = x_sol[6]
        #println(θ2)
        # 6. スラック母線(Bus 1)とPV母線(Bus 2)の電力を計算
        V_final = [case_data.V1, case_data.V2, V3, V4]
        θ_final = [0.0, θ2, θ3, θ4]
        B_final = build_B_matrix(r, params, case_data)
        
        P1_calc = 0.0
        Q1_calc = 0.0
        Q2_calc = 0.0

        for j in 1:4
            θ_1j = θ_final[1] - θ_final[j]
            P1_calc -= V_final[1] * V_final[j] * B_final[1, j] * sin(θ_1j)
            Q1_calc += V_final[1] * V_final[j] * B_final[1, j] * cos(θ_1j)

            θ_2j = θ_final[2] - θ_final[j]
            Q2_calc += V_final[2] * V_final[j] * B_final[2, j] * cos(θ_2j)
        end


        
        g!(res,x) = solve_th!(res,x,V_final[1],P1_calc,Q1_calc)
        
        init_x = [1.0,0] #Eth,θth
        th_sol = nlsolve(g!,init_x,method =:newton)
        th_sols = th_sol.zero
        
        #println(sols[1])

        h!(res,x) = solve_V2o!(res,x,case_data.V2,θ_final[2],case_data.P2,Q2_calc)
        init_xs = [0,1.0] #δ,V2o
        sol = nlsolve(h!,init_xs,method =:newton)
        sols = sol.zero

        #println(solutions[1])
        
        for i in 1:4
            θ_final[i] -= th_sols[2]
        
        end
        
        return Dict(
            "status" => "Converged",
            "V1" => V_final[1], "θ1_deg" => θ_final[1],
            "V2" => V_final[2], "θ2_deg" => θ_final[2],
            "V3" => V_final[3], "θ3_deg" => θ_final[3],
            "V4" => V_final[4], "θ4_deg" => θ_final[4],
            "LTC_tap_r" => r,
            "P1_MW" => P1_calc * baseMVA,
            "Q1_Mvar" => Q1_calc * baseMVA,
            "Q2_Mvar" => Q2_calc * baseMVA,
            "Eth" => th_sols[1], "θth" => th_sols[2],
            "δ" => sols[1] - th_sols[2], "V2o" =>sols[2]
        )

    catch e
        return Dict("error" => "An exception occurred during solving.", "details" => string(e))
    end
end

end # module VanCutsemPowerFlow