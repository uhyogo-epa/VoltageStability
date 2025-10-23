module VoltageStabilityExample
"""
Voltage Stability Example – time-domain simulator (Julia)
Chapter 6 §6.3 Example
from Van Cutsem & Vournas, Voltage Stability of Electric Power Systems.

Design goals
- DAE formulation (semi-explicit):  ẋ = f(x,y,z),  0 = g(x,y,z)
- Discrete LTC logic via callbacks
- Clean API for Python (juliacall/pyjulia) and RL: reset!, observe, step!(action), simulate!

Dependencies
  add DifferentialEquations Sundials NLsolve Parameters StaticArrays Random

This file provides a minimal, *working* scaffold with numerics in per-unit.
Fill/adjust the numerical parameters (network B’s, machine constants, load mix) to your case.
"""

using Printf
using DifferentialEquations
using Sundials
using NLsolve
using Parameters
using StaticArrays
using Random
using SciMLBase

#########################################################################
# Types and Parameters
#########################################################################
@with_kw mutable struct LineParams
    # Parameters from Section 6.7
    Xth::Float64 = 0.01    
    X14::Float64 = 0.0277
    X24::Float64 = 0.016
    X34::Float64 = 0.004
    n::Float64   = 1.04
    Eth::Float64 = 1.10
end

@with_kw mutable struct GenParams
    # Parameters from Section 6.7 (500-MVA -> 100-MVA base)
    Sbase_gen::Float64 = 500.0 # MVA
    Sbase_sys::Float64 = 100.0 # MVA
    Xd::Float64      = 2.1 * (Sbase_sys / Sbase_gen)
    Xq::Float64      = 2.1 * (Sbase_sys / Sbase_gen)
    Xd_p::Float64    = 0.4 * (Sbase_sys / Sbase_gen)
    Tdo_p::Float64   = 8.0       # s
    ω0::Float64      = 2π * 50.0 # System synchronous speed
    H::Float64       = 3.5       # Inertia constant [s]
    D::Float64       = 4.0 * (Sbase_sys / Sbase_gen)
    Pm::Float64      = 200.0 / Sbase_sys # Mechanical power
    # AVR parameters
    G::Float64       = 50.0    # AVR gain
    T::Float64       = 0.1     # AVR time constant [s]
    Vref::Float64    = 1.01    # AVR setpoint
    Vfd_max::Float64 = 5.0
    Vfd_min::Float64 = 0.0
    # OXL parameters
    I_lim::Float64    = 2.825 # Eq_lim in pu
    S1::Float64       = 1.0
    S2::Float64       = 2.0
    K1::Float64       = 20.0
    K2::Float64       = 0.1
    Kr::Float64       = 1.0
    Ki::Float64       = 0.1
    OXL_on::Bool      = true
end

@with_kw mutable struct LTCParams
    r::Float64    = 1.0
    rmin::Float64 = 0.8
    rmax::Float64 = 1.1
    Δr::Float64   = 0.01
    V3o::Float64  = 1.0
    deadband::Float64 = 0.01
    Tjo::Float64 = 20.0 # Initial fixed delay
    Tj::Float64 = 10.0  # Subsequent fixed delay
    Tm::Float64 = 0.0   # Mechanical delay (assuming it's included in Tj)
    mode_sequential::Bool = true
    # internal timers (mutated)
    armed::Bool     = false
    t_next::Float64 = Inf
end

@with_kw mutable struct LoadParams
    # Parameters for Case 1 from Table 6.8
    α::Float64   = 1.5
    β::Float64   = 2.5
    P0::Float64  = 1500.0 / 100
    Q0::Float64  = P0 / 2
    V0::Float64  = 1.0
    Bs3::Float64 = 6.0 # pu on 100-MVA base
end

@with_kw mutable struct MotorParams
    # Parameters from Section 6.7, on 100-MVA system base
    Sbase_motor::Float64 = 800.0 # MVA
    Sbase_sys::Float64   = 100.0 # MVA
    Xs::Float64 = 0.1   * (Sbase_sys / Sbase_motor)
    Xr::Float64 = 0.18  * (Sbase_sys / Sbase_motor)
    Xm::Float64 = 3.2   * (Sbase_sys / Sbase_motor)
    Rr::Float64 = 0.018 * (Sbase_sys / Sbase_motor)
    H::Float64  = 0.5   # s
    TM::Float64 = 0.0   # Mechanical Torque
end

#--------------------------------
# Case Parameters
#---------------------------------
@with_kw mutable struct CaseParams
    line::LineParams = LineParams()
    gen::GenParams = GenParams()
    ltc::LTCParams = LTCParams()
    load::LoadParams = LoadParams()
    motor::MotorParams = MotorParams()
    P3::Float64 = 1500 / gen.Sbase_sys
    Q3::Float64 =  150 / gen.Sbase_sys
end

#-----------------------------
# Utility
#-----------------------------
clamp_(x, lo, hi) = x < lo ? lo : (x > hi ? hi : x)

#########################################################################
# Device models ("injection" currents)
#########################################################################

""" Generator injected currents at Bus 2: Eq.(6.24)"""
function gen_currents(Vx2, Vy2, δ, Eqp, gen::GenParams)
    Xdp, Xq = gen.Xd_p, gen.Xq
    invXdp = 1 / gen.Xd_p
    invXq  = 1 / gen.Xq
    sinδ, cosδ, sin2δ = sin(δ), cos(δ), sin(2δ)
    Ix2 =  ((sin2δ/2) * (invXq - invXdp) * (Vx2 - Eqp * cosδ) 
            - (cosδ^2 / Xq + sinδ^2 / Xdp) * (Vy2 - Eqp * sinδ))
    Iy2 =  ((sinδ^2 / Xq + cosδ^2 / Xdp) * (Vx2 - Eqp * cosδ) 
            + (sin2δ/2) * (invXdp - invXq) * (Vy2 - Eqp * sinδ))
    #println("---------------------------") 
    #println("Ix2 = $(Ix2), Iy2 = $(Iy2)")
    #Vd, Vq = Vx2 * sin(δ) - Vy2 * cos(δ), Vx2 * cos(δ) + Vy2 * sin(δ) 
    #id, iq = (Eqp  - Vq) / gen.Xd_p, Vd / gen.Xq
    #Ix, Iy =  id * sin(δ) + iq * cos(δ), -id * cos(δ) + iq * sin(δ)
    #println("Ix2 = $(Ix), Iy2 = $(Iy)")
    return Ix2, Iy2
end

"""Exponential static load currents at Bus3."""
function exp_load_currents(Vx3, Vy3, lp::LoadParams)
    Vsq  = Vx3^2 + Vy3^2
    VsqS = max(Vsq, 1e-8) # zero-division guard
    V    = sqrt(VsqS)
    P = lp.P0 * (V / lp.V0)^(lp.α)
    Q = lp.Q0 * (V / lp.V0)^(lp.β)
    # current "injected" into the network
    ix = -(P * Vx3 + Q * Vy3) / VsqS
    iy =  (Q * Vx3 - P * Vy3) / VsqS
    return ix, iy
end

"""Induction motor currents at Bus3"""
function motor_currents(Vx3, Vy3, s, mp::MotorParams)
    Xs, Xm, Xr, Rr = mp.Xs, mp.Xm, mp.Xr, mp.Rr

    # Re + jXe per Eq. (4.14)
    denom = Rr/s + im*(Xm + Xr)
    ReXe  = im * Xm * (Rr/s + im*Xr) / denom
    Re, Xe = real(ReXe), imag(ReXe)

    den2 = Re^2 + (Xs + Xe)^2
    IxM  = (-Vx3*Re + Vy3*(Xs + Xe)) / den2
    IyM  = (-Vy3*Re - Vx3*(Xs + Xe)) / den2
    return IxM, IyM
end

"""Shunt capacitor current at Bus3"""
function shunt_comp_currents(Vx3, Vy3, lp::LoadParams)
    # opposite direction to the text (p192) 
    ixC =  lp.Bs3 * Vy3
    iyC = -lp.Bs3 * Vx3
    return ixC, iyC
end


#########################################################################
# Differential algebraic equations
# States order: δ=1, ω=2, Eqp=3, Vfd=4, Xt=5, Xoxl=6, s=7
# y = [Vx1,Vy1, Vx2,Vy2, Vx3,Vy3, Vx4,Vy4]  (real/imag bus voltages)
#########################################################################
const IDX = (δ=1, ω=2, Eqp=3, Vfd=4, Xt=5, Xoxl=6, s=7)
const IDY = (Vx1=1, Vy1=2, Vx2=3, Vy2=4, Vx3=5, Vy3=6, Vx4=7, Vy4=8)

"""Compute susceptances (B_ij and B_sij) according to Table 6.2 (p.187)"""
function calc_network_params(lp::LineParams, r::Float64)
    @unpack X14, X24, X34, n = lp
    # Series susceptances (Table 6.2)
    B14 = -1 / X14
    B24 = -1 / (n * X24)
    B34 = -1 / (r * X34)
    # Shunt susceptances referred to bus 2,3 and bus 4
    Bs24 = (1 - n) / (n * X24)       # seen from bus 2
    Bs42 = (n - 1) / (n^2 * X24)     # seen from bus 4
    Bs34 = (1 - r) / (r * X34)       # seen from bus 3
    Bs43 = (r - 1) / (r^2 * X34)     # seen from bus 4
    return B14, B24, B34, Bs24, Bs42, Bs34, Bs43
end

"""Residuals of network algebraic equations g(y) = 0: Eqs. (6.16a)–(6.16h)"""
function network_residual!(res, y, x, cp::CaseParams)

    # Network parameters
    r = cp.ltc.r
    B14, B24, B34, Bs24, Bs42, Bs34, Bs43 = calc_network_params(cp.line, r)
    Eth, Xth = cp.line.Eth, cp.line.Xth
    θth = 0.0
    
    # Differential states
    δ   = x[IDX.δ]
    Eqp = x[IDX.Eqp]
    s   = x[IDX.s]

    # Bus voltages (real/imag)
    Vx1, Vy1, Vx2, Vy2, Vx3, Vy3, Vx4, Vy4 = y

    # Device current injections
    Ix2, Iy2 = gen_currents(Vx2, Vy2, δ, Eqp, cp.gen)
    ixE, iyE = exp_load_currents(Vx3, Vy3, cp.load)
    ixM, iyM = (cp.motor.TM > 0.0) ? motor_currents(Vx3, Vy3, s, cp.motor) : 0.0, 0.0
    ixC, iyC = shunt_comp_currents(Vx3, Vy3, cp.load)
    Ix3, Iy3 = ixE + ixM + ixC, iyE + iyM + iyC

    # Network equations Eq.(6.16)
    invXth = 1 / Xth
    res[1] =  (Eth*invXth)*sin(θth) - invXth*Vy1 + B14*(Vy1 - Vy4)
    res[2] = -(Eth*invXth)*cos(θth) + invXth*Vx1 - B14*(Vx1 - Vx4)
    res[3] = Ix2 + Bs24*Vy2 + B24*(Vy2 - Vy4)
    res[4] = Iy2 - Bs24*Vx2 - B24*(Vx2 - Vx4)
    res[5] = Ix3 + Bs34*Vy3 + B34*(Vy3 - Vy4)
    res[6] = Iy3 - Bs34*Vx3 - B34*(Vx3 - Vx4)
    res[7] = (Bs42 + Bs43)*Vy4 + B14*(Vy4 - Vy1) + B24*(Vy4 - Vy2) + B34*(Vy4 - Vy3)
    res[8] = -(Bs42 + Bs43)*Vx4 - B14*(Vx4 - Vx1) - B24*(Vx4 - Vx2) - B34*(Vx4 - Vx3)

    return nothing
end

""" Dynamic equation dxdt=f(x,y) """
function f!(dx, x, y, cp::CaseParams, t)
    # --- Unpack parameters and varibles 
    @unpack gen, motor = cp
    δ, ω, Eqp, Vfd, Xt, Xoxl, s = x
    Vx2, Vy2 = y[IDY.Vx2], y[IDY.Vy2]
    Vx3, Vy3 = y[IDY.Vx3], y[IDY.Vy3]

    # --- Generator dynamics: δ,ω (Eq. 6.19, 6.20)
    ω0 = gen.ω0
    Ix2, Iy2 = gen_currents(Vx2, Vy2, δ, Eqp, gen)
    P2 = Vx2 * Ix2 + Vy2 * Iy2
    dx[IDX.δ] = ω
    dx[IDX.ω] = (ω0 / (2 * gen.H)) * (gen.Pm - P2 - (gen.D / ω0) * ω)

    # --- Field voltage dynamics: Eq' (Eq. 6.21) 
    id = Ix2 * sin(δ) - Iy2 * cos(δ)     # Eq. (6.23)
    Eq = Eqp + (gen.Xd - gen.Xd_p) * id  # Eq. (3.12)
    dx[IDX.Eqp] = (Vfd - Eq) / gen.Tdo_p

    # --- AVR dynamics: Vfd (Eq. 6.25)
    V2  = hypot(Vx2, Vy2)
    cmd = gen.G * (gen.Vref - V2 - Xoxl)
    dx[IDX.Vfd] = (cmd - Vfd) / gen.T
    # Apply AVR limits
    if (Vfd >= gen.Vfd_max && dx[IDX.Vfd] > 0) || (Vfd <= gen.Vfd_min && dx[IDX.Vfd] < 0)
        dx[IDX.Vfd] = 0.0
    end

    # --- OXL dynamics: Xt, Xoxl (Fig 3.12, Eq. 6.26-6.28)
    if gen.OXL_on
        # Integrator for timer Xt
        X2 = (Eq >= gen.I_lim) ? gen.S1 * (Eq - gen.I_lim) : gen.S2 * (Eq - gen.I_lim)
        dx[IDX.Xt] = X2
        if (Xt >= gen.K2 && X2 >= 0) || (Xt <= -gen.K1 && X2 < 0)
            dx[IDX.Xt] = 0.0
        end

        # Integrator for output Xoxl
        X3 = (Xt > 0) ? (Eq - gen.I_lim) : -gen.Kr
        dx[IDX.Xoxl] = gen.Ki * X3
        if (Xoxl <= 0 && X3 < 0)
             dx[IDX.Xoxl] = 0.0
        end
    else
        dx[IDX.Xt] = 0.0
        dx[IDX.Xoxl] = 0.0
    end

    # --- Motor dynamics ---
    # Slip dynamics based on torque balance (Eq. 4.22)
    ds = 0.0
    if motor.TM > 0.0  # Activate only if motor is configured
        ixM, iyM = motor_currents(Vx3, Vy3, s, cp.motor)
        Te = -Vx3*ixM -Vy3*iyM
        ds = (motor.TM - Te) / (2 * motor.H)
        println("TM = $(motor.TM), Te = $(Te)")
        
    end
    dx[IDX.s] = ds

    return nothing
end

""" DAE wrapper for Sundials.IDA """
struct DAEWrapper{P}
    p::P
end
# Make DAEWrapper callable
function (F::DAEWrapper)(res, du, u, p, t)
    # u = [x; y]
    # x: δ, ω, Eqp, Vfd, Xt, Xoxl, s  (7 states)
    # y: Vx1...Vy4                   (8 vars)
    x = @view u[1:7]
    y = @view u[8:15]

    dx = similar(x)
    f!(dx, x, y, F.p, t)

    # differential part residuals: du - f = 0
    @inbounds for i in 1:7
        res[i] = du[i] - dx[i]
    end

    # algebraic: g(y)=0
    network_residual!(@view(res[8:15]), y, x, F.p)
    return nothing
end

""" Mass matrix pattern via differential_vars """
differential_vars_mask() = SVector{15}((true, true, true, true, 
                                        true, true, true,
                                        false, false, false, false, 
                                        false, false, false, false))

#########################################################################
# LTC discrete logic 
#########################################################################
function ltc_update!(cp::CaseParams, V3, t)
    ltc = cp.ltc
    err = V3 - ltc.V3o
    
    # Determine deadband (wider for first move if hysteresis is on)
    db = ltc.deadband
    if ltc.mode_sequential && !ltc.armed
        # For sequential mode, first move might have a different delay logic
        # Here we simplify and just use one set of delays
    end

    outside = (err > db && ltc.r < ltc.rmax) || (err < -db && ltc.r > ltc.rmin)

    if !ltc.armed && outside
        # Arm the LTC for the first move
        ltc.armed = true
        ltc.t_next = t + ltc.Tjo # Use initial delay
    elseif ltc.armed
        if t >= ltc.t_next
            # Perform tap change
            s = (err > 0) ? -1.0 : 1.0 # Lower r to increase V2 (tap on HV side)
            ltc.r = clamp(ltc.r + s * ltc.Δr, ltc.rmin, ltc.rmax)
            
            # Check if still outside deadband to schedule next move
            err_new = V3 - ltc.V3o # Re-evaluate error (approximate)
            still_outside = (err_new > db && ltc.r < ltc.rmax) || (err_new < -db && ltc.r > ltc.rmin)

            if ltc.mode_sequential && still_outside
                ltc.t_next = t + ltc.Tj # Schedule next move with subsequent delay
            else
                ltc.armed = false # Disarm
                ltc.t_next = Inf
            end
        end
    end
    return nothing
end


#########################################################################
# Bulid initial state
#########################################################################
function build_initial_state(cp::CaseParams; verbose::Bool=false)
    io = verbose ? stdout : devnull
    @unpack line, gen, ltc, load, motor = cp
    #----------------------------
    # 1) Set initial estimates 
    # ---------------------------
    u0 = [
        -0.1,        # δ: Machine angle relative to reference
        1.08, -0.1,  # V1 ,θ1 
        1.01, -0.2,  # V2 ,θ2
        1.0,  -0.5,  # V3, θ3
        1.0,  -0.5,  # V4 ,θ4
        1.0,  # r
    ]

    #----------------------------
    # 2) Define the single residual function F(u) for the solver
    #----------------------------
    function f_init!(res, u)

        # Unpack state variables
        δ, V1, θ1, V2, θ2, V3, θ3, V4, θ4, r = u
        
        # Generator currents
        Eq = gen.G * (gen.Vref - V2)        
        Xd, Xq = gen.Xd, gen.Xq
        IP2 = (Eq / Xd) * sin(δ - θ2) + (V2 / 2) * (1/Xq - 1/Xd) * sin(2 * (δ - θ2))
        IQ2 = (Eq / Xd) * cos(δ - θ2) - V2 * (sin(δ - θ2)^2 / Xq + cos(δ - θ2)^2 / Xd)        
        Pm = gen.Pm
        
        # Load currents
        IP3 = - cp.P3 / cp.ltc.V3o
        IQ3 = - cp.Q3 / cp.ltc.V3o

        # Network_equations
        Eth, Xth = line.Eth, line.Xth
        θth = 0.0
        B14, B24, B34, Bs24, Bs42, Bs34, Bs43 = calc_network_params(line, r)        
        res[1] = IP2 * V2 - gen.Pm
        res[2] = -(Eth/Xth) * sin(θ1 - θth) + B14 * V4 * sin(θ1 - θ4)
        res[3] =  (Eth/Xth) * cos(θ1 - θth) - V1 / Xth + B14 * V1 - B14 * V4 * cos(θ1 - θ4)
        res[4] = IP2 + B24 * V4 * sin(θ2 - θ4)
        res[5] = IQ2 + (Bs24 + B24) * V2 - B24 * V4 * cos(θ2 - θ4)
        res[6] = IP3 + B34 * V4 * sin(θ3 - θ4)
        res[7] = IQ3 + (Bs34 + B34) * V3 - B34 * V4 * cos(θ3 - θ4)
        res[8] = B14*V1*sin(θ4 - θ1) + B24*V2*sin(θ4 - θ2) + B34*V3*sin(θ4 - θ3)
        res[9] = (B14 + Bs42 + B24 + Bs43 + B34)*V4 - B14*V1*cos(θ4 - θ1) - B24*V2*cos(θ4 - θ2) - B34*V3*cos(θ4 - θ3)
        res[10] = V3 - cp.ltc.V3o

        return nothing
    end

    #----------------------------    
    # 3) Solve F(u) = 0 using a nonlinear solver
    #----------------------------
    println(io,"##################################################")
    println(io,"Solving for equilibrium point...")
    sol = nlsolve(
        u -> (res = similar(u); f_init!(res, u); res),
        u0;
        method = :newton,
        ftol = 1e-10,
        show_trace = verbose,
    )
    if !converged(sol)
        error("Nonlinear solver for equilibrium point did not converge.")
    end
    
    println(io,"Equilibrium found")
    u_ss = sol.zero
    δ, V1, θ1, V2, θ2, V3, θ3, V4, θ4, r = u_ss
    δw = rem2pi(δ, RoundNearest)
    println(io,"δ = $(round(δw, digits=4)),  r = $(round(r, digits=4)) pu")
    println(io,"V1=$(round(V1,digits=4))∠$(round(θ1,digits=4)),  ",
               "V2=$(round(V2,digits=4))∠$(round(θ2,digits=4)),  ",
               "V3=$(round(V3,digits=4))∠$(round(θ3,digits=4)),  ",
               "V4=$(round(V4,digits=4))∠$(round(θ4,digits=4))")    

    #----------------------------
    # 4) Construct initial state for DAE
    #----------------------------
    println(io,"##################################################")           
    println(io,"Constructing DAE initial conditions...")   
    y0 = [
        V1 * cos(θ1), V1 * sin(θ1), # Vx1, Vy1
        V2 * cos(θ2), V2 * sin(θ2), # Vx2, Vy2
        V3 * cos(θ3), V3 * sin(θ3), # Vx3, Vy3
        V4 * cos(θ4), V4 * sin(θ4)  # Vx4, Vy4
    ]
    
    # Vfd, Eq
    Vfd = gen.G * (gen.Vref - V2)
    Eq = Vfd
    Vx2, Vy2 = V2 * cos(θ2), V2 * sin(θ2)
    Vq = Vx2 * cos(δ) + Vy2 * sin(δ)
    id = (Eq - Vq) / gen.Xd
    Eqp = Eq - (gen.Xd - gen.Xd_p) * id
    println(io,"Additional states:")
    println(io,"  Vfd = $(round(Vfd, digits=4)),  Eq' = $(round(Eqp, digits=4)) pu")    
    
    # s (motor slip)
    if motor.TM > 0.0  # Activate only if motor is configured
        @unpack Xs, Xr, Xm, Rr, TM = motor
        X1 = Xm * Xs / (Xm + Xs)
        f!(F, s) = (F[1] = (Rr/s[1]*Xm^2*V3^2)/(((Rr/s[1])^2+(X1+Xr)^2)*(Xs+Xm)^2) - TM)
        s = nlsolve(f!, [0.02]).zero[1]
        println(io,"  s (motor slip) = $(round(s, digits=6))")
    else
        s = 0
        println(io,"  s (motor slip) = $(round(s, digits=6)) (Motor inactive)")
    end
    
    x0 = [
        δw,      # δ
        0.0,     # ω
        Eqp,     # Eqp
        Vfd,     # Vfd
        -gen.K1, # Xt
        0.0,     # Xoxl
        s        # s
    ]
    cp.ltc.r = r
    u0 = vcat(x0, y0)
    du0 = zeros(length(u0))   

    # Check load power    
    Vx3, Vy3 = V3 * cos(θ3), V3 * sin(θ3)
    ixE, iyE = exp_load_currents(Vx3, Vy3, load)
    ixM, iyM = (cp.motor.TM > 0.0) ? motor_currents(Vx3, Vy3, s, motor) : 0.0, 0.0
    ixC, iyC = shunt_comp_currents(Vx3, Vy3, load)

    PE, QE = Vx3 * ixE + Vy3 * iyE, Vy3 * ixE - Vx3 * iyE
    PM, QM = Vx3 * ixM + Vy3 * iyM, Vy3 * ixM - Vx3 * iyM
    PC, QC = Vx3 * ixC + Vy3 * iyC, Vy3 * ixC - Vx3 * iyC
    println(io,"Load powers (injeceted into network):")
    @printf(io,"  PE = %8.4f,  QE = %8.4f\n", PE, QE)
    @printf(io,"  PM = %8.4f,  QM = %8.4f\n", PM, QM)
    @printf(io,"  PC = %8.4f,  QC = %8.4f\n", PC, QC)

    # Check residual      
    wrapper = DAEWrapper(cp)
    res_check = similar(u0)
    wrapper(res_check, du0, u0, cp, 0.0)
    println(io,"DAE Residual vector:")  
    println(io,"  x[1:7] (differential part): [",
        join([@sprintf("% .3e", v) for v in res_check[1:7]], " "), "]")
    println(io,"  y[1:8] (algebraic part): [",
        join([@sprintf("% .3e", v) for v in res_check[8:15]], " "), "]")
        
    return u0, du0, r
end


#########################################################################
# Simulator and RL-style environment
#########################################################################
mutable struct SimEnv
    p::CaseParams
    u0::Vector{Float64}
    du0::Vector{Float64}
    t::Float64
    prob::SciMLBase.AbstractDAEProblem
    sol::Union{Nothing, SciMLBase.AbstractODESolution}
end

function SimEnv(cp::CaseParams; seed::Int=0)
    u0, du0, r = build_initial_state(cp; verbose=true)
    wrapper = DAEWrapper(cp)
    diffs = differential_vars_mask()
    prob = DAEProblem(wrapper, du0, u0, (0.0, 1.0); differential_vars=diffs)
    return SimEnv(cp, u0, du0, 0.0, prob, nothing)
end


"""Observation vector: [V2, V3, V4, δ, Eq', Vfd, r]."""
function observe(env::SimEnv)
    u = env.u0
    x = @view u[1:7]
    y = @view u[8:15]
    Vx2, Vy2 = y[IDY.Vx2], y[IDY.Vy2]
    Vx3, Vy3 = y[IDY.Vx3], y[IDY.Vy3]
    Vx4, Vy4 = y[IDY.Vx4], y[IDY.Vy4]
    V2 = hypot(Vx2, Vy2)
    V3 = hypot(Vx3, Vy3)
    V4 = hypot(Vx4, Vy4)
    return [V2, V3, V4, x[IDX.δ], x[IDX.Eqp], x[IDX.Vfd], env.p.ltc.r]
end
const IDO = (V2=1, V3=2, V4=3, δ=4, Eqp=5, Vfd=6, r=7)


"""Reset environment (optionally randomize load level). Returns initial observation."""
function reset!(env::SimEnv; P0_scale::Float64=1.0)
    # Initial condition
    env.u0, env.du0, r = build_initial_state(env.p)
    
    # LTC
    env.p.ltc.armed  = false
    env.p.ltc.t_next = Inf
    env.p.ltc.r = r
    
    env.t = 0.0
    env.prob = remake(env.prob; u0=env.u0, du0=env.du0, tspan=(0.0, 0.0))
    env.sol = nothing
    return observe(env)
end


"""One integration step of Δt seconds with optional control action.
`action` can set: AVR setpoint shift ΔVref, LTC block flag, or emergency load shed.
Return (obs, reward, done, info).
"""
# --------- helper for Python dict / Julia NamedTuple actions ----------
# Both NamedTuple with Symbol keys and Python dict with string keys are accepted.
_has(a, k::Symbol) = (try
        haskey(a, k) || haskey(a, String(k))
    catch
        false
    end)

_get(a, k::Symbol, default) = (try
        haskey(a, k) ? a[k] : (haskey(a, String(k)) ? a[String(k)] : default)
    catch
        default
    end)

voltage_deviation_penalty(obs) = (abs(obs[1]-1.0) + abs(obs[2]-1.0) + 0.2*abs(obs[3]-1.0))


function step!(env::SimEnv; Δt::Float64=0.1, action=NamedTuple())

    p = env.p
    # apply action (accept Symbol or String keys)
    if _has(action, :ΔVref)
        p.gen.Vref += _get(action, :ΔVref, 0.0)
    end
    if _get(action, :ltc_block, false)
        p.ltc.armed = false
        p.ltc.t_next = Inf
    end
    if _has(action, :shed_frac)
        s = clamp_(_get(action, :shed_frac, 0.0), 0.0, 1.0)
        p.load.P0 *= (1.0 - s)
        p.load.Q0 *= (1.0 - s)
    end

    # Evolve LTC timer based on V3
    V3 = observe(env)[IDO.V3]
    ltc_update!(p, V3, env.t)

    # Integrate DAE from t to t+Δt (Sundials.IDA)
    tspan = (env.t, env.t + Δt)
    wrapper = DAEWrapper(p)
    diffs = differential_vars_mask()

    env.prob = DAEProblem(wrapper, env.du0, env.u0, tspan; differential_vars=diffs)
    sol = solve(env.prob, IDA(); reltol=1e-6, abstol=1e-8)

    #--- Update time and state to last time point
    env.sol = sol
    env.t = sol.t[end]
    env.u0 = Array(sol.u[end])
    env.du0 = Array(sol.du[end])

    obs = observe(env)
    reward = -voltage_deviation_penalty(obs)
    done = false
    info = (;)
    return obs, reward, done, info
end


"""Run open-loop simulation for horizon T with fixed step saveat.
Returns times, matrix of observations.
"""
function simulate!(env::SimEnv; T::Float64=60.0, dt::Float64=0.1)
    n = Int(floor(T/dt))
    obs_dim = length(observe(env))
    obs_hist = zeros(obs_dim, n+1)
    t_hist = zeros(n+1)
    obs_hist[:,1] .= observe(env)
    t_hist[1] = env.t
    for k in 1:n
        obs, r, done, info = step!(env; Δt=dt)
        obs_hist[:,k+1] .= obs
        t_hist[k+1] = env.t
    end
    return t_hist, obs_hist
end

# ---------------------------
# Exports
# ---------------------------

# Python-friendly wrappers (no ! in the name)
reset(env::SimEnv; kwargs...) = reset!(env; kwargs...)
step(env::SimEnv; kwargs...) = step!(env; kwargs...)
simulate(env::SimEnv; kwargs...) = simulate!(env; kwargs...)


export CaseParams, SimEnv, reset!, step!, observe, simulate!, reset, step, simulate, build_initial_state

end # module

