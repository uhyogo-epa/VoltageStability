# -*- coding: utf-8 -*-
"""
# case1_simulation.py
"""
from juliacall import Main as jl
import numpy as np
import matplotlib.pyplot as plt

# Load julia module
jl.include("VoltageStabilityExample.jl")
jl.seval("using .VoltageStabilityExample")
VSE = jl.VoltageStabilityExample

# Parameter settings (see Table 6.8)
params = VSE.CaseParams()
params.gen.Pm   =  200 / params.gen.Sbase_sys
params.load.P0  = 1500 / params.gen.Sbase_sys
params.load.Q0  = 0.5 * params.load.P0
params.load.α   = 1.5
params.load.β   = 2.5
params.load.Bs3 = 6.0
params.motor.TM =   0 / params.gen.Sbase_sys
params.P3 = params.load.P0 + params.motor.TM
params.Q3 = 150 / params.gen.Sbase_sys

# Initalization
env = VSE.SimEnv(params)

#######################################
# Simulation
#######################################
dt = 0.1
T_pre = 10.0
T_all = 60.0      

# Pre-contingency
obs_A = VSE.reset(env)
t_pre, obs_pre = VSE.simulate(env, T=T_pre, dt=dt)

# Post-contingency
env.p.line.X14 *= 2 # one circuit tripping
t_post, obs_post = VSE.simulate(env, T=T_all-T_pre, dt=dt)
t, obs = jl.vcat(t_pre, t_post), jl.hcat(obs_pre, obs_post)

# Plot results
V2 = np.array(obs[0,:])
plt.figure(figsize=(10, 6))
plt.plot(t, V2)
plt.xlabel("Time (s)")
plt.ylabel("Bus 2 Voltage (pu)")
plt.title("Voltage Stability Simulation (Case 1)")
plt.grid(True)
plt.show()
