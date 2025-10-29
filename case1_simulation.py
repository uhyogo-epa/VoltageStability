# -*- coding: utf-8 -*-
"""
# case1_simulation.py
"""
from juliacall import Main as jl
import pandas as pd
import matplotlib.pyplot as plt

# Load julia module
jl.include("VoltageStabilityExample.jl")
jl.seval("using .VoltageStabilityExample")
VSE = jl.VoltageStabilityExample

# Parameter settings (see Table 6.8)
params = VSE.CaseParams()
params.gen.Pm   =  300 / params.gen.Sbase_sys
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
T_all = 200.0      

# Pre-contingency
VSE.reset(env)
obs_pre = VSE.simulate(env, T=T_pre, dt=dt)
obs_pre = pd.DataFrame(dict(obs_pre)).set_index("t")

# Post-contingency
env.p.line.X14 *= 2 # one circuit tripping
obs_post = VSE.simulate(env, T=T_all-T_pre, dt=dt)
obs_post = pd.DataFrame(dict(obs_post)).set_index("t")
obs      = pd.concat([obs_pre, obs_post])

######################################
# Plot result 
######################################
#selected_keys = ["V2", "V3", "r", "δ", "Eq"]

selected_keys = {
    "V2": "auto", 
    "V3": "auto", 
    "r":  "auto", 
    "δ":  (-3.14, 3.14), 
    "Eq": "auto"
    }


def plot_time_series(obs, selected_keys):

    keys = list(selected_keys.keys())
    n = len(keys)
    
    fig, axes = plt.subplots(
        n, 1, figsize=(8, 1.8 * n), sharex=True, constrained_layout=True
    )

    if n == 1:
        axes = [axes]  # Ensure iterable

    for ax, key in zip(axes, selected_keys):
        if key not in obs.columns:
            print(f"[Warning] Key '{key}' not found in DataFrame.")
            continue
        # plot
        ax.plot(obs.index, obs[key], lw=1.5)
        ax.set_ylabel(key, rotation=0, labelpad=25, fontsize=11)
        ax.grid(True, linestyle="--", alpha=0.4)
        
        # range
        yr = selected_keys[key]
        if isinstance(yr, tuple) and len(yr) == 2:
           ax.set_ylim(*yr)
        elif isinstance(yr, str) and yr.lower() == "auto":
           pass

    axes[-1].set_xlabel("Time [s]", fontsize=11)
    fig.suptitle("Voltage Stability Time Series (case1)", fontsize=13, y=1.02)
    plt.show()

plot_time_series(obs, selected_keys)

