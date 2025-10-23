# Voltage Stability Simulation (Section 6.3 & 8.5 — Van Cutsem & Vournas)

This repository reproduces the dynamic simulations of voltage instability phenomena described in:

> T. Van Cutsem and C. Vournas, *Voltage Stability of Electric Power Systems*, Springer, 1998.  
> § 6.3 “Example System” and § 8.5 “Voltage Instability Phenomena.”

The implementation couples **Julia** (for Differential-Algebraic Equation modeling) and **Python** (for scenario management, plotting, and potential RL integration).

---

## :file_folder: Repository Structure
```
.
├── VoltageStabilityExample.jl    # Julia DAE simulator (core model)
├── case1_simulation.py           # Case 1: base system (Table 6.8, §6.3)
├── case2_simulation.py           # Case 2: with motor load or LTC variation
├── case3_simulation.py           # Case 3: combined contingencies
├── power_flow.jl                 # Julia-based load-flow initialization
├── power_flow.py                 # Python wrapper for power-flow calculation
├── README.md
└── /figures /results             # optional folders for outputs
```

---

## :gear: Environment Setup

### 1️. Install Julia Packages
Open Julia REPL and run:
```julia
using Pkg
Pkg.add([
    "DifferentialEquations",
    "Sundials",
    "NLsolve",
    "Parameters",
    "StaticArrays",
    "Random",
    "SciMLBase"
])
```

### 2️. Python Environment
```bash
pip install juliacall numpy matplotlib
```

Confirm that `julia` is available in your system path.  
`juliacall` will automatically handle inter-language communication.

---

## 3. Run Simulation

Example (Case 1 — single-line tripping):
```bash
python case1_simulation.py
```

This will:
1. Load and compile the Julia module `VoltageStabilityExample.jl`.
2. Build the equilibrium state from § 6.3.
3. Simulate 10 s pre-fault and 50 s post-fault with one-circuit tripping (`X14 × 2`).
4. Plot the bus-2 voltage trajectory showing slow voltage collapse (similar to Fig. 8.26).

---

## Extending Cases

- **Case 2:** Add induction-motor dynamics (`motor.TM > 0`) to represent motor stalling.
- **Case 3:** Enable LTC action (`ltc.mode_sequential=true`) to illustrate long-term instability.
- Each case script reuses the same `VoltageStabilityExample` module with modified `CaseParams`.

---

## Outputs

- Time series of bus voltages (`V₂,V₃,V₄`), generator states (`δ,E′_q,V_fd`), and LTC tap ratio `r`.
- Optional CSV/Matplotlib export for comparison with textbook waveforms (Fig. 8.24 – 8.27).

---

## Power-Flow Initialization

`power_flow.jl` and `power_flow.py` compute the steady-state voltage magnitudes and angles corresponding to the given load levels before transient simulation, ensuring consistency with Table 6.8 data.

---

## Research Use

The simulator provides a clean API suitable for:
- Reinforcement-learning agents controlling AVR/LTC actions.
- Co-simulation with Python for parameter sweeps.
- Sensitivity and bifurcation analysis of voltage stability margins.

---

## Reference
Van Cutsem & Vournas, *Voltage Stability of Electric Power Systems*, Springer, 1998.  
