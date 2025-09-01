# Data-Driven Computational Mechanics (DDCM) for Electric Circuits ⚡

This repository contains MATLAB implementations of a **data-driven computational framework** for the simulation of electric circuits.  
The approach replaces traditional constitutive models with **finite sets of data points** (synthetic or experimental), enforcing them through a feedback operator within a **variational time-integration scheme**.

---

## ✨ Features
- Circuit assembly from **graph representations** with inductors and capacitors.  
- Non-holonomic constraints handling via reduced incidence matrices.  
- **Alternating Direction Method (ADM)** solver for DDCM applied to RLC-type networks.  
- Synthetic capacitor data generation with nonlinear effects.  
- Tools for post-processing:
  - Phase-space plots (charge–voltage pairs vs. data)  
  - Time evolution of charges  
  - Residuum, error cost, and global cost function histories  

---

## 📂 Repository Structure
- `main.m` – Example script that assembles and simulates a circuit using DDCM.  
- `circuit_assembly.m` – Builds circuit matrices (L, C, constraints).  
- `KKT_Time.m` – Constructs KKT matrices for time integration.  
- `phiESOperator.m` – Feedback operator: projects trial states to closest data points.  
- `DDCM_ADM_Electric.m` – Core solver implementing DDCM with ADM.  

---

## 🚀 Getting Started

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/DDCM-ElectricCircuits.git



