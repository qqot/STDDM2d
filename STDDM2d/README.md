# Domain Decomposition Method (DDM) for Large-Scale Electromagnetic Simulations
This repository provides the source codes accompanying the paper:

> **"Parallel overlapping-domain decomposition FDFD for large-scale complex nanostructures modeling"**  
> Z. Wang *et al.*, 2025.  
> (Manuscript currently under review.)

## 🔍 Overview

This project implements a parallel overlapping domain decomposition method (DDM) based on the finite-difference frequency-domain (FDFD) formulation to model the electromagnetic response of large-scale complex nanostructures. The global computational domain is partitioned into multiple overlapping subdomains terminated with perfectly matched layers (PMLs), enabling seamless source transfer between adjacent subdomains.

## ⚙️ Requirements

- MATLAB **R2022a** or later  
- **Parallel Computing Toolbox**  
- (Optional) Multi-core or multi-node cluster for `spmd` parallel runs  

## 🧠 Citation and Usage Policy
If you use, adapt, or build upon this code in your research, please cite the following paper:

Z. Wang, et al.
Parallel overlapping-domain decomposition FDFD for large-scale complex nanostructures modeling, 2025.