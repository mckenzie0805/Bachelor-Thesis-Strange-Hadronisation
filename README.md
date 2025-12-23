# Strange Hadron Production in Proton–Proton Collisions

This repository contains the simulation and analysis code used for my Bachelor thesis at Maastricht University:

**“Strange Hadron Production in Proton–Proton Collisions: Evaluating PYTHIA 8 Hadronisation Models Against ALICE Measurements.”**

The project investigates strange meson and baryon production in proton–proton collisions at √s = 7 TeV, comparing the default PYTHIA 8 Monash tune with the junction-enhanced colour reconnection model. Subsequently, I also compare both models to published data from the LHC at CERN.

---

## Scientific Context

Hadronisation is a non-perturbative process in QCD and remains one of the least understood aspects of particle physics.  
This work focuses on strange hadrons (K, Λ, Ξ, Ω) as sensitive probes of hadronisation mechanisms.

The analysis compares:
- Particle yields as a function of charged-particle multiplicity
- Transverse momentum (pT) spectra
- Baryon-to-meson ratios
- Two-particle angular correlations

Simulation results are validated against experimental measurements from the ALICE collaboration.

---

## Tools and Frameworks

- **PYTHIA 8** for Monte Carlo event generation  
- **ROOT (C++ macros)** for data analysis and visualisation  
- **HTCondor** for large-scale event generation on a computing cluster  

---

## Repository Structure

- `*.C`, `*.cpp`  
  ROOT C++ analysis macros for:
  - pT spectra
  - Yield extraction
  - Baryon-to-meson ratios
  - Multiplicity-dependent studies
  - Model-to-ALICE comparisons

- `*.sh`  
  Shell scripts used for compilation, job submission, and analysis automation.

- `*.cmnd`  
  PYTHIA configuration files defining simulation parameters and hadronisation models.

---

## Notes on Code Organisation

This repository reflects an active research workflow:
- Multiple versions of macros are retained to track analysis development.
- Scripts are designed for batch execution on a computing cluster.
- Code prioritises transparency and reproducibility over software abstraction.

The intent is to document and preserve the full analysis pipeline used in the thesis.

---

## Thesis

The full Bachelor thesis PDF is included in this repository for scientific context and interpretation of the results.

---

## Author

McKenzie Pedro  
Bachelor of Science, Maastricht University  

