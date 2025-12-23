# Particle Ratios

This directory contains ROOT macros used to analyse strange hadron particle ratios in proton–proton collisions. The ratios provide a sensitive probe of hadronisation mechanisms and are used to compare different PYTHIA 8 models with experimental ALICE measurements.

## Files

- `baryon_to_meson_ratios.C` – Computes and plots baryon-to-meson ratios (e.g. Λ/K⁰_S, Ξ/π, Ω/π) from PYTHIA simulations to compare hadronisation models.

- `baryon_to_meson_vs_ALICE.C` – Compares PYTHIA 8 baryon-to-meson ratios directly to published ALICE data, enabling validation of model predictions against experimental measurements.

## Notes

- Ratios are used to reduce sensitivity to overall event normalisation and systematic effects.
- Comparisons focus on the default Monash tune and the junction-enhanced colour reconnection model.
- Experimental reference data are taken from ALICE measurements at √s = 7 TeV.
- These macros form the basis of the particle-ratio results discussed in the Bachelor thesis.

