# Transverse Momentum (pT) Spectra

This directory contains ROOT macros used to analyse and visualise the transverse momentum (pT) spectra of strange hadrons produced in proton–proton collisions. The spectra are used to compare different PYTHIA 8 hadronisation models and to study strange hadron production mechanisms.

## Files

- `pt_spectra_k0s.C` – Produces pT spectra for K⁰_S mesons from PYTHIA simulations.

- `pt_spectra_lambda.C` – Produces pT spectra for Λ baryons, enabling comparison between Monash and junction-enhanced models.

- `pt_spectra_xi.C` – Produces pT spectra for Ξ baryons to study strange baryon production mechanisms.

- `pt_spectra_omega.C` – Produces pT spectra for Ω baryons, focusing on multi-strange baryon production.

- `pt_spectra_focused.C` – Generates focused pT spectrum comparisons in selected kinematic ranges or multiplicity classes for detailed model comparison.

## Notes

- All spectra are obtained from PYTHIA 8 simulations at √s = 7 TeV.
- Comparisons between the default Monash tune and the junction-enhanced colour reconnection model are a central focus.
- These macros form the basis for the pT-dependent results presented in the Bachelor thesis.
