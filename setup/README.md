# Setup

Scripts for setting up the environment and compiling the analysis code.

- setupEnv.sh – Sets up the required software environment (ROOT, PYTHIA paths, and environment variables) for running the analysis.
- compile.sh – Compiles the C/C++ analysis code and ROOT macros used in the project.

- main07.cmnd – Base PYTHIA 8 configuration file used as the starting point for proton–proton collision simulations.
- ssbar_monash.cmnd – PYTHIA 8 configuration using the default Monash tune for strange hadron production.
- ssbar_junctions.cmnd – PYTHIA 8 configuration enabling junction-enhanced colour reconnection for strange baryon production.
- ssbar_ropes.cmnd – PYTHIA 8 configuration using the rope hadronisation model to study collective effects.
- ssbar_ee.cmnd – PYTHIA configuration for reference or control simulations (e.g. e⁺e⁻ or simplified baseline setup).

- generate.sh – Runs PYTHIA simulations using the selected configuration files to generate event data.
- analyse.sh – Executes the ROOT-based analysis macros on generated simulation output.
- runCondorJob.sh – Submits large-scale simulation or analysis jobs to an HTCondor batch system.

- ssbar_generate.cpp – Generates strange hadron events with PYTHIA 8 and writes the output used for subsequent analysis.
- ssbar_CorrelationsAnalysis.cpp – Analyses strange hadron correlations and angular distributions using ROOT.
