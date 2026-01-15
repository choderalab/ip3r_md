# Production

To continue simulations following energy minimization and equilibration:

1. Download CHARMM topology and parameter files from: https://charmm-gui.org/?doc=toppar and ensure that simulation scripts (`*.py` in this directory) can properly point to a directory containing these files.
2. For the ligand-free system, run `/4production/ligand-free/sim_zn_500ns.py`. For systems with ATP, ADP, or AMP bound, run `/4production/liganded_ATP/sim_ATP_500ns.py`. For systems with adenosine or guanosine bound, run `/4production/liganded_adenosine_guanosine/sim_nuc_500ns.py`. For the cAMP-bound system, run `/4production/ligand_cAMP/sim_cAMP_500ns.py`. Ensure that the proper dependencies are installed and that the script(s) point to the proper directories/files.
   
