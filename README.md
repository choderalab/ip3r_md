# ip3r_md

This project is focused on preparing and simulating a minimal models of the human type 2 inositol trisphosphate receptor (hIP3R2) containing the juxtamembrane domain (JD), a portion of armadillo repeat 3 (ARM3), and a portion of the S6 helix with or without various bound nucleotides. 

# IP3R

IP3R is an endoplasmic reticulum (ER) resident calcium (Ca2+) channel which plays a central role in intracellular Ca2+ signaling. IP3R is activated by binding of inositol trisphosphate (IP3) and low Ca2+ concentrations, and is inhibited by high Ca2+ concentrations. IP3R activity is potentiated modulated by adenine-containing molecules, but the precise structural mechanism of potentiation has not been well described. The purpose of this project is to use molecular dynamics (MD) simulations to predict impact of adenine-nucleotide binding on the structure of the hIP3R2 domain that binds the adenine nucleotides. 

The general simulation and analysis scheme for this project is as follows. A reduced monomeric portion of hIP3R2 in the resting state containing the JD, a portion of ARM3, and a portion of the S6 helix having either an adenine nucleotide bound or in a ligand-free state was placed in a water box with neutralizing ions using CHARMMGUI (https://charmm-gui.org/). The system was simulated using OpenMM version 7.7.0 (https://openmm.org/). using the CHARMM36m forcefield (https://academiccharmm.org/showcase/natmeth_2016_14_71). 

The impact of ATP binding on JD stability was measured by finding the absolute difference in RMSF of the ATP bound system from the ligand-free system. Based on the RMSF analysis, the euclidean distance between the centers of mass of the backbone atoms of helices α103 and α108 was chosen as a collective variable (CV) for tracking the effect of various nucleotides on JD dynamics. 


# Simulation system setup


