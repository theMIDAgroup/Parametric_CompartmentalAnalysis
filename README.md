# Parametric_CompartmentalAnalysis
Matlab Graphic User Interface (GUI) for parametric compartmental analysis of kidneys and tumor for the FDG uptake in murine models 

The code is based on the following publication:

* 2017 Scussolini M, Garbarino S, Sambuceti S, Caviglia G and Piana M "A physiology--based parametric imaging method for FDG--PET data", Inverse Problems 33 125010 


A GUI is provided for the following compartmental systems:

    a standard 2D compartmental model for reversible FDG uptake. See for instance [2002 Gunn, R.N., Gunn, S.R., Turkheimer, F.E., Aston, J.A. and Cunningham, V.J. "Positron emission tomography compartmental models: a basis pursuit strategy for kinetic modeling", Journal of Cerebral Blood Flow & Metabolism, 22(12), pp.1425-1439] for oncological applications.

    a (3+1)D compartmental model for FDG kinetic in the kidneys [2014 Garbarino S, Caviglia G, Sambuceti G, Benvenuto F and Piana M "A novel description of FDG excretion in the renal system: application to metformin-treated models" Physics in Medicine and Biology, 59, 2469-2484]

# Usage:

Code is written in Matlab R2015b and tested with versions up to R2019b.

run compartmental_gui.m file launches a Matlab GUI for tumor or kidneys parametric compartmental modelling. It requires 4D .img.hdr data and .voistat for the IF. 
Test data can be found in the "ALBIRA data" folder. Pictorial representation of the 2 compartmental models implemented (2D and (3+1)D) are available in the "models" folder, and accessible from the GUI itself.
