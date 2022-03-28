# Cross-Section Based Data Generation

This folder contains the generated pseudodata using the BKM02 and BKM10 DVCS cross section formulations at JLab Hall-A and Hall-B DVCS kinematic sets and phi values which are read from the data file */DataFiles/dvcs_Jlabdata.root*.

Pseudodata generation code is *genpseudoBKM.C*. The resulting cross section values are saved in .root and .csv files where the name contains the following information:

*pseudo_inputCFFsmodel_DVCSformulation_KinematicsDetails_twist-approx_xserror*

The pseudodata generation code can be run in ROOT as which takes 4 input arguments:

`root -l 'genpseudoBKM.C(1, "BKM10", "Jlab_all", "KM15")'`

##### Input arguments:
