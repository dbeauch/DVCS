# Cross-Section Based Data Generation

This folder contains the generated pseudodata using the BKM02 and BKM10 DVCS cross section formulations at JLab Hall-A and Hall-B DVCS kinematic sets and phi values which are read from the data file */DataFiles/dvcs_Jlabdata.root*.

Pseudodata generation code is *genpseudoBKM.C*. The resulting cross section values are saved in .root and .csv files where the name contains the following information:

*pseudo_inputCFFsmodel_DVCSformulation_KinematicsDetails_twist-approx_xserror*

The pseudodata generation code takes 4 input arguments. It can be run in ROOT as:

`root -l 'genpseudoBKM.C(1, "BKM10", "Jlab_all", "KM15")'`

###### Input arguments:

* First:
```diff
+ 0
- 1
```
```ruby 0```CFFs are taken from ANN global fit results which are read from the file ./ANN_GlobalFit_CFFs/BKMXX_ModelFromData.txt.
   - 0   CFFs are taken from ANN global fit results which are read from the file ./ANN_GlobalFit_CFFs/BKMXX_ModelFromData.txt.
         BKMXX is either BKM02 or BKM10 determined by the second argument.
   - 1  If taking the CFFs from a model.
* Second:

   Sets the DVCS cross section formulation to be used. At this point we have available the following formulations:   
   - "BKM02"
   - "BKM10"
* Third:

   Experimetal kinematic sets and phi points to be used. The following options can be used:
   - "E00-110 Hall-A" sets 1 to 20 of the data file (20 sets)	--> [https://arxiv.org/abs/1504.05453](https://arxiv.org/abs/1504.05453)
   - "E12-06-114 Hall-A" sets 21 to 65 of the data file (45 sets) --> [https://arxiv.org/pdf/2201.03714.pdf](https://arxiv.org/pdf/2201.03714.pdf)
   - "Hall-B" sets 66 to 175 of the data file (110 sets) --> [https://arxiv.org/pdf/1504.02009.pdf](https://arxiv.org/pdf/1504.02009.pdf)
   - "E07–007 HAll-A" sets 176 to 195 of the data file (20 sets) --> [https://arxiv.org/pdf/1703.09442.pdf](https://arxiv.org/pdf/1703.09442.pdf)
