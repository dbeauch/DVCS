# Cross-Section Based Data Generation

This folder contains the generated pseudodata using the BKM02 and BKM10 DVCS cross section formulations at JLab Hall-A and Hall-B DVCS kinematic sets and phi values which are read from the data file `/DataFiles/dvcs_Jlabdata.root`.

Pseudodata generation code is `genpseudoBKM.C`. The resulting cross section values are saved in `.root` and `.csv` files where the name contains the following information:

```ruby
pseudo_inputCFFsmodel_DVCSformulation_KinematicsDetails_twist-approx_xserror
```

The pseudodata generation code takes 4 input arguments. It can be run in ROOT as:

`root -l 'genpseudoBKM.C(1, "BKM10", "Jlab_all", "KM15")'`

###### Input arguments:

* First:

   - [0]() - If CFFs are taken from ANN global fit results which are read from the file `./ANN_GlobalFit_CFFs/BKMXX_ModelFromData.txt`.
         BKMXX is either BKM02 or BKM10 determined by the second argument.
   - [1]() - If taking the CFFs from a model.
* Second:

   Sets the DVCS cross section formulation to be used. At this point we have available the following formulations:   
   - "BKM02"
   - "BKM10"
* Third:

   Experimetal kinematic sets and phi points to be used. The following options can be used:
   - "Jlab_all": All available Hall-A and Hall-B data kinematics. Sets 1 to 195 of the datafile.
   - "HallA_E00-110": Sets 1 to 20 of the data file (20 sets)	--> [https://arxiv.org/abs/1504.05453](https://arxiv.org/abs/1504.05453)
   - "HallA_E12-06-114": Sets 21 to 65 of the data file (45 sets) --> [https://arxiv.org/pdf/2201.03714.pdf](https://arxiv.org/pdf/2201.03714.pdf)
   - "HallB": Sets 66 to 175 of the data file (110 sets) --> [https://arxiv.org/pdf/1504.02009.pdf](https://arxiv.org/pdf/1504.02009.pdf)
   - "HallA_E07–007": Sets 176 to 195 of the data file (20 sets) --> [https://arxiv.org/pdf/1703.09442.pdf](https://arxiv.org/pdf/1703.09442.pdf)

* Fourth:

   If using CFFs from a model, it can be defined in this argument. There are two models implemented:
   - "model1" This is an out-of-the-hat model.
   - "KM15" Modified version of the KM15 model. It is coded up on the file `GPD_Models/TGPDModels.h`.


Resulting cross section distributions are saved on the `./Graph` directory.
