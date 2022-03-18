# DVCS Experimental Data

This folder contains the available JLab Hall-A and Hall-B unpolarized dvcs data.

Data is stored in a .csv file with following format: ###### #Set,index,k,QQ,x_b,t,phi_x,F,sigmaF,varF,F1,F2
In this file each set has 24 angles where for the sets that do not have data for a particular \phi value, F is set to zero.

`gentree_Jlabdata.C` saves the .csv datafile of all Jlab data into a `TTree` in a .root file. Every `TTree` entry corresponds to one kimenatic set.

There are 195 kinematic set for the unpolarized data; 85 sets from Hall-A and 110 sets from Hall-B.

###### Data references:
* E00-110 Hall-A sets 1 - 20 (20 sets)	--> [https://arxiv.org/abs/1504.05453](https://arxiv.org/abs/1504.05453)
* E12-06-114 Hall-A sets 21 - 65 (45 sets) --> [https://arxiv.org/pdf/2201.03714.pdf] (https://arxiv.org/pdf/2201.03714.pdf )
* Hall-B sets 66 - 175 (110 sets) --> [https://arxiv.org/pdf/1504.02009.pdf](https://arxiv.org/pdf/1504.02009.pdf)
* E07–007 HAll-A sets 176 - 195 (20 sets) --> [https://arxiv.org/pdf/1703.09442.pdf](https://arxiv.org/pdf/1703.09442.pdf)	
