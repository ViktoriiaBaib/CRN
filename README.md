# CRN

Text-mining analysis is provided as a python notebook in CRN_Text_mining. Synthesis recipe dataset is also released.

1. CRN generation of intermediate graphs

Install SCINE Molassembler for python.

Go to CRN_Molassembler_frag_recomb and run a python script.

2. CRN generation of 3d conformers with Architector (on a cluster)

Install miniconda.
Check that you have file ~/miniconda3/etc/profile.d/conda.sh

Set up a conda environment arch (on a cluster):
```chmod +x set_arch.sh```
```./set_arch.sh```

Go to CRN_Architector and run a python script.

3. Calculations

If you want to run calculations with QuAcc, you can benefit from our installation and example given in CRN_Quacc.

Install the environment with:
```chmod +x set_newquacc.sh```
```./set_newquacc.sh```

Use build_test.py as an example.