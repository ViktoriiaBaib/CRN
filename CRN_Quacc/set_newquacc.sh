#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda create --name newquacc python=3.10 -y

conda activate newquacc

#ASE
pip install --no-cache-dir https://gitlab.com/ase/ase/-/archive/master/ase-master.zip

# Install development version of quacc
pip install git+https://github.com/quantum-accelerators/quacc.git

#QChem
conda install -c conda-forge openbabel -y

pip install quacc[sella]

pip install quacc[tblite]

pip install quacc[jobflow]

#custodian from Samuel Blau
cd
#git clone https://github.com/samblau/custodian.git
cd custodian
git checkout qchem
git pull
python setup.py install

#atomate from Samuel Blau
cd 
#git clone https://github.com/samblau/atomate.git
cd atomate
git checkout qchem
python setup.py install
cd