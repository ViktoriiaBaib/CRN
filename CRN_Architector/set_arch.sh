#!/bin/bash

# List of packages and their versions
declare -a packages=(
"numpy==1.24.4"
"pandas==2.0.3"
"scipy==1.11.1"
"tqdm==4.65.0"
"openbabel==3.1.1"
"numba==0.57.1"
"pynauty==2.8.6"
"jupyter==1.0.0"
"ase==3.22.1"
"xtb-python==22.1"
"py3Dmol==2.0.1"
"mendeleev==0.14.0"
"pympipool==0.7.0"
)

source ~/miniconda3/etc/profile.d/conda.sh

conda create --name arch python=3.10 -y

conda activate arch

# Loop through packages and install them
for pkg in "${packages[@]}"
do
    conda install -c conda-forge $pkg -y
done

pip install Sella

conda install -c conda-forge architector -y

cd 

# Switch to developer version of Architector to optimize with Sella
git clone https://github.com/lanl/Architector.git

cd Architector
git checkout Secondary_Solvation_Shell
conda env update --file environment.yml --name arch
pip install -e .

#Switch to dev version of pympipoll
cd
git clone https://github.com/pyiron/pympipool.git
git checkout parallel_executable
pip install -e .