#!/usr/bin/env sh
#PBS -A mnh5174_c_g_sc_default
#PBS -l nodes=2:ppn=20
#PBS -l walltime=168:00:00
#PBS -j oe
#PBS -M nth7@psu.edu
#PBS -m abe
#PBS -W group_list=mnh5174_collab
#PBS -N flanker_small


export G=/gpfs/group/mnh5174/default

module use $G/sw/modules

#env
cd $PBS_O_WORKDIR

command -v deactivate >/dev/null 2>&1 && deactivate #exit existing virtual environment if active
module unload python #make sure no system python modules are loaded
module use /gpfs/group/mnh5174/default/sw/modules
#activate python 2.7.9 environment
source /gpfs/group/mnh5174/default/lab_resources/lab_python/bin/activate

python /gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/run_flanker_models.py
