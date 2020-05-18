#!/usr/bin/env sh
#PBS -A open
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -M nth7@psu.edu
#PBS -m abe
#PBS -N assess_fit


export G=/gpfs/group/mnh5174/default

module use $G/sw/modules

#env
cd $PBS_O_WORKDIR

command -v deactivate >/dev/null 2>&1 && deactivate #exit existing virtual environment if active
module unload python #make sure no system python modules are loaded
module use /gpfs/group/mnh5174/default/sw/modules
#activate python 2.7.9 environment
source /gpfs/group/mnh5174/default/lab_resources/lab_python/bin/activate



######
## sample flanker models
#python /gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/run_flanker_models.py

######
## sample recent probes models
#python /gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/run_recent_probes_models.py

#######
## gelman-rubin checks, chain export, DIC, export 
python /gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/assess_fit.py




