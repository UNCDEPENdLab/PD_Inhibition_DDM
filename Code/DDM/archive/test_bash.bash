#!/bin/bash
export G=/gpfs/group/mnh5174/default

module use $G/sw/modules

#env
cd $G/Nate/PD_Inhibition_DDM/Code/DDM

command -v deactivate >/dev/null 2>&1 && deactivate #exit existing virtual environment if active
module unload python #make sure no system python modules are loaded
module use /gpfs/group/mnh5174/default/sw/modules
#activate python 2.7.9 environment
source /gpfs/group/mnh5174/default/lab_resources/lab_python/bin/activate



######
## run HDDM models with inputs from PBS args
#python /gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/run_hddm.py ${DF} ${OUTDIR} ${TASK} -m ${MODELS} -nc ${NCHAINS} -nb ${NBURN} -ns ${NSAMP}

python run_hddm.py '/gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Data/preprocessed/flanker_accCode.csv' /gpfs/group/mnh5174/default/Nate/HDDM_outputs_PD_Inhibition/samp1000/flanker/clean_sample flanker -m v
