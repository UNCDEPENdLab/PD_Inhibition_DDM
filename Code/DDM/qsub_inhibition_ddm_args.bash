#!/usr/bin/env sh
#PBS -A mnh5174_c_g_sc_default
#PBS -l nodes=${NODES}:ppn=${PPN}
#PBS -l walltime=${WT}:00:00
#PBS -j oe
#PBS -M nth7@psu.edu
#PBS -m abe
#PBS -W group_list=mnh5174_collab
#PBS -N ${JOB_ID_NH}


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
## run HDDM models with inputs from PBS args
python /gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/run_hddm.py ${DF} ${OUTDIR} ${TASK} -m ${MODELS} -nc ${NCHAINS} -nb ${NBURN} -ns ${NSAMP} 
