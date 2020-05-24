#!/usr/bin/env sh
#PBS -A wff3_a_g_hc_default
#PBS -l nodes=1:ppn=5
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -o test.out
#PBS -W group_list=mnh5174_collab
#PBS -N test_job

export G=/gpfs/group/mnh5174/default

module use $G/sw/modules

#env
cd $PBS_O_WORKDIR

command -v deactivate >/dev/null 2>&1 && deactivate #exit existing virtual environment if active
module unload python #make sure no system python modules are loaded
module use /gpfs/group/mnh5174/default/sw/modules

module load python/3.6.3

#activate python 3 environment containing checkout of developer hddm
source /gpfs/group/mnh5174/default/Nate/PD_Inhibition_DDM/Code/DDM/gng_sim_test/hddm_local/bin/activate

# check that we are in fact using the right version of python
which python

echo "testing script run"
python testing_script.py
