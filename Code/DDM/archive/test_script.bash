#!/usr/bin/env sh
#PBS -A mnh5174_c_g_sc_default
#PBS -l nodes=${NODES}:ppn=${PPN}
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -M nth7@psu.edu
#PBS -m abe
#PBS -W group_list=mnh5174_collab
#PBS -N test_script

export G=/gpfs/group/mnh5174/default

module use $G/sw/modules

#env
cd $PBS_O_WORKDIR

echo ${NODES}
echo ${PPN}
