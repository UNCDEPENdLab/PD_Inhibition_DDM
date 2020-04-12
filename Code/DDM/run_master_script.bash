#!/usr/bin/env sh

R CMD BATCH hddm_pipelines_qsub_master.R "pbs_outputs/$(date +"%m-%d-%y-%T.Rout")"

