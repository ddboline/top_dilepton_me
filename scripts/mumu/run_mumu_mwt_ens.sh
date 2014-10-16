#!/bin/bash
# . /etc/bashrc

# . setup_command.sh 

TIME="72:00:00"
MEM="768mb"
# TYPES="os:false ss:true"
# COMMAND=" cafe"
COMMAND="echo runme -cabsrv1 -cput=${TIME} -mem=${MEM} -outdir=/rooms/snapper/projects/ddboline/results_mwt/mwt_weights/mumu -tar=/rooms/snapper/projects/ddboline/me_p18.10.00.tar.gz"
COMMAND=" ./top_dilepton_me/bin/top_dilepton_me_x"

RCP="-rcp=top_dilepton_me/rcp/mumu/mwt_wmt_mumu.rcp"
NPERENS=40

FINAL_STATE="mumu"
# SUFFIX="sigonly"
SUFFIX="sigbkg"

## Ensemble test without systematics

${COMMAND} ${RCP} -cmd=mwt_ens -in_sample=samples_mwt_p20/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ens -out_root=tt_${FINAL_STATE}_ens_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SUFFIX}
