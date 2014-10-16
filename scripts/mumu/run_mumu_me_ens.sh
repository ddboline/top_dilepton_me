#!/bin/bash
# . /etc/bashrc

# . setup_command.sh 

TIME="72:00:00"
MEM="768mb"
# TYPES="os:false ss:true"
# COMMAND=" cafe"
COMMAND=" ./top_dilepton_me/bin/top_dilepton_me_x"

# RCP="-rcp=top_dilepton_me/rcp/mumu/mwt_wmt_mumu.rcp"
RCP="-rcp=top_dilepton_me/rcp/mumu/me_mt_mumu.rcp"
NPERENS=40

FINAL_STATE="mumu"
# SUFFIX="sigonly"
SUFFIX="sigbkg"

## normalization
${COMMAND} ${RCP} -cmd=norm -in_samples=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_norm_all_me -out_root=tt_${FINAL_STATE}_norm_all_me_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_norm_all_me_${SUFFIX}.txt -title=tt_${FINAL_STATE}_norm_all_me_${SUFFIX}

## ftop ensemble test

# ${COMMAND} ${RCP} -cmd=me_ft_ens -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ens -out_root=tt_${FINAL_STATE}_ft_ens_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ft_ens_all_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SUFFIX}

## Ensemble test without systematics

# ${COMMAND} ${RCP} -cmd=me_ens -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ens -out_root=tt_${FINAL_STATE}_ens_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SUFFIX}
