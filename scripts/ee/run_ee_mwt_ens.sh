#!/bin/bash
# . /etc/bashrc

# . setup_command.sh 

TIME="72:00:00"
MEM="768mb"
# TYPES="os:false ss:true"
# COMMAND=" cafe"
COMMAND="echo runme -cabsrv1 -cput=${TIME} -mem=${MEM} -outdir=/rooms/snapper/projects/ddboline/results_mwt/mwt_weights/ee -tar=/rooms/snapper/projects/ddboline/me_p18.10.00.tar.gz"
COMMAND=" ./top_dilepton_me/bin/top_dilepton_me_x"

RCP="-rcp=top_dilepton_me/rcp/ee/mwt_wmt_ee.rcp"
NPERENS=40

FINAL_STATE="ee"
# SUFFIX="sigonly"
SUFFIX="sigbkg"

## Ensemble test without systematics

${COMMAND} ${RCP} -cmd=mwt_ens -in_sample=samples_mwt_p20/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ens -out_root=tt_${FINAL_STATE}_ens_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SUFFIX}

# nohup sh top_dilepton_me/scripts/ee/run_ee_mwt_ens.sh > cafe_ee_ens.out 2> cafe_ee_ens.err &
# nohup sh top_dilepton_me/scripts/emu/run_emu_mwt_ens.sh > cafe_emu_ens.out 2> cafe_emu_ens.err &
# nohup sh top_dilepton_me/scripts/mumu/run_mumu_mwt_ens.sh > cafe_mumu_ens.out 2> cafe_mumu_ens.err &
# nohup sh top_dilepton_me/scripts/etrk/run_etrk_mwt_ens.sh > cafe_etrk_ens.out 2> cafe_etrk_ens.err &
# nohup sh top_dilepton_me/scripts/mutrk/run_mutrk_mwt_ens.sh > cafe_mutrk_ens.out 2> cafe_mutrk_ens.err &
# nohup sh top_dilepton_me/scripts/run_comb_mwt_ens.sh > cafe_comb_ens.out 2> cafe_comb_ens.err &
