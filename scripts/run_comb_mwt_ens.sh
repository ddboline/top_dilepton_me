#!/bin/bash
# . /etc/bashrc
# . setup_command.sh 

TIME="72:00:00"
MEM="768mb"
# TYPES="os:false ss:true"
# COMMAND=" cafe"
COMMAND="echo runme -clued0 -cput=${TIME} -mem=${MEM} -outdir=/rooms/snapper/projects/ddboline/results_mwt/mwt_weights/mumu -tar=/rooms/snapper/projects/ddboline/me_p18.10.00.tar.gz"
COMMAND="./top_dilepton_me/bin/top_dilepton_me_x"

RCP="-rcp=top_dilepton_me/rcp/mwt_base/mwt_wmt_comb.rcp"
NPERENS=40

# FINAL_STATE="mumu"
# SUFFIX="sigonly"
SUFFIX="sigbkg"
SUFFIX_FINAL="${SUFFIX}"
# SUFFIX_FINAL="${SUFFIX}_ee_emu_mumu"

rm comb_ens_sample.txt
for FINAL_STATE in ee emu mumu etrk mutrk;
# for FINAL_STATE in ee emu mumu;
do
    echo ensemble combination tt_${FINAL_STATE}_ens_all_${SUFFIX}.txt >> comb_ens_sample.txt
#     echo data     tt_${FINAL_STATE}_data_mwt_${SUFFIX}.txt >> comb_ens_sample.txt
done

${COMMAND} ${RCP} -cmd=comb_ens -in_sample=comb_ens_sample.txt -name=comb_all_ens -out_root=tt_comb_ens_all_${SUFFIX_FINAL}.root -out_ascii=tt_comb_all_ens_${SUFFIX_FINAL}.txt -nensembles=1000 -title=comb_all_ens
