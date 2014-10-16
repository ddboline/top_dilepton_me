#!/bin/bash
# . /etc/bashrc
# . setup_command.sh 

COMMAND=" ./top_dilepton_me/bin/top_dilepton_me_x"

RCP="-rcp=top_dilepton_me/rcp/mwt_base/mwt_wmt_comb.rcp"
NPERENS=40

# FINAL_STATE="mumu"
# SUFFIX="sigonly"
SUFFIX="sigbkg"
SUFFIX_FINAL="${SUFFIX}"
# SUFFIX_FINAL="${SUFFIX}_ee_emu_mumu"

FINAL_STATES="ee emu mumu etrk mutrk"
# FINAL_STATES="ee emu mumu mutrk"
# FINAL_STATES="ee emu mumu"

rm comb_ens_sample.txt
for FINAL_STATE in ${FINAL_STATES};
do
    echo ensemble combination tt_${FINAL_STATE}_ens_all_${SUFFIX}.txt >> comb_ens_sample.txt
    echo data     tt_${FINAL_STATE}_data_mwt_${SUFFIX}.txt >> comb_ens_sample.txt
done

${COMMAND} ${RCP} -cmd=comb_data -in_sample=comb_ens_sample.txt -name=comb_all_data -out_root=tt_comb_data_${SUFFIX_FINAL}.root -out_ascii=tt_comb_all_data.txt -nensembles=1000 -title=comb_all_data > data_values_comb_${SUFFIX_FINAL}.txt

# ${COMMAND} ${RCP} -cmd=comb_ens -in_sample=comb_ens_sample.txt -name=comb_all_ens -out_root=tt_comb_ens_all_${SUFFIX_FINAL}.root -out_ascii=tt_comb_all_ens_${SUFFIX_FINAL}.txt -nensembles=1000 -title=comb_all_ens
# 
# ### Correlated systematics
# for SYSTEMATIC in lumi alp_signal color_recomb_dw bjesp_migration bjesn_migration jesp_migration jesn_migration btagbp btagbn btagcp btagcn btaglp btagln bfragn bfragp bfragbn bfragbp bfragcn bfragcp gluonradn gluonradp temp0 temp1 temp2 temp3 emcorrn emcorrp emcorr_preseln emcorr_preselp muonidn muonidp muonison muonisop muontrackn muontrackp tracktrackn tracktrackp em_resn em_resp muon_resn muon_resp em_scalen em_scalep em_offsetn em_offsetp `awk 'BEGIN{for(i=0;i<=40;i++) ITEM=ITEM" pdf_"i; print ITEM}'` `awk 'BEGIN{for(i=0;i<=40;i++) ITEM=ITEM" stat_"i; print ITEM}'`;
# for SYSTEMATIC in alp_signal;
# do
#     rm comb_ens_sample_${SYSTEMATIC}.txt
#     for FINAL_STATE in ${FINAL_STATES};
#     do
#         echo ensemble combination tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.txt >> comb_ens_sample_${SYSTEMATIC}.txt
#     done
#     ${COMMAND} ${RCP} -cmd=comb_ens -rcp=top_dilepton_me/rcp/mwt_wmt_${SYSTEMATIC}.rcp -in_sample=comb_ens_sample_${SYSTEMATIC}.txt -name=comb_all_ens_${SYSTEMATIC}_${SUFFIX} -out_root=tt_comb_ens_all_${SYSTEMATIC}_${SUFFIX}.root -out_ascii=tt_comb_ens_all_${SYSTEMATIC}_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_comb_ens_all_${SYSTEMATIC}_${SUFFIX};
# done

# ### Uncorrelated systematics
# for SYSTEMATIC in sb_up sb_down pyt_bkg triggern triggerp;
# do
#     for FS in ${FINAL_STATES};
#     do
#       rm comb_ens_sample_${SYSTEMATIC}_${FS}.txt
#       for FINAL_STATE in ${FINAL_STATES};
#       do
#           echo $FS $FINAL_STATE
#           if [ "$FS" == "$FINAL_STATE" ]; then
#             echo ensemble combination tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.txt >> comb_ens_sample_${SYSTEMATIC}_${FS}.txt
#           else
#             echo ensemble combination tt_${FINAL_STATE}_ens_all_${SUFFIX}.txt >> comb_ens_sample_${SYSTEMATIC}_${FS}.txt
#           fi
#       done
#       ${COMMAND} ${RCP} -cmd=comb_ens -rcp=top_dilepton_me/rcp/mwt_wmt_${SYSTEMATIC}.rcp -in_sample=comb_ens_sample_${SYSTEMATIC}_${FS}.txt -name=comb_all_ens_${SYSTEMATIC}_${FS}_${SUFFIX} -out_root=tt_comb_ens_all_${SYSTEMATIC}_${FS}_${SUFFIX}.root -out_ascii=tt_comb_ens_all_${SYSTEMATIC}_${FS}_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_comb_ens_all_${SYSTEMATIC}_${FS}_${SUFFIX};
#     done
# done
# 
# 
