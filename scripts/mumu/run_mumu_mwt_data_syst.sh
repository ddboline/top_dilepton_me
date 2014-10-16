#!/bin/bash
# . /etc/bashrc

# . setup_command.sh 

COMMAND=" ./top_dilepton_me/bin/top_dilepton_me_x"

RCP="-rcp=top_dilepton_me/rcp/mumu/mwt_wmt_mumu.rcp"
NPERENS=40

FINAL_STATE="mumu"
# SUFFIX="sigonly"
SUFFIX="sigbkg"

## obligatory run over the data

${COMMAND} ${RCP} -cmd=mwt_data -in_sample=samples_mwt_p20/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_data_mwt -out_root=tt_${FINAL_STATE}_data_mwt_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_data_mwt_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_data_mwt_${SUFFIX} > data_values_${FINAL_STATE}.txt

# ## Ensemble test without systematics
# 
# ${COMMAND} ${RCP} -cmd=mwt_ens -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ens -out_root=tt_${FINAL_STATE}_ens_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SUFFIX}
# 
# ## S/B systematics -- different naming convention
# 
# for SYSTEMATIC in sb_up sb_down;
# do
#     ${COMMAND} ${RCP} -cmd=mwt_ens -rcp=top_dilepton_me/rcp/systematics/mwt_wmt_${FINAL_STATE}_${SYSTEMATIC}.rcp -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ens_${SYSTEMATIC}_${SUFFIX} -out_root=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX};
# done
# 
# ### Resolution systematics -- different sample files, no efficiency correction
# for SYSTEMATIC in em_resn em_resp muon_resn muon_resp em_scalen em_scalep em_offsetn em_offsetp;
# do
#     ${COMMAND} ${RCP} -cmd=mwt_ens -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}_${SYSTEMATIC}.txt -name=${FINAL_STATE}_all_ens_${SYSTEMATIC} -out_root=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}
# done
# 
# ### Shape systematics -- different sample files
# 
# for SYSTEMATIC in pyt_bkg alp_signal color_recomb_dw;
# for SYSTEMATIC in alp_signal;
# do
#     ${COMMAND} ${RCP} -cmd=ftop_mwt -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}_${SYSTEMATIC}.txt -name=${FINAL_STATE}_all_ftop -out_root=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.txt -title=tt_${FINAL_STATE}_ftop_all_${SUFFIX} > ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt
#     echo ${SYSTEMATIC}
#     echo ftop_true `awk '/mass 175/ && /ftop_exp/ {if($6>0 && $6<=1)FTOP=FTOP" "$6 ; if($6<0) FTOP=FTOP" 0.01" ; if($6>1) FTOP=FTOP" 1.0"} END{print FTOP}' ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt`
#     echo ftop_true `awk '/mass 175/ && /ftop_exp/ {if($6>0 && $6<=1)FTOP=FTOP" "$6 ; if($6<0) FTOP=FTOP" 0.01" ; if($6>1) FTOP=FTOP" 1.0"} END{print FTOP}' ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt` > temp_sb_${FINAL_STATE}.rcp
# 
#     ${COMMAND} ${RCP} -cmd=mwt_ens -rcp=temp_sb_${FINAL_STATE}.rcp -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}_${SYSTEMATIC}.txt -name=${FINAL_STATE}_all_ens_${SYSTEMATIC} -out_root=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}
# done

# ### jes , bjes
# 
# for SYSTEMATIC in bjesp_migration bjesn_migration jesp_migration jesn_migration;
# do
#     SYSTEMATIC_SHORT=`echo $SYSTEMATIC | sed 's:_migration::g'`
#     ${COMMAND} ${RCP} -cmd=ftop_mwt -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}_${SYSTEMATIC}.txt -name=${FINAL_STATE}_all_ftop -out_root=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.txt -title=tt_${FINAL_STATE}_ftop_all_${SUFFIX} > ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt
#     echo ${SYSTEMATIC}
#     echo ftop_true `awk '/mass 170/ && /ftop_exp/ {if($6>0 && $6<=1)FTOP=FTOP" "$6 ; if($6<0) FTOP=FTOP" 0.01" ; if($6>1) FTOP=FTOP" 1.0"} END{print FTOP}' ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt`
#     echo ftop_true `awk '/mass 170/ && /ftop_exp/ {if($6>0 && $6<=1)FTOP=FTOP" "$6 ; if($6<0) FTOP=FTOP" 0.01" ; if($6>1) FTOP=FTOP" 1.0"} END{print FTOP}' ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt` > temp_sb_${FINAL_STATE}.rcp
# 
#     ${COMMAND} ${RCP} -cmd=mwt_ens -rcp=temp_sb_${FINAL_STATE}.rcp -rcp=top_dilepton_me/rcp/systematics/mwt_wmt_${SYSTEMATIC}.rcp -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}_${SYSTEMATIC_SHORT}.txt -name=${FINAL_STATE}_all_ens_${SYSTEMATIC} -out_root=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}
# done
# 
# 
# ### btag , bfrag , gluon radiation , template , efficiency systematics
# 
# for SYSTEMATIC in btagbp btagbn btagcp btagcn btaglp btagln bfragn bfragp bfragbn bfragbp bfragcn bfragcp gluonradn gluonradp temp0 temp1 temp2 temp3;
# do
#     ${COMMAND} ${RCP} -cmd=ftop_mwt -rcp=top_dilepton_me/rcp/systematics/mwt_wmt_${SYSTEMATIC}.rcp -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ftop -out_root=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.txt -title=tt_${FINAL_STATE}_ftop_all_${SUFFIX} > ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt
#     echo ${SYSTEMATIC}
#     echo ftop_true `awk '/mass 170/ && /ftop_exp/ {if($6>0 && $6<=1)FTOP=FTOP" "$6 ; if($6<0) FTOP=FTOP" 0.01" ; if($6>1) FTOP=FTOP" 1.0"} END{print FTOP}' ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt`
#     echo ftop_true `awk '/mass 170/ && /ftop_exp/ {if($6>0 && $6<=1)FTOP=FTOP" "$6 ; if($6<0) FTOP=FTOP" 0.01" ; if($6>1) FTOP=FTOP" 1.0"} END{print FTOP}' ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt` > temp_sb_${FINAL_STATE}.rcp
# 
#     ${COMMAND} ${RCP} -cmd=mwt_ens -rcp=temp_sb_${FINAL_STATE}.rcp -rcp=top_dilepton_me/rcp/systematics/mwt_wmt_${SYSTEMATIC}.rcp -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ens_${SYSTEMATIC}_${SUFFIX} -out_root=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX};
# done
# 
# ### Efficiency Systematics
# for SYSTEMATIC in emcorrn emcorrp emcorr_preseln emcorr_preselp triggern triggerp muonidn muonidp muonison muonisop muontrackn muontrackp tracktrackn tracktrackp lumi;
# for SYSTEMATIC in lumi;
# do
#     ${COMMAND} ${RCP} -cmd=ftop_mwt -rcp=top_dilepton_me/rcp/systematics/mwt_wmt_${SYSTEMATIC}_${FINAL_STATE}.rcp -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ftop -out_root=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.txt -title=tt_${FINAL_STATE}_ftop_all_${SUFFIX} > ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt
#     echo ${SYSTEMATIC}
#     echo ftop_true `awk '/mass 170/ && /ftop_exp/ {if($6>0 && $6<=1)FTOP=FTOP" "$6 ; if($6<0) FTOP=FTOP" 0.01" ; if($6>1) FTOP=FTOP" 1.0"} END{print FTOP}' ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt`
#     echo ftop_true `awk '/mass 170/ && /ftop_exp/ {if($6>0 && $6<=1)FTOP=FTOP" "$6 ; if($6<0) FTOP=FTOP" 0.01" ; if($6>1) FTOP=FTOP" 1.0"} END{print FTOP}' ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt` > temp_sb_${FINAL_STATE}.rcp
# 
#     ${COMMAND} ${RCP} -cmd=mwt_ens -rcp=temp_sb_${FINAL_STATE}.rcp -rcp=top_dilepton_me/rcp/systematics/mwt_wmt_${SYSTEMATIC}_${FINAL_STATE}.rcp -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ens_${SYSTEMATIC}_${SUFFIX} -out_root=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX};
# done

# ### PDF systematics
# 
# for NUM in `awk 'BEGIN{for(i=0;i<=40;i++) ITEM=ITEM" "i; print ITEM}'`;
# do
#     echo pdf_syst            ${NUM} > temp_pdf_syst_${FINAL_STATE}.rcp
#     SYSTEMATIC="pdf_${NUM}"
# 
#     ${COMMAND} ${RCP} -cmd=ftop_mwt -rcp=temp_pdf_syst.rcp -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ftop -out_root=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.txt -title=tt_${FINAL_STATE}_ftop_all_${SUFFIX} > ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt
#     echo ${SYSTEMATIC}
#     echo ftop_true `awk '/mass 170/ && /ftop_exp/ {if($6>0 && $6<=1)FTOP=FTOP" "$6 ; if($6<0) FTOP=FTOP" 0.01" ; if($6>1) FTOP=FTOP" 1.0"} END{print FTOP}' ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt`
#     echo ftop_true `awk '/mass 170/ && /ftop_exp/ {if($6>0 && $6<=1)FTOP=FTOP" "$6 ; if($6<0) FTOP=FTOP" 0.01" ; if($6>1) FTOP=FTOP" 1.0"} END{print FTOP}' ftop_values_${FINAL_STATE}_${SYSTEMATIC}.txt` > temp_pdf_syst_${FINAL_STATE}.rcp
# 
#     ${COMMAND} ${RCP} -cmd=mwt_ens -rcp=temp_pdf_syst_${FINAL_STATE}.rcp -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ens_${SYSTEMATIC}_${SUFFIX} -out_root=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX};
# done
# 
# ### Statistical systematic
# 
# for NUM in `awk 'BEGIN{for(i=0;i<=40;i++) ITEM=ITEM" "i; print ITEM}'`;
# do
#     SYSTEMATIC="stat_${NUM}"
#     ${COMMAND} ${RCP} -cmd=mwt_ens -in_sample=samples_mwt/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ens_${SYSTEMATIC}_${SUFFIX} -out_root=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX}.txt -perensemble=${NPERENS} -nensembles=1000 -title=tt_${FINAL_STATE}_ens_all_${SYSTEMATIC}_${SUFFIX};
# done
