#!/bin/bash
# . /etc/bashrc

# . setup_command.sh 

TIME="72:00:00"
MEM="768mb"
# TYPES="os:false ss:true"
# COMMAND=" cafe"
COMMAND=" ./top_dilepton_me/bin/top_dilepton_me_x"

RCP="-rcp=top_dilepton_me/rcp/ee/mwt_wmt_ee.rcp"
NPERENS=40

FINAL_STATE="ee"
# SUFFIX="sigonly"
SUFFIX="sigbkg"

### Get Systematic output

for SYSTEMATIC in "lumi,lumi" "alp_signal,color_recomb_dw" "bjesn_migration,bjesp_migration" "jesn_migration,jesp_migration" "btagbp,btagbn" "btagcp,btagcn" "btaglp,btagln" "gluonradn,gluonradp" "bfragn,bfragp" "bfragbn,bfragbp" "bfragcn,bfragcp" "emcorrn,emcorrp" "emcorr_preseln,emcorr_preselp" "muonidn,muonidp" "muonison,muonisop" "muontrackn,muontrackp" "tracktrackn,tracktrackp" "sb_up,sb_down" "pyt_bkg,pyt_bkg" "triggern,triggerp";
do 
  PDF1=`echo ${SYSTEMATIC} | awk '/./ {split($1,A,","); print A[1]}'` PDF2=`echo ${SYSTEMATIC} | awk '/./ {split($1,A,","); print A[2]}'`;
  python get_systematic.py tt_${FINAL_STATE}_ens_all_sigbkg.root tt_${FINAL_STATE}_ens_all_${PDF1}_sigbkg.root tt_${FINAL_STATE}_ens_all_${PDF2}_sigbkg.root ${FINAL_STATE}_${PDF1} > temp.txt;
  echo $PDF1 `awk '/cal_0/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_0/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_1/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_1/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_2/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_2/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_3/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_3/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_4/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_4/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_5/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_5/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_6/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_6/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_7/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_7/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt`;
done

echo "" > temp.txt
for SYSTEMATIC in "temp0,temp1" "temp2,temp3";
do 
  PDF1=`echo ${SYSTEMATIC} | awk '/./ {split($1,A,","); print A[1]}'` PDF2=`echo ${SYSTEMATIC} | awk '/./ {split($1,A,","); print A[2]}'`;
  python get_systematic.py tt_${FINAL_STATE}_ens_all_sigbkg.root tt_${FINAL_STATE}_ens_all_${PDF1}_sigbkg.root tt_${FINAL_STATE}_ens_all_${PDF2}_sigbkg.root ${FINAL_STATE}_${PDF1} >> temp.txt;
done
echo $PDF1 `awk '/cal_0/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_0/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_1/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_1/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_2/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_2/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_3/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_3/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_4/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_4/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_5/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_5/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_6/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_6/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_7/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_7/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt`;

echo "" > temp.txt
for SYSTEMATIC in `awk 'BEGIN{for(i=0;i<=40;i+=2) ITEM=ITEM" pdf_"i",pdf_"i+1; print ITEM}'`;
do 
  PDF1=`echo ${SYSTEMATIC} | awk '/./ {split($1,A,","); print A[1]}'` PDF2=`echo ${SYSTEMATIC} | awk '/./ {split($1,A,","); print A[2]}'`;
  python get_systematic.py tt_${FINAL_STATE}_ens_all_sigbkg.root tt_${FINAL_STATE}_ens_all_${PDF1}_sigbkg.root tt_${FINAL_STATE}_ens_all_${PDF2}_sigbkg.root ${FINAL_STATE}_${PDF1} >> temp.txt;
done
echo $PDF1 `awk '/cal_0/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_0/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_1/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_1/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_2/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_2/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_3/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_3/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_4/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_4/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_5/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_5/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_6/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_6/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_7/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_7/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt`;

echo "" > temp.txt
for SYSTEMATIC in `awk 'BEGIN{for(i=0;i<=40;i+=2) ITEM=ITEM" stat_"i",stat_"i+1; print ITEM}'`;
do 
  PDF1=`echo ${SYSTEMATIC} | awk '/./ {split($1,A,","); print A[1]}'` PDF2=`echo ${SYSTEMATIC} | awk '/./ {split($1,A,","); print A[2]}'`;
  python get_systematic.py tt_${FINAL_STATE}_ens_all_sigbkg.root tt_${FINAL_STATE}_ens_all_${PDF1}_sigbkg.root tt_${FINAL_STATE}_ens_all_${PDF2}_sigbkg.root ${FINAL_STATE}_${PDF1} >> temp.txt; 
done
echo $PDF1 `awk '/cal_0/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_0/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_1/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_1/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_2/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_2/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_3/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_3/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_4/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_4/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_5/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_5/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_6/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_6/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_7/ {DOWN+=sqrt($3**2);ITS+=1} END{print DOWN/ITS}' temp.txt` `awk '/cal_7/ {DOWN+=sqrt($4**2);ITS+=1} END{print DOWN/ITS}' temp.txt`;

# nohup sh top_dilepton_me/scripts/ee/run_ee_mwt_syst_final.sh > cafe_ee_ens_syst.out 2> cafe_ee_ens_syst.err &
# nohup sh top_dilepton_me/scripts/emu/run_emu_mwt_syst_final.sh > cafe_emu_ens_syst.out 2> cafe_emu_ens_syst.err &
# nohup sh top_dilepton_me/scripts/mumu/run_mumu_mwt_syst_final.sh > cafe_mumu_ens_syst.out 2> cafe_mumu_ens_syst.err &
# nohup sh top_dilepton_me/scripts/etrk/run_etrk_mwt_syst_final.sh > cafe_etrk_ens_syst.out 2> cafe_etrk_ens_syst.err &
# nohup sh top_dilepton_me/scripts/mutrk/run_mutrk_mwt_syst_final.sh > cafe_mutrk_ens_syst.out 2> cafe_mutrk_ens_syst.err &
# nohup sh top_dilepton_me/scripts/run_comb_mwt_syst_final.sh > cafe_comb_ens_syst.out 2> cafe_comb_ens_syst.err &
