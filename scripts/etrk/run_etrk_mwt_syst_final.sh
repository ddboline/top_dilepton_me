#!/bin/bash
# . /etc/bashrc

# . setup_command.sh 

TIME="72:00:00"
MEM="768mb"
# TYPES="os:false ss:true"
# COMMAND=" cafe"
COMMAND=" ./top_dilepton_me/bin/top_dilepton_me_x"

RCP="-rcp=top_dilepton_me/rcp/etrk/mwt_wmt_etrk.rcp"
NPERENS=40

FINAL_STATE="etrk"
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
