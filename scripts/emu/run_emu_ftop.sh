#!/bin/bash
# . /etc/bashrc

# . setup_command.sh 

COMMAND=" ./top_dilepton_me/bin/top_dilepton_me_x"

RCP="-rcp=top_dilepton_me/rcp/emu/mwt_wmt_emu.rcp"
# RCP="-rcp=top_dilepton_me/rcp/emu/me_mt_emu.rcp"
NPERENS=40

FINAL_STATE="emu"
# SUFFIX="sigonly"
SUFFIX="sigbkg"

## Get ftop values for both nominal and various systematics

${COMMAND} ${RCP} -cmd=ftop_mwt -in_sample=samples_mwt_p20/ensemble_all_${FINAL_STATE}.txt -name=${FINAL_STATE}_all_ftop -out_root=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.root -out_ascii=tt_${FINAL_STATE}_ftop_all_${SUFFIX}.txt -title=tt_${FINAL_STATE}_ftop_all_${SUFFIX} > ftop_values_${FINAL_STATE}.txt

echo expected from MC
echo ftop_input `awk '/mass 170/ && /ftop_exp/ {FTOP=FTOP" "$6} END{print FTOP}' ftop_values_${FINAL_STATE}.txt`
echo ftop_true `awk '/mass 170/ && /ftop_exp/ {FTOP=FTOP" "$6} END{print FTOP}' ftop_values_${FINAL_STATE}.txt`
echo sb up
echo ftop_true `awk '/mass 170/ && /ftop_exp/ {FTOP=FTOP" "$10} END{print FTOP}' ftop_values_${FINAL_STATE}.txt`
echo sb down
echo ftop_true `awk '/mass 170/ && /ftop_exp/ {FTOP=FTOP" "$11} END{print FTOP}' ftop_values_${FINAL_STATE}.txt`

echo from data
echo ftop_input `awk '/mass 170/ && /ftop_data/ {FTOP=FTOP" "$6} END{print FTOP}' ftop_values_${FINAL_STATE}.txt`
echo ftop_true `awk '/mass 170/ && /ftop_data/ {FTOP=FTOP" "$6} END{print FTOP}' ftop_values_${FINAL_STATE}.txt`
echo sb up
echo ftop_true `awk '/mass 170/ && /ftop_data/ {FTOP=FTOP" "$10} END{print FTOP}' ftop_values_${FINAL_STATE}.txt`
echo sb down
echo ftop_true `awk '/mass 170/ && /ftop_data/ {FTOP=FTOP" "$11} END{print FTOP}' ftop_values_${FINAL_STATE}.txt`

echo sig `awk '/mass 170/ && / bkg / {ITEM=ITEM" "$6" "$8} END{print ITEM}' ftop_values_${FINAL_STATE}.txt`
echo bkg `awk '/mass 170/ && / bkg / {ITEM=ITEM" "$10" "$12} END{print ITEM}' ftop_values_${FINAL_STATE}.txt`
echo fake `awk '/mass 170/ && / bkg / {ITEM=ITEM" "$14" "$16} END{print ITEM}' ftop_values_${FINAL_STATE}.txt`
echo data `awk '/mass 170/ && / data / {ITEM=ITEM" "$10" "$12} END{print ITEM}' ftop_values_${FINAL_STATE}.txt`
