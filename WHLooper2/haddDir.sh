#!/bin/bash

if [ -z "$1" ]
  then
    echo "usage: haddDir.sh <dir>"
fi

# first make wlight dd background
root -b -q -l makeWLightBg.C\(\"$1\"\)

#hadd -f $1/ttbar_histos.root $1/ttbar_0l_histos.root $1/ttbar_1l_histos.root $1/ttbar_2l_histos.root 
#hadd -f $1/ttbar_histos.root $1/ttbar_1l_histos.root $1/ttbar_2l_histos.root 
#hadd -f $1/ttbar_mg_histos.root $1/ttbar_mg_1l_histos.root $1/ttbar_mg_2l_histos.root 
#hadd -f $1/single_top_histos.root $1/single_top_1l_histos.root $1/single_top_2l_histos.root 
#hadd -f $1/wjets_histos.root $1/wjets_nobb_histos.root $1/wbb_histos.root 
#hadd -f $1/wjets_histos.root $1/wjets_nobb_histos.root $1/wjets_onlybb_histos.root 
hadd -f $1/wbb_all_histos.root $1/wbb_histos.root $1/wzbb_histos.root 
hadd -f $1/rare_histos.root $1/VV_histos.root $1/zjets_histos.root $1/ttV_histos.root $1/VVV_histos.root
hadd -f $1/rare_all_histos.root $1/rare_histos.root $1/whbb_histos.root
#hadd -f $1/allbg_histos.root $1/ttbar_mg_histos.root $1/wjets_comb_histos.root $1/single_top_histos.root $1/rare_histos.root 
#hadd -f $1/allbg_histos.root $1/ttbar_mg_histos.root $1/wjets_comb_histos.root $1/single_top_histos.root $1/rare_histos.root $1/whbb_histos.root
hadd -f $1/allbg_histos.root $1/ttbar_mg_1l_histos.root $1/ttbar_mg_2l_histos.root $1/wjets_nobb_histos.root $1/wbb_all_histos.root $1/single_top_1l_histos.root $1/single_top_2l_histos.root $1/rare_histos.root $1/whbb_histos.root
#hadd -f $1/allbg_histos.root $1/ttbar_mg_histos.root $1/wjets_histos.root $1/single_top_histos.root $1/rare_histos.root 
hadd -f $1/allbg_dd_histos.root $1/ttbar_mg_1l_histos.root $1/ttbar_mg_2l_histos.root $1/wlight_dd_histos.root $1/wbb_all_histos.root $1/single_top_1l_histos.root $1/single_top_2l_histos.root $1/rare_histos.root $1/whbb_histos.root 

hadd -f $1/bg_1lb_histos.root $1/ttbar_mg_1l_histos.root $1/wbb_all_histos.root $1/single_top_1l_histos.root 
hadd -f $1/top1l_histos.root $1/ttbar_mg_1l_histos.root $1/single_top_1l_histos.root 
hadd -f $1/bg_others_histos.root $1/ttbar_mg_2l_histos.root $1/wjets_nobb_histos.root $1/single_top_2l_histos.root $1/rare_histos.root $1/whbb_histos.root
hadd -f $1/bg_others_dd_histos.root $1/ttbar_mg_2l_histos.root $1/wlight_dd_histos.root $1/single_top_2l_histos.root $1/rare_histos.root $1/whbb_histos.root
hadd -f $1/diltop_histos.root $1/ttbar_mg_2l_histos.root $1/single_top_2l_histos.root 


#hadd -f $1/allbg_histos.root $1/ttbar_mg_histos.root $1/wjets_comb_histos.root $1/single_top_histos.root $1/rare_histos.root $1/whbb_histos.root
#hadd -f $1/allbg_histos.root $1/ttbar_mg_histos.root $1/wjets_histos.root $1/single_top_histos.root $1/rare_histos.root $1/whbb_histos.root
#hadd -f $1/allbg_histos.root $1/ttbar_histos.root $1/wjets_histos.root $1/single_top_histos.root $1/rare_histos.root $1/whbb_histos.root
#hadd -f $1/allbg_histos.root $1/ttbar_histos.root $1/wjets_histos.root $1/single_top_histos.root $1/rare_histos.root 