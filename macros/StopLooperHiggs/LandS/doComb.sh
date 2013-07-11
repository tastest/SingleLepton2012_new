#! /bin/bash

function process_card
{
cardname=$1
# card.txt
# lands.exe
/home/users/vimartin/CMSSW/CMSSW_5_3_2_patch4/LandS/test/lands.exe -d ${cardname}.txt -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --rMin 0 --rMax 1 --freq --nToysForCLsb 10000 --nToysForCLb 5000 --seed 123 -n $cardname &
#/home/users/vimartin/CMSSW/CMSSW_5_3_2_patch4/LandS/test/lands.exe -d ${cardname}.txt -M Hybrid --freq --ExpectationHints Asymptotic --scanRs 1 --rMin 0 --rMax 20 --freq --nToysForCLsb 10000 --nToysForCLb 5000 --seed 123 -n $cardname &
}

process_card m350
process_card m450

