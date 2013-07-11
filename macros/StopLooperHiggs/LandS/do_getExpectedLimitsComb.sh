#! /bin/bash

function run_limit
{
cardname=$1

echo $cardname > limitlogs/${cardname}_log.txt
echo $cardname >> limitlogs/Summary_log_comb.txt

# root -b -q "getExpectedLimits.C (\"$cardname\")" | grep "limit:" >> limitlogs/${cardname}_log.txt &
# root -b -q "getExpectedLimits.C (\"$cardname\")" | grep "limit:" >> limitlogs/Summary_log.txt

root -b -q "getExpectedLimits.C (\"$cardname\")" | grep -i "obs" >> limitlogs/${cardname}_log.txt &
root -b -q "getExpectedLimits.C (\"$cardname\")" | grep -i "obs" >> limitlogs/Summary_log_comb.txt

echo "" >> limitlogs/Summary_log_comb.txt

}

echo "Summary:" > limitlogs/Summary_log_comb.txt

run_limit m350
run_limit m450

