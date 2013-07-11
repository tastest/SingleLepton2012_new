#! /bin/bash

function run_limit
{
cardname=$1

echo $cardname > limitlogs/${cardname}_log.txt
echo $cardname >> limitlogs/Summary_log.txt

# root -b -q "getExpectedLimits.C (\"$cardname\")" | grep "limit:" >> limitlogs/${cardname}_log.txt &
# root -b -q "getExpectedLimits.C (\"$cardname\")" | grep "limit:" >> limitlogs/Summary_log.txt

root -b -q "getExpectedLimits.C (\"$cardname\")" | grep -i "obs" >> limitlogs/${cardname}_log.txt &
root -b -q "getExpectedLimits.C (\"$cardname\")" | grep -i "obs" >> limitlogs/Summary_log.txt

echo "" >> limitlogs/Summary_log.txt

}

echo "Summary:" > limitlogs/Summary_log.txt

run_limit all
run_limit 1l3b
run_limit 1l4b
run_limit 2l3b
run_limit 2l4b

