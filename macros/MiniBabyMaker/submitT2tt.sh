#!/bin/bash -f


echo "Submitting merged_T2tt_mStop-500to650_mLSP-0to225"
root -b doAll.C\(\"merged_T2tt_mStop-500to650_mLSP-0to225\"\)    > output_V00-03-06/merged_T2tt_mStop-500to650_mLSP-0to225.log 2>&1 &
sleep 20

echo "Submitting merged_T2tt_mStop-500to650_mLSP-250to550"
root -b doAll.C\(\"merged_T2tt_mStop-500to650_mLSP-250to550\"\)  > output_V00-03-06/merged_T2tt_mStop-500to650_mLSP-250to550.log 2>&1 &
sleep 20

echo "Submitting merged_T2tt_mStop-150to350_mLSP-0to250"
root -b doAll.C\(\"merged_T2tt_mStop-150to350_mLSP-0to250\"\)    > output_V00-03-06/merged_T2tt_mStop-150to350_mLSP-0to250.log 2>&1 &
sleep 20

echo "Submitting merged_T2tt_mStop-150to475_mLSP-1"
root -b doAll.C\(\"merged_T2tt_mStop-150to475_mLSP-1\"\)         > output_V00-03-06/merged_T2tt_mStop-150to475_mLSP-1.log 2>&1 &
sleep 20

echo "Submitting merged_T2tt_mStop-500to800_mLSP-1"
root -b doAll.C\(\"merged_T2tt_mStop-500to800_mLSP-1\"\)         > output_V00-03-06/merged_T2tt_mStop-500to800_mLSP-1.log 2>&1 &
sleep 20

echo "Submitting merged_T2tt_mStop-375to475_mLSP-0to375"
root -b doAll.C\(\"merged_T2tt_mStop-375to475_mLSP-0to375\"\)    > output_V00-03-06/merged_T2tt_mStop-375to475_mLSP-0to375.log 2>&1 &
sleep 20

echo "Submitting merged_T2tt_mStop-675to800_mLSP-0to275"
root -b doAll.C\(\"merged_T2tt_mStop-675to800_mLSP-0to275\"\)    > output_V00-03-06/merged_T2tt_mStop-675to800_mLSP-0to275.log 2>&1 &
sleep 20

echo "Submitting merged_T2tt_mStop-675to800_mLSP-300to700"
root -b doAll.C\(\"merged_T2tt_mStop-675to800_mLSP-300to700\"\)  > output_V00-03-06/merged_T2tt_mStop-675to800_mLSP-300to700.log 2>&1 &

