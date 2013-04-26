#!/bin/bash -f


echo "Submitting merged_T2bw_coarse.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse\"\)    > output_V00-03-10/merged_T2bw_coarse.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_coarse_1.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse_1\"\)    > output_V00-03-10/merged_T2bw_coarse_1.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_coarse_2.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse_2\"\)    > output_V00-03-10/merged_T2bw_coarse_2.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_coarse_3.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse_3\"\)    > output_V00-03-10/merged_T2bw_coarse_3.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_coarse_4.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse_4\"\)    > output_V00-03-10/merged_T2bw_coarse_4.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_coarse_5.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse_5\"\)    > output_V00-03-10/merged_T2bw_coarse_5.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_coarse_6.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse_6\"\)    > output_V00-03-10/merged_T2bw_coarse_6.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_coarse_7.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse_7\"\)    > output_V00-03-10/merged_T2bw_coarse_7.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_coarse_8.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse_8\"\)    > output_V00-03-10/merged_T2bw_coarse_8.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_coarse_9.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse_9\"\)    > output_V00-03-10/merged_T2bw_coarse_9.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_coarse_10.root"
root -b doAll_ben.C\(\"merged_T2bw_coarse_10\"\)    > output_V00-03-10/merged_T2bw_coarse_10.log 2>&1 &
sleep 10

