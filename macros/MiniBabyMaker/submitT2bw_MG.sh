#!/bin/bash -f

echo "Submitting merged_T2bw_mStop-100to475_mLSP-0to375_x025.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop-100to475_mLSP-0to375_x025\"\)    > output_V00-03-11/merged_T2bw_mStop-100to475_mLSP-0to375_x025.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_mStop-500to800_mLSP-0to700_x025.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop-500to800_mLSP-0to700_x025\"\)    > output_V00-03-11/merged_T2bw_mStop-500to800_mLSP-0to700_x025.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_mStop_100to475_mLSP_0to375_x050.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop_100to475_mLSP_0to375_x050\"\)    > output_V00-03-11/merged_T2bw_mStop_100to475_mLSP_0to375_x050.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_mStop_500to800_mLSP_0to700_x050.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop_500to800_mLSP_0to700_x050\"\)    > output_V00-03-11/merged_T2bw_mStop_500to800_mLSP_0to700_x050.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_mStop-100to475_mLSP-0to375_x075.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop-100to475_mLSP-0to375_x075\"\)    > output_V00-03-11/merged_T2bw_mStop-100to475_mLSP-0to375_x075.log 2>&1 &
sleep 10

echo "Submitting merged_T2bw_mStop-500to800_mLSP-0to700_x075.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop-500to800_mLSP-0to700_x075\"\)    > output_V00-03-11/merged_T2bw_mStop-500to800_mLSP-0to700_x075.log 2>&1 &
sleep 10

echo "merged_T2bw_mStop-500to650_mLSP-150to300_lepFilter_x025.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop-500to650_mLSP-150to300_lepFilter_x025\"\)    > output_V00-03-11/merged_T2bw_mStop-500to650_mLSP-150to300_lepFilter_x025.log 2>&1 &
sleep 10

echo "merged_T2bw_mStop-500to650_mLSP-150to300_lepFilter_x050.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop-500to650_mLSP-150to300_lepFilter_x050\"\)    > output_V00-03-11/merged_T2bw_mStop-500to650_mLSP-150to300_lepFilter_x050.log 2>&1 &
sleep 10

echo "merged_T2bw_mStop-500to650_mLSP-150to300_lepFilter_x075.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop-500to650_mLSP-150to300_lepFilter_x075\"\)    > output_V00-03-11/merged_T2bw_mStop-500to650_mLSP-150to300_lepFilter_x075.log 2>&1 &
sleep 10

echo "merged_T2bw_mStop-100to475_mLSP-0to300_lepFilter_x025.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop-100to475_mLSP-0to300_lepFilter_x025\"\)    > output_V00-03-11/merged_T2bw_mStop-100to475_mLSP-0to300_lepFilter_x025.log 2>&1 &
sleep 10

echo "merged_T2bw_mStop-100to475_mLSP-0to300_lepFilter_x050.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop-100to475_mLSP-0to300_lepFilter_x050\"\)    > output_V00-03-11/merged_T2bw_mStop-100to475_mLSP-0to300_lepFilter_x050.log 2>&1 &
sleep 10

echo "merged_T2bw_mStop-100to475_mLSP-0to300_lepFilter_x075.root"
root -b doAll_ben.C\(\"merged_T2bw_mStop-100to475_mLSP-0to300_lepFilter_x075\"\)    > output_V00-03-11/merged_T2bw_mStop-100to475_mLSP-0to300_lepFilter_x075.log 2>&1 &
sleep 10

