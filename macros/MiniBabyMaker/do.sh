XPATH=/nfs-7/userdata/stop/output_V00-02-21_2012_4jskim/
#XPATH=/nfs-7/userdata/stop/output_V00-02-21_2012_4jskim/altttbar

#XPATH=/nfs-3/userdata/stop/Train/V00-02-20__V00-03-01_4jetsMET50_bkg/
#XPATH=/nfs-3/userdata/stop/Train/V00-02-20__V00-03-01_4jetsMET50_bkg/altttbar
#XPATH=/nfs-3/userdata/stop/Train/V00-02-18__V00-03-01_4jetsMET50_T2tt/

#XPATH=/nfs-3/userdata/stop/cms2V05-03-25_stoplooperV00-02-18/T2tt_mad/
#XPATH=/nfs-3/userdata/stop/cms2V05-03-25_stoplooperV00-02-18/T2bw_coarse/


#XPATH=/nfs-3/userdata/stop/cms2V05-03-18_stoplooperV00-02-07/crabT2bw_3/

for tag in `ls -1 $XPATH | grep root |  cut -d'.' -f1`;
#do echo root -b -q -l doFile.C\(\"$XPATH\",\"$tag\"\)
 do nohup root -b -q -l doFile.C\(\"${XPATH}\",\"${tag}\"\) > output/log/${tag}.log &
#  sleep 10
done

#let "FAIL=0"

#for job in `jobs -p`
#do echo $job
#   wait $job || let "FAIL+=1"
#done

#echo $FAIL

#if [ "$FAIL" == "0" ];
# then echo "YAY!"
# else echo "FAIL! ($FAIL)"
#fi
