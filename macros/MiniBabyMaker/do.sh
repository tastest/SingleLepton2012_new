#XPATH=/nfs-3/userdata/stop/output_V00-02-16_2012_4jskim/
XPATH=/nfs-3/userdata/stop/output_V00-02-18_2012_4jskim/
#XPATH=/nfs-3/userdata/stop/output_V00-02-16_2012/

#XPATH=/nfs-3/userdata/stop/cms2V05-03-18_stoplooperV00-02-07/crabT2tt_3
#XPATH=/nfs-3/userdata/stop/cms2V05-03-18_stoplooperV00-02-07/crabT2bw_3/

#XPATH=/nfs-3/userdata/stop/MiniBabies/V00-02-16_2012_4jskim

for tag in `ls -1 $XPATH | grep root |  cut -d'.' -f1`;
  do  root -b -q -l doFile.C\(\"$XPATH\",\"$tag\"\) | tee $tag.log &
#  sleep 10
done

let "FAIL=0"

for job in `jobs -p`
do echo $job
   wait $job || let "FAIL+=1"
done

echo $FAIL

if [ "$FAIL" == "0" ];
 then echo "YAY!"
 else echo "FAIL! ($FAIL)"
fi

