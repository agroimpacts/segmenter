AMID="ami id here"
ITYPE="t2.xlarge"
KEYNAME="aws key name"
SECURITY="aws security group name"
INAME="instance name (not id, should be named instance)"
OWNER="aws account number"
SDASIZE="100"
IAM="activemapper_planet_readwriteS3"
VALIDUNTIL="2020-03-01T23:00:00"  # update date here
MAXPRICE="0.08"  # set max price for spot instance
aws ec2 run-instances --image-id $AMID --count 1 --instance-type $ITYPE \
--iam-instance-profile Name=$IAM --key-name $KEYNAME --security-groups $SECURITY \
--monitoring Enabled=true --block-device-mappings \
"[ { \"DeviceName\": \"/dev/sda1\", \"Ebs\": { \"VolumeSize\": $SDASIZE } } ]" \
--tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value='$INAME'}]' \
'ResourceType=volume,Tags=[{Key=Owner,Value='$OWNER'}]' \
--instance-market-options 'MarketType=spot, SpotOptions={MaxPrice='$MAXPRICE',SpotInstanceType=persistent,ValidUntil='$VALIDUNTIL', InstanceInterruptionBehavior=stop}'
