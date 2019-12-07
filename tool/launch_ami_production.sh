AMID="***REMOVED***"
ITYPE="t2.xlarge"
KEYNAME="***REMOVED***"
SECURITY="***REMOVED***"
INAME="segmenter_6"
OWNER="***REMOVED***"
SDASIZE="100"
IAM="activemapper_planet_readwriteS3"
aws ec2 run-instances --image-id $AMID --count 1 --instance-type $ITYPE --iam-instance-profile Name=$IAM --key-name $KEYNAME --security-groups $SECURITY --monitoring Enabled=true --block-device-mappings "[ { \"DeviceName\": \"/dev/sda1\", \"Ebs\": { \"VolumeSize\": $SDASIZE } } ]" --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value='$INAME'}]' 'ResourceType=volume,Tags=[{Key=Owner,Value='$OWNER'}]' --instance-market-options 'MarketType=spot, SpotOptions={MaxPrice=0.08,SpotInstanceType=persistent,ValidUntil=2019-12-16T23:00:00, InstanceInterruptionBehavior=stop}'
