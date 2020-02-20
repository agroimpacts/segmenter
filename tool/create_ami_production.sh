IID=***REMOVED***
instance_name="segmenter_V2"
description="UsedforParallelizationSegmentation"
aws ec2 create-image --instance-id $IID --name $instance_name --description $description
