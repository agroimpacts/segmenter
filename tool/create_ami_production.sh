IID=i-0560f279478c59af9
instance_name="segmenter_V2"
description="UsedforParallelizationSegmentation"
aws ec2 create-image --instance-id $IID --name $instance_name --description $description
