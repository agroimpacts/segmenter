IID=<instance id>
instance_name="instance name here"
description="UsedforParallelizationSegmentation"
aws ec2 create-image --instance-id $IID --name $instance_name --description $description
