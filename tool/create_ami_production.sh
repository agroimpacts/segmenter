IID=i-0051f2d8ab522963a
instance_name="segmenter_ami"
description="UsedforParallelizationSegmentation"
aws ec2 create-image --instance-id $IID --name $instance_name --description $description
