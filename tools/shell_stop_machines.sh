#!/bin/bash


#mini-worker
aws ec2 stop-instances --instance-ids i-03df2ae895498b6da
#large-worker-2
aws ec2 stop-instances --instance-ids i-05baefdd869a266c5
#large-worker-1
aws ec2 stop-instances --instance-ids i-0155f06ee284b92c7
#large-worker-3
aws ec2 stop-instances --instance-ids i-0e97157f6412b89d8
