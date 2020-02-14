#/bin/bash
$(aws ecr get-login --no-include-email --region us-west-2)
docker tag webiks/igome-profile:latest 223455578796.dkr.ecr.us-west-2.amazonaws.com/igome-profile:latest
docker push 223455578796.dkr.ecr.us-west-2.amazonaws.com/igome-profile:latest

