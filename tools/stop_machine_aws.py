import subprocess
import paramiko
import boto3 

def stop_machines():
    print('Stop all machines...')
    region='us-west-2'
    #ec2 = boto3.resource('ec2',region)
    #ids = ['i-03df2ae895498b6da']
    #ec2.instances.filter(InstanceIds = ids).stop() #for stopping an ec2 instance
    
    ec2 = boto3.resource('ec2')
    instance = ec2.Instance('i-03df2ae895498b6da')
    response = instance.stop()
    print(response)
    
    #client = boto3.client('ec2',region_name=region, endpoint_url=f'https://sts.{region}.amazonaws.com')
    #response = client.stop_instances(InstanceIds=['i-03df2ae895498b6da'])
    #Hibernate=True|False,
    #DryRun=True|False,
    #Force=True|False
    #instances = ec2.instances.filter(Filters=[{'Name': 'instance-state-name', 'Values': ['running']}])
    #for instance in instances:
    #    print(instance.id, instance.instance_type)

if __name__ == "__main__":    
    stop_machines()
