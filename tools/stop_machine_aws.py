import subprocess
import paramiko
import boto3 
import os
import sys

def stop(instance,f_log):
    if instance.state['Name']!='stopped':
            f_log.write("Id: {0}, Platform: {1},Type: {2}, Public IPv4: {3}, AMI: {4}, State: {5}\n".format(instance.id, instance.platform, instance.instance_type, instance.public_ip_address, instance.image.id, instance.state))
            instance.stop()
            f_log.write('After stop function:\n')
            f_log.write("Id: {0}, Platform: {1}, Type: {2}, Public IPv4: {3}, AMI: {4}, State: {5}\n".format(instance.id, instance.platform, instance.instance_type, instance.public_ip_address, instance.image.id, instance.state))

def stop_machines(log_path,type_machine):
    print('Stop all machines...')
    f_log=open(log_path,'w')
    ec2 = boto3.resource('ec2')
    for instance in ec2.instances.all():
        #print("Id: {0}\nPlatform: {1}\nType: {2}\nPublic IPv4: {3}\nAMI: {4}\nState: {5}\n".format(instance.id, instance.platform, instance.instance_type, instance.public_ip_address, instance.image.id, instance.state))    
        if type_machine: 
            if instance.state['Name']!='stopped' and instance.instance_type==type_machine:
                stop(instance,f_log)
        else:
            if instance.state['Name']!='stopped':
                stop(instance,f_log)

                
    f_log.close()

if __name__ == "__main__":    
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('log_path', help='A path to a txt file for log of the run')
    parser.add_argument('--type_machine', help='stop only one type of machine')
    args = parser.parse_args()
    stop_machines(args.log_path,args.type_machine)
