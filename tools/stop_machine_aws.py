import subprocess
import paramiko
import boto3 
import os
import sys
import logging
import argparse

def stop(instance, verbose):
    if instance.state['Name'] != 'stopped':
            if verbose:
                logger.info("Id: {0}, Platform: {1},Type: {2}, Public IPv4: {3}, AMI: {4}, State: {5}\n".format(instance.id, instance.platform, instance.instance_type, instance.public_ip_address, instance.image.id, instance.state))
            instance.stop()
            if verbose:
                logger.info('After stop function:\n')
                logger.info("Id: {0}, Platform: {1}, Type: {2}, Public IPv4: {3}, AMI: {4}, State: {5}\n".format(instance.id, instance.platform, instance.instance_type, instance.public_ip_address, instance.image.id, instance.state))

def stop_machines(type_machine, verbose):
    print('Stop machines...')
    ec2 = boto3.resource('ec2')
    if type_machine:
        list_instances = [instance for instance in ec2.instances.all() if instance.instance_type==type_machine ]
    else:
        list_instances = ec2.instances.all()
    for instance in list_instances:
        if instance.state['Name'] != 'stopped':
            stop(instance, verbose)

                

if __name__ == "__main__":    
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')
    parser = argparse.ArgumentParser()
    parser.add_argument('--type_machine', help='stop only one type of machine')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)
    logger = logging.getLogger('main')

    stop_machines(args.type_machine, args.verbose)
