import subprocess
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

def stop_machines(done_path, type_machine, verbose, argv):
    print('Stop machines...')
    ec2 = boto3.resource('ec2', region_name='us-west-2')
    list_instances = ec2.instances.all()  
    types=type_machine.split('_')
    for instance in list_instances:
        if not instance.state['Name'] == 'stopped' and  instance.instance_type in types:
            stop(instance, verbose)

    with open(done_path, 'w') as f:
            f.write(' '.join(argv) + '\n')    

if __name__ == "__main__":    
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')
    parser = argparse.ArgumentParser()
    parser.add_argument('done_path', type=str, help='A path to file that signals that the script finished running successfully.')
    parser.add_argument('--type_machines', default='t2.medium_t2.2xlarge_m5a.24xlarge', help='stop only parts of types machines')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)
    logger = logging.getLogger('main')

    stop_machines(args.done_path, args.type_machines, args.verbose, sys.argv)
