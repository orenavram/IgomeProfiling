import subprocess
import boto3 
import os
import sys
import logging
import argparse
import datetime

def stop(instance, logger):
    if instance.state['Name'] != 'stopped':
            instance.stop()
            logger.info("Stop the instance - Id: {0}, Type: {1}, Public IPv4: {2}, State: {3}".format(instance.id, instance.instance_type, instance.public_ip_address, instance.state))

def stop_machines(type_machine, name_machines, logger):
    logger.info('Stop machines...')
    try:
        ec2 = boto3.resource('ec2', region_name='us-west-2')
    except:
        logger.info(f'{datetime.datetime.now()}: Can\'t connect to ec2 with boto3 for stop the machines')
        return

    list_instances = ec2.instances.all()

    types =type_machine.split('_') if type_machine else '' 
    names =name_machines.split('_') if name_machines else ''
    
    for instance in list_instances:
        if not instance.state['Name'] == 'stopped' and (types[0]=='all' or instance.instance_type in types) and (names[0]=='all' or instance.tags[0]['Value'] in names):
            stop(instance, logger)
        else:
            logger.info(f'{datetime.datetime.now()}: Skipping on instance -  {instance.id} and not stop it.')
