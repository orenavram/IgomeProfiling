import boto3 
import datetime
import re

def stop(instance, logger):
    if instance.state['Name'] != 'stopped':
            instance.stop()
            logger.info("Stop the instance - Id: {0}, Type: {1}, Public IPv4: {2}, State: {3}".format(instance.id, instance.instance_type, instance.public_ip_address, instance.state))

def is_match_name_instance(names, instance):
    for name in names:
        if re.search(name, instance.tags[0]['Value']):
            return True
    return False

def stop_machines(type_machine, name_machines, logger):
    logger.info('Stop machines...')
    try:
        ec2 = boto3.resource('ec2', region_name='us-west-2')
    except botocore.exceptions.ClientError as error:
        logger.error(f'{datetime.datetime.now()}: Error {error}  - Can\'t connect to ec2 with boto3 for stop the machines')

    list_instances = ec2.instances.all()

    types = type_machine.split(',') if type_machine else '' 
    names = name_machines.split(',') if name_machines else ''

    for instance in list_instances:
        if not instance.state['Name'] == 'stopped' and (((types and names) and (instance.instance_type in types or is_match_name_instance(names, instance))) or ((not types or instance.instance_type in types) and (not names or is_match_name_instance(names, instance)))):
            stop(instance, logger)
        else:
            logger.info(f'{datetime.datetime.now()}: Skipping on instance -  {instance.id} and not stop it.')
