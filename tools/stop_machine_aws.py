import subprocess
import paramiko

def stop_machines(shell_script):
    print('Stop all machines...')
    sec=subprocess.call(['sh', shell_script])
    if sec==0:
        print('the machines are close')
    else:
        print('can\'t close the machines')    

def stop_manager_machines(private_key):
    print('Connecting to manager machine...')
     # Connect
    key = paramiko.RSAKey.from_private_key_file(private_key)
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect('35.160.134.140', username='ec2-user', password='',pkey=key)
    command='sudo shutdown -h now'
    stdin, stdout, stderr = ssh.exec_command(command)
    lines = stdout.readlines()
    print(lines)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('shell_script', type=str, help='A path for the shell script that close all the machine')
    #parser.add_argument('private_key', type=str, help='A path for the private key pem file')
    args = parser.parse_args()

    stop_machines(args.shell_script)
    #stop_manager_machines(args.private_key)