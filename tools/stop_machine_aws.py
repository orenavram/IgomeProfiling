def connect_ssh(host_ssh,username_ssh,password_ssh,private_key):
    # Connect
    print('connecting ssh...')
    key = paramiko.RSAKey.from_private_key_file(private_key)
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host_ssh, username=username_ssh, password=password_ssh,pkey=key)
    print('connected')

    return client

def stop_machine(ssh_client):
    print('stop all machine...')
    


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('host', type=str, help='The host name for connect to ssh')
    parser.add_argument('user_name', type=str, help='The user name for connect ssh')
    parser.add_argument('password', type=str, help='The password foe connect ssh')
    parser.add_argument('key', type=str, help='A path for pem file of the private key')
    args = parser.parse_args()

    client=connect_ssh(args.host,args.user_name,args.password,args.pkey)
    stop_machine(client)