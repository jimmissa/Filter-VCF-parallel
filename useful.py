import sys
import re
import subprocess

url=sys.argv[1]

def get_ssh(x):
    first="ssh -L 8888:"
    middle=re.search(r'http://(.*)/.*', x).group(1)
    third=" jimmissa@narval.computecanada.ca"
    return first + middle + third

def get_url(x):
    token=re.search(r'.*token=(.*)', x).group(1)
    print("http://localhost:8888/?token=" + token)

get_url(url)

script_name = 'run_jup.sh'

with open(script_name,'w') as file:
    file.write('#!/bin/bash')
    file.write('\n')
    file.write('set -eou pipefail')
    file.write('\n')
    file.write(get_ssh(url))

subprocess.run(['bash', script_name])




