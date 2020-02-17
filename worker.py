# export CELERY_BROKER_URL="pyamqp://guest@localhost//"
# export CELERY_TASK_ACKS_LATE=True
# export CELERY_WORKER_PREFETCH_MULTIPLIER=1
# celery -A worker worker -P eventlet --loglevel=info # -c 1000 -n worker1@%h
from celery import Celery
from time import sleep
from eventlet.green.subprocess import Popen

app = Celery('tasks')

@app.task
def submit(args, shell = False):
    print(f'got task: args={args}, shell={shell}')
    proc = Popen(args, shell=shell)
    result = proc.wait()
    print(f'end task: args={args}, shell={shell}, result=${result}')
    return result
