#!/bin/bash
source .venv/bin/activate
celery -A worker worker -P eventlet --loglevel=info

