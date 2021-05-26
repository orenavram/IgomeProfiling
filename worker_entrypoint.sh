#!/bin/bash
source .venv/bin/activate
celery -A worker worker -P eventlet --prefetch-multiplier=1 --loglevel=info

