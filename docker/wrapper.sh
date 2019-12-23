#!/usr/bin/env bash

# gather time and memory statistics, save everything to log
echo "/usr/bin/time -f '\\\\nDEBUG_MAX_MEM:%M\\\\nDEBUG_RUNTIME:%E\\\\n' /opt/MarginPolish/build/marginPolish $@\n" > /data/marginPolish.log
eval "/usr/bin/time -f '\\nDEBUG_MAX_MEM:%M\\nDEBUG_RUNTIME:%E\\n' /opt/MarginPolish/build/marginPolish $@" 2>&1 | tee -a /data/marginPolish.log
