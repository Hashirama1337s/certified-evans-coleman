#!/bin/bash
cd <LAB>
MC="py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/count_b4v3.py"
timeout 1200 $MC smoke 1/3 C1 > runs7/ctl_C1.log 2>&1; echo "rc=$?" >> runs7/ctl_C1.log
grep -q "C1a.*PASS" runs7/ctl_C1.log && grep -q "C1b.*PASS" runs7/ctl_C1.log || { echo "STOP: C1 control failed" > runs7/STATUS; exit 1; }
timeout 1800 $MC count 3 runs7/ctl_R3.jsonl > runs7/ctl_R3.log 2>&1; echo "rc=$?" >> runs7/ctl_R3.log
grep -q "N(R1=3) = 1 " runs7/ctl_R3.log || { echo "STOP: R1=3 control failed" > runs7/STATUS; exit 1; }
echo "controls PASS; count running" > runs7/STATUS
timeout 12240 $MC count 40 runs7/count_R40.jsonl > runs7/count_R40.log 2>&1; echo "rc=$?" >> runs7/count_R40.log
echo "count done; kstar running" > runs7/STATUS
timeout 3600 $MC kstar > runs7/kstar.log 2>&1; echo "rc=$?" >> runs7/kstar.log
echo "ALL DONE" > runs7/STATUS
