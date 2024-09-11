#!/bin/bash
for i in {1..200} ; do cp grep_PG.py $i/ ; cd $i ; PG_op="$( python3 grep_PG.py fp_opt.vasp | grep symmetry | awk -F ',' '{print $2}' | awk '{print $1}' )" ; echo $i $PG_op ; cd .. ; done > "$(pwd | awk -F '/' '{print $NF }')"_PG_op_list
