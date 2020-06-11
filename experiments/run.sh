#!/bin/bash

FILE=$1
if [ "$FILE" = "" ]
then
    echo "argument 1 missing: input application name"
    exit 0
fi

qq=$2
if [ "$qq" = "" ]
then
    echo "argument 2 missing: number of qubits"
    exit 0
fi

for i in {1..5}; do
    echo "=== Iter $i ===" >> ${FILE}_res/q${qq}.out
    python3 frequency_simulate.py -i ${FILE} -q $qq -m qiskit -s qiskit -f layer -x 1 -d flexible -c 0 -v 1 -n 0 > ${FILE}_res/q$qq-layer.circ
    echo "--- ${FILE} layer q${qq} ---">> ${FILE}_res/q${qq}.out
    tail -n 8 ${FILE}_res/q$qq-layer.circ >> ${FILE}_res/q${qq}.out
    python3 frequency_simulate.py -i ${FILE} -q $qq -m qiskit -s qiskit -f opt -x 1 -d flexible -c 0 -v 1 -n 0 > ${FILE}_res/q$qq-opt.circ
    echo "--- ${FILE} opt q${qq} ---">> ${FILE}_res/q${qq}.out
    tail -n 8 ${FILE}_res/q$qq-opt.circ >> ${FILE}_res/q${qq}.out
    python3 frequency_simulate.py -i ${FILE} -q $qq -m qiskit -s qiskit -f full -x 1 -d flexible -c 0 -v 1 -n 0 > ${FILE}_res/q$qq-full.circ
    echo "--- ${FILE} full q${qq} ---">> ${FILE}_res/q${qq}.out
    tail -n 8 ${FILE}_res/q$qq-full.circ >> ${FILE}_res/q${qq}.out
    python3 frequency_simulate.py -i ${FILE} -q $qq -m qiskit -s qiskit -f google -x 1 -d flexible -c 0 -v 1 -n 0 -r 0.0 > ${FILE}_res/q$qq-google.circ
    echo "--- ${FILE} google q${qq} ---">> ${FILE}_res/q${qq}.out
    tail -n 8 ${FILE}_res/q$qq-google.circ >> ${FILE}_res/q${qq}.out
done
