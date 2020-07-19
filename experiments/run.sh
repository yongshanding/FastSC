#!/bin/bash

FILE=$1
if [ "$FILE" == "" ]
then
    echo "argument 1 missing: input application name"
    exit 0
fi

qq=$2
if [ "$qq" == "" ]
then
    echo "argument 2 missing: number of qubits"
    exit 0
fi

APP=$FILE
pp=0
if [ "$FILE" == "xeb_iswap_barrier_p5" ]
then
    pp=5
    APP="xeb_iswap_barrier"
fi
if [ "$FILE" == "xeb_iswap_barrier_p10" ]
then
    pp=10
    APP="xeb_iswap_barrier"
fi
if [ "$FILE" == "xeb_iswap_barrier_p15" ]
then
    pp=15
    APP="xeb_iswap_barrier"
fi

for i in {1..2}; do
    echo "=== Iter $i ===" >> ${FILE}_res/q${qq}.out
    python3 frequency_simulate.py -i ${APP} -q $qq -p $pp -m qiskit -s qiskit -f layer -x 1 -d flexible -c 0 -v 1 -n 0 -t grid> ${FILE}_res/q$qq-layer.circ 2> ${FILE}_res/q${qq}.err
    echo "--- ${FILE} layer q${qq} ---">> ${FILE}_res/q${qq}.out
    tail -n 8 ${FILE}_res/q$qq-layer.circ >> ${FILE}_res/q${qq}.out

    python3 frequency_simulate.py -i ${APP} -q $qq -p $pp -m qiskit -s qiskit -f opt -x 1 -d flexible -c 0 -v 1 -n 0 -t grid > ${FILE}_res/q$qq-opt.circ 2> ${FILE}_res/q${qq}.err
    echo "--- ${FILE} opt q${qq} ---">> ${FILE}_res/q${qq}.out
    tail -n 8 ${FILE}_res/q$qq-opt.circ >> ${FILE}_res/q${qq}.out

    python3 frequency_simulate.py -i ${APP} -q $qq -p $pp -m qiskit -s qiskit -f full -x 1 -d flexible -c 0 -v 1 -n 0 -t grid > ${FILE}_res/q$qq-full.circ 2> ${FILE}_res/q${qq}.err
    echo "--- ${FILE} full q${qq} ---">> ${FILE}_res/q${qq}.out
    tail -n 8 ${FILE}_res/q$qq-full.circ >> ${FILE}_res/q${qq}.out

    python3 frequency_simulate.py -i ${APP} -q $qq -p $pp -m qiskit -s tiling -f google -x 1 -d flexible -c 0 -v 1 -n 0 -r 0.0 -t grid > ${FILE}_res/q$qq-google.circ 2> ${FILE}_res/q${qq}.err
    echo "--- ${FILE} google q${qq} ---">> ${FILE}_res/q${qq}.out
    tail -n 8 ${FILE}_res/q$qq-google.circ >> ${FILE}_res/q${qq}.out

    python3 frequency_simulate.py -i ${APP} -q $qq -p $pp -m qiskit -s greedy -f layer -x 1 -d flexible -c 0 -v 1 -n 0 -r 0.0 -t grid > ${FILE}_res/q$qq-uniform.circ 2> ${FILE}_res/q${qq}.err
    echo "--- ${FILE} uniform q${qq} ---">> ${FILE}_res/q${qq}.out
    tail -n 8 ${FILE}_res/q$qq-uniform.circ >> ${FILE}_res/q${qq}.out

    python3 frequency_simulate.py -i ${APP} -q $qq -p $pp -m qiskit -s qiskit -f full -x 1 -d iswap -c 0 -v 1 -n 0 -r 0.0 -u 1 -t grid > ${FILE}_res/q$qq-naive.circ 2> ${FILE}_res/q${qq}.err
    echo "--- ${FILE} naive q${qq} ---">> ${FILE}_res/q${qq}.out
    tail -n 8 ${FILE}_res/q$qq-naive.circ >> ${FILE}_res/q${qq}.out
done
