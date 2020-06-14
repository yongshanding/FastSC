#!/bin/bash

test_array=( "bv" "qgan" "qaoa" "ising" "xeb_iswap_barrier_p5" "xeb_iswap_barrier_p10" "xeb_iswap_barrier_p15" )
#test_array=( "xeb_iswap_barrier_p5" )

# Clean up
for arg in "$@"
do
    if [ "$arg" == "--clean" ] || [ "$arg" == "-c" ]
    then
        echo "Clean argument detected..."
        for i in $(seq 0 $(((${#test_array[*]}-1))))
        do
            rm  -rf ${test_array[$i]}_res/
        done
    fi
done

# Set up directories
for i in $(seq 0 $(((${#test_array[*]}-1))))
do
    echo "Checking output directories..."
    if [ ! -d "${test_array[$i]}_res/" ]; then
        mkdir ${test_array[$i]}_res/
    fi
done


echo "Begin simulation..."
# Simulate each circuit
p=0
for i in $(seq 0 $(((${#test_array[*]}-1))))
do
    #if [ "${test_array[$i]}" == "parallel_cnot" ] ; then

    for j in {2..5}; do
        qq=$(($j*$j))
        ./run.sh ${test_array[$i]} ${qq} &
        pids[${p}]=$!
        p=$((p+1))

    done

done

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done

echo "End simulation."
