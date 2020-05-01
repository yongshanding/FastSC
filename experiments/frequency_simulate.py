
import os, sys, getopt, math, re, random
import numpy as np

import time
from datetime import datetime
#import config


print("Importing qiskit...")
import qiskit
print("Importing fqc...")
import fastsc
from fastsc.coloring import success_rate_rand_coloring, success_rate_full_coloring, success_rate_layer_coloring, success_rate_google_like

from fastsc.coloring import compute_decoherence, compute_crosstalk_by_layer
from fastsc.models import Device
from fastsc.benchmarks import get_circuit

#data_path = config.DATA_PATH
#file_name = datetime.today().strftime("%h%d")

###########################################################
# Simulation Setup
###########################################################

### Device Parameters ###
#side_length = 4
omega_max = 5.0 #GHz
delta_int = 1.0
delta_ext= 0.5
delta_park = 1.0
alpha = -0.2
ejs = 8
ejl = 20
ec = 0.3
cqq = 0.012

###########################################################
# Simulation
###########################################################

#def reschedule(circ, scheduler):
#    # scheduler = qiskit, local, global
#    res = circ
#    if (scheduler == 'local'):
#        print("Local crosstalk-aware scheduler: Not yet implemented.")
#        # Local crosstalk-aware scheduler
#
#    elif (scheduler == 'global'):
#        print("Global crosstalk-aware scheduler: Not yet implemented.")
#        # Global crosstalk-aware scheduler
#
#    elif (scheduler != 'qiskit'):
#        print("Scheduler %s not recognized." % scheduler)
#    return res


def simulate(device, circuit, mapper, scheduler, freq, dist, decomp, outputfile=None, depth=0, lim_colors=0,verbose=0):
    circ = get_circuit(device.side_length * device.side_length, circuit, dep=depth)
    # Crosstalk-aware mapping yet to be implemented.
    #scheduled = reschedule(circ, scheduler)

    #if (freq == 'random'):
        # random coloring
    #    sr, avg, worst, d_before, d_after, t, c, t_act, t_2q = success_rate_rand_coloring(device, circ, scheduler, dist, decomp)

    if (freq == 'full'):
        # Full coloring
        static_coloring()
        compute_crosstalk_by_layer()
        #sr, avg, worst, d_before, d_after, t, c, t_act, t_2q = success_rate_full_coloring(device, circ, scheduler, dist, decomp, outputfile, verbose)
    elif (freq == 'layer'):
        # Layered coloring
        sr, avg, worst, d_before, d_after, t, c, t_act, t_2q = success_rate_layer_coloring(device, circ, scheduler, dist, decomp, outputfile, lim_colors, verbose)
    elif (freq == 'google'):
        # with (Google-like) tunable coupling
        sr, avg, worst, d_before, d_after, t, c, t_act, t_2q = success_rate_google_like(device, circ, scheduler, dist, decomp, outputfile, verbose)
    else:
        sr = 0.0
        avg = 0.0
        worst = 0.0
        d_before = 0
        d_after = 0
        t = 0.0
        c = 0
        t_act = 0.0
        t_2q = 0.0
    return sr, avg, worst, d_before, d_after, t, c, t_act, t_2q





###########################################################
# Main procedure
###########################################################

def main():
    random.seed(60615)
    np.random.seed(60615)
    circuit = None
    qubits = 0
    mapper = None
    scheduler = None
    freq = None
    dist = None
    decomp = None
    outputfile = None
    depth = 0
    lim_colors = 0 # when lim_colors=0 we don't limit the number of colors
    verbose = 0 # 0 - verbose, 1 - less verbose
    try:
        opt, args = getopt.getopt(sys.argv[1:], "hi:p:m:s:f:x:d:o:q:c:v:", ["help", "input=", "depth=", "mapper=", "scheduler=", "frequency=", "crosstalk=", "decomposition=", "output=", "qubits=","colors=","verbose="])
    except getopt.GetOptError as err:
        print(err)
        print("Usage: frequency_simulate.py -i <input circuit=bv> -q <num qubits (square)> -p <depth of supremacy circuit> -m <mapper=qiskit> -s <scheduler=qiskit> -f <frequency assignment=full> -x <crosstalk distance=1> -d <circuit decomposition=iswap> -o <output file> -c <max colors> -v <verbosity>")
        sys.exit(2)
    usage = True
    for o,a in opt:
        usage = False
        if o in ("-h", "--help"):
            print("Usage: frequency_simulate.py -i <input circuit=bv> -q <num qubits (square)> -p <depth of supremacy circuit> -m <mapper=qiskit> -s <scheduler=qiskit> -f <frequency assignment=full> -x <crosstalk distance=1> -d <circuit decomposition=iswap> -o <output file> -c <max colors> -v <verbosity>")
            sys.exit()
        elif o in ("-i", "--input"): # bv, qft,
            circuit = a
        elif o in ("-q", "--qubits"): # 4, 16
            qubits = int(a)
        elif o in ("-p", "--depth"):
            depth = int(a)
        elif o in ("-m", "--mapper"): # qiskit
            mapper = a
        elif o in ("-s", "--scheduler"): # qiskit, greedy, tiling
            scheduler = a
        elif o in ("-f", "--frequency"): # (qoc,) random, full, layer, google
            freq = a
        elif o in ("-x", "--crosstalk"): # 1, 2, 3
            dist = int(a)
        elif o in ("-d", "--decomposition"): # iswap, cphase
            decomp = a
        elif o in ("-o", "--output"):
            outputfile = a
        elif o in ("-c", "--colors"):
            lim_colors = int(a)
        elif o in ("-v", "--verbose"):
            verbose = int(a)
        else:
            print("Usage: frequency_simulate.py -i <input circuit=bv> -q <num qubits (square)> -p <depth of supremacy circuit> -m <mapper=qiskit> -s <scheduler=qiskit> -f <frequency assignment=full> -x <crosstalk distance=1> -d <circuit decomposition=iswap> -o <output file> -c <max colors>")
            sys.exit(2)

    if (usage):
        print("------")
        print("Usage: frequency_simulate.py -i <input circuit=bv> -q <num qubits (square)> -dep <depth of circuit> -m <mapper=qiskit> -s <scheduler=qiskit> -f <frequency assignment=full> -x <crosstalk distance=1> -d <circuit decomposition=iswap> -o <output file> -c <max colors>")
        print("------")

    if (mapper == None): mapper = 'qiskit'
    if (scheduler == None): scheduler = 'qiskit'
    if (freq == None): freq = 'full'
    if (dist == None): dist = 1
    if (decomp == None): decomp = 'iswap'
    #if (outputfile == None): outputfile = file_name

    side_length = int(np.sqrt(qubits))

    device = Device(side_length, omega_max, delta_int, delta_ext, delta_park, cqq, alpha, ejs, ejl, ec)
    start = time.time()
    # success, avg, worst, d_before, d_after, t, c, t_act, t_2q = simulate(device, circuit, mapper, scheduler, freq, dist, decomp, outputfile, depth=depth, lim_colors=lim_colors, verbose=verbose)
    success, swap_err, leak_err, d_before, d_after, t, c, t_act, t_2q = simulate(device, circuit, mapper, scheduler, freq, dist, decomp, depth=depth, lim_colors=lim_colors, verbose=verbose)
    # calculate decoherence
    decoh = compute_decoherence(device, t_act, t_2q)
    success *= decoh
    end = time.time()
    print("======================")
    print("Avg(worst) success rate per timestep: %12.10f(%12.10f)" % (avg, worst))
    print("Final success rate: %12.10f" % success)
    print("Circuit depth: %d, %d" % (d_before, d_after)) # before and after decomp 2-q gates
    print("Circuit execution time: %12.10f ns" % t)
    print("Decoherence factor: ", decoh)
    print("Compilation time: %12.10f s" % (end-start))
    print("Colors: %d" % c)


if __name__ == "__main__":
    main()
