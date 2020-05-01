"""
coloring.py -  Functions for frequency assignment via graph coloring
"""
import numpy as np
from qiskit import Aer, BasicAer, QuantumCircuit, QuantumRegister, execute
from qiskit.extensions.standard import *
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.transpiler import PassManager, CouplingMap
from qiskit.compiler import transpile
from qiskit.transpiler.passes import (BasicSwap, CXCancellation)
from qiskit.quantum_info.operators import Operator

#from fqc.util import get_unitary
from fastsc.util import get_connectivity_graph, get_aug_line_graph, get_map_circuit, get_layer_circuits, get_nearest_neighbor_coupling_list
import networkx as nx

import h5py

#GREEDY = 1 # Greedy algorithm parameter: distance to be considered as collision
#CQQ = 0.012*2*np.pi
ALPHA = -0.2 # aharmonicity - change this.
EJS = 8
EJL = 20
EC = 0.3
CQQ = 0.012
sigma = 0.0 # standard deviation of the Gaussian noise
# goog_uni = True # whether google_like uses uniform interaction frequency
full_uni = True # whether full_coloring uses uniform interaction frequency
res_coupling = 0. # strength of residual coupling (for google_like), 0~1
GATETIMES = {'unitary': 55,'rz': 30, 'z':30, 'u1': 30, 's': 30, 't': 30, 'rx': 30, 'x': 30, 'u2': 30, 'ry': 30, 'y': 30, 'u3': 30, 'h': 30, 'measure': 0.0, 'barrier': 0.0} # all in ns

def compute_crosstalk(iswaps, coupling, qubit_freqs):
    # returns error rate based on a qubit-frequency configuration
    # ISWAP interaction
    residual_success = 1.0
    Cqq = CQQ
    for (q1, q2) in coupling:
        if (q1 > 1 or q2 > 1):
            continue
        print((q1,q2))
        J = 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq
        tau = np.pi/(2*J)
        epsilon_leak = (16*J**2)/(ALPHA**2+16*J**2) * np.sin(0.5*np.sqrt(ALPHA**2+16*J**2)*tau)**2
        print('leakage: '+ str(epsilon_leak))
        if (q1,q2) in iswaps or (q2,q1) in iswaps:
            print('q' + str(q1) + ': ' + str(qubit_freqs[q1]) + ', q' + str(q2) + ': ' + str(qubit_freqs[q2]))
            assert(qubit_freqs[q1] == qubit_freqs[q2])
            # assume stay on resonant for pi/2g time
            # Compute leakage
            residual_success *= 1- epsilon_leak
        else:
            assert(qubit_freqs[q1] != qubit_freqs[q2])
            #print("q%d q%d: %f,%f" % (q1, q2, qubit_freqs[q1], qubit_freqs[q2]))
            delta_omega = np.absolute(qubit_freqs[q1] - qubit_freqs[q2])
            residual_coupling = J**2 / (delta_omega * 4 * np.pi)
            #print("here")
            #print(J)
            #print(residual_coupling)
            residual_success *= 1- epsilon_leak
            tau = np.pi/(2 * 0.5 * np.sqrt(max(qubit_freqs[q1], qubit_freqs[q2])**2) * Cqq)
            residual_success *= 1 - np.sin(residual_coupling * tau)**2
            #print(residual_success)
    return 1 - residual_success

def swap_channel_google(coupling, active, qubit_freqs, taus):
    # print("inside swap_channel_google")
    # print("active tiling:", active)
    Cqq = CQQ
    res = []
    for (i, (q1,q2)) in enumerate(coupling):
        J = 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq
        delta_omega = np.absolute(qubit_freqs[q1] - qubit_freqs[q2])
        if delta_omega == 0:
            residual_coupling = J
        else:
            residual_coupling = J**2 / (delta_omega * 4 * np.pi)
        tau = taus[i]
        if (q1,q2) in active or (q2,q1) in active:
            res.append(np.sin(residual_coupling * tau)**2)
        else:
            res.append(res_coupling*np.sin(residual_coupling * tau)**2)
    # print("couplings:",coupling)
    # print("swap chances:",res)
    return res

def swap_channel(coupling, qubit_freqs, taus):
    # return prob of rabi oscillation between 01 and 10 (i.e. iswap)
    # for each pair of coupled qubits
    # taus is list of hold durations: max swap happens at 2pi/g
    Cqq = CQQ
    res = []
    for (i, (q1,q2)) in enumerate(coupling):
        J = 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq
        delta_omega = np.absolute(qubit_freqs[q1] - qubit_freqs[q2])
        if delta_omega == 0:
            residual_coupling = J
        else:
            residual_coupling = J**2 / (delta_omega * 4 * np.pi)
        tau = taus[i]
        res.append(np.sin(residual_coupling * tau)**2)
    return res

def leak_channel_google(coupling, active, qubit_freqs, alphas, taus):
    # print("inside leak_channel_google")
    # print("active tiling:", active)
    Cqq = CQQ
    res = []
    for (i, (q1,q2)) in enumerate(coupling):
        # 02
        f1 = qubit_freqs[q1]
        f2 = qubit_freqs[q2] + alphas[q2]
        J = 0.5 * np.sqrt(f1 * f2) * Cqq
        delta_omega = np.absolute(f1 - f2)
        if delta_omega == 0:
            residual_coupling = np.sqrt(2) * J
        else:
            residual_coupling = (np.sqrt(2)*J)**2 / (delta_omega * 4 * np.pi)
        tau = taus[i]
        epsilon_leak02 = np.sin(residual_coupling * tau)**2
        ## 20
        #f1 = qubit_freqs[q1] + alphas[q2]
        #f2 = qubit_freqs[q2]
        #J = 0.5 * np.sqrt(f1 * f2) * Cqq
        #delta_omega = np.absolute(f1 - f2)
        #if delta_omega == 0:
        #    residual_coupling = np.sqrt(2) * J
        #else:
        #    residual_coupling = (np.sqrt(2)*J)**2 / (delta_omega * 4 * np.pi)
        #tau = taus[i]
        #epsilon_leak20 = np.sin(residual_coupling * tau)**2
        #res.append((epsilon_leak02 + epsilon_leak20) / 2)
        if (q1,q2) in active or (q2,q1) in active:
            res.append(epsilon_leak02)
        else:
            res.append(epsilon_leak02*res_coupling)
    #print("leak chances:",res)
    return res

def leak_channel(coupling, qubit_freqs, alphas, taus):
    # return prob of rabi oscillation between 11 and (20+02)/sqrt(2) (i.e. leakage)
    # for each pair of coupled qubits
    # taus is list of hold durations
    Cqq = CQQ
    res = []
    for (i, (q1,q2)) in enumerate(coupling):
        # 02
        f1 = qubit_freqs[q1]
        f2 = qubit_freqs[q2] + alphas[q2]
        J = 0.5 * np.sqrt(f1 * f2) * Cqq
        delta_omega = np.absolute(f1 - f2)
        if delta_omega == 0:
            residual_coupling = np.sqrt(2) * J
        else:
            residual_coupling = (np.sqrt(2)*J)**2 / (delta_omega * 4 * np.pi)
        tau = taus[i]
        epsilon_leak02 = np.sin(residual_coupling * tau)**2
        ## 20
        #f1 = qubit_freqs[q1] + alphas[q2]
        #f2 = qubit_freqs[q2]
        #J = 0.5 * np.sqrt(f1 * f2) * Cqq
        #delta_omega = np.absolute(f1 - f2)
        #if delta_omega == 0:
        #    residual_coupling = np.sqrt(2) * J
        #else:
        #    residual_coupling = (np.sqrt(2)*J)**2 / (delta_omega * 4 * np.pi)
        #tau = taus[i]
        #epsilon_leak20 = np.sin(residual_coupling * tau)**2
        #res.append((epsilon_leak02 + epsilon_leak20) / 2)
        res.append(epsilon_leak02)
    return res

def compute_crosstalk_iswaps(iswaps, coupling, qubit_freqs, taus, gt):
    # returns error rate of simultaneous iswaps
    success = 1.0
    swap_success = 1.0
    leak_success = 1.0
    Cqq = CQQ
    # one complete swap: tau = pi/2g
    all_taus = []
    for (q1,q2) in coupling:
        if (q1,q2) in taus:
            all_taus.append(taus[(q1,q2)])
        else:
            all_taus.append(gt)
        #taus.append(c * np.pi/(2 * 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq))
    #taus = [np.pi/(2 * 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2])) * Cqq for (q1,q2) in coupling]
    #taus = [tau * tau_special[(q1,q2)] if (q1,q2) in tau_special or (q2,q1) in tau_special else tau for (i,tau) in taus]
    prob_swap = swap_channel(coupling, qubit_freqs, all_taus)
    alphas = [ALPHA for f in qubit_freqs] #TODO
    prob_leak = leak_channel(coupling, qubit_freqs, alphas, all_taus)
    #print("swap: ", prob_swap)
    #print("leak: ", prob_leak)
    for (i, (q1,q2)) in enumerate(coupling):
        if not((q1,q2) in iswaps or (q2,q1) in iswaps):
            success *= 1 - prob_swap[i]
            swap_success *= 1 - prob_swap[i]
        ## sqrt iswap
        #elif ((q1,q2) in tau_special) or ((q2,q1) in tau_special):
        #    print("sqrt iswap 0.5-prob_swap ",abs(0.5 - prob_swap[i]))
        #    success *= 1 - abs(0.5 - prob_swap[i])
        #    swap_success *= 1 - abs(0.5 - prob_swap[i])
        ## iswap
        #else:
        #    #print("2nd loop:")
        #    print("iswap prob_swap ",prob_swap[i])
        #    #print("qubit freq 1",qubit_freqs[q1])
        #    #print("qubit freq 2",qubit_freqs[q2])
        #    #print("taus[i]:",taus[i])
        #    #print("g",np.pi/(2*0.5*np.sqrt(qubit_freqs[q1]*qubit_freqs[q2])*Cqq))
        #    #print("j", 0.5*np.sqrt(qubit_freqs[q1]*qubit_freqs[q2])*Cqq)
        #    success *= prob_swap[i]
        #    swap_success *= prob_swap[i]
        #    print("iswap success ",success)
        #    print("iswap swap_success ",swap_success)

        success *= 1 - prob_leak[i]
        leak_success *= 1 - prob_leak[i]
    return 1 - success, 1-swap_success, 1-leak_success

def compute_crosstalk_cphase(iswaps, coupling, qubit_freqs, taus, gt):
    # returns error rate of simultaneous cphases
    success = 1.0
    swap_success = 1.0
    leak_success = 1.0
    Cqq = CQQ
    all_taus = []
    alphas = [ALPHA for f in qubit_freqs] #TODO
    for (q1,q2) in coupling:
        if (q1, q2) in taus:
            all_taus.append(taus[(q1,q2)])
        else:
            all_taus.append(gt)
        #    c = 1.0
        #    f1 = qubit_freqs[q1]
        #    f2 = qubit_freqs[q2]+alphas[q2]
        #elif (q1,q2) in tau_special:
        #    c = tau_special[(q1,q2)]
        #    f1 = qubit_freqs[q1]
        #    f2 = qubit_freqs[q2]
        #elif (q2,q1) in tau_special:
        #    c = tau_special[(q2,q1)]
        #    f1 = qubit_freqs[q1]
        #    f2 = qubit_freqs[q2]
        #else:
        #    c = 1.0
        #    f1 = qubit_freqs[q1]
        #    f2 = qubit_freqs[q2]
        #taus.append(c * np.pi/(np.sqrt(2) * 0.5 * np.sqrt(f1 * f2) * Cqq))
    # one complete swap and back: tau = 2 * (pi/(2*sqrt(2)*g) = pi/sqrt(2)g
    #taus = [np.pi/(np.sqrt(2) * 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2])) * Cqq for (q1,q2) in coupling]
    #taus = [tau * tau_special[(q1,q2)] if (q1,q2) in tau_special or (q2,q1) in tau_special else tau for tau in taus]
    prob_swap = swap_channel(coupling, qubit_freqs, all_taus)
    prob_leak = leak_channel(coupling, qubit_freqs, alphas, all_taus)
    #print("swap: ", prob_swap)
    #print("leak: ", prob_leak)
    for (i, (q1,q2)) in enumerate(coupling):
        if not((q1,q2) in iswaps or (q2,q1) in iswaps):
            success *= 1 - prob_swap[i]
            swap_success *= 1 - prob_swap[i]
        success *= 1 - prob_leak[i]
        leak_success *= 1 - prob_leak[i]
    return 1 - success, 1-swap_success, 1-leak_success

def compute_crosstalk_flexible(cphases, iswaps, coupling, qubit_freqs, taus, gt):
    success = 1.0
    swap_success = 1.0
    leak_success = 1.0
    Cqq = CQQ
    all_taus = []
    alphas = [ALPHA for f in qubit_freqs]
    for (q1,q2) in coupling:
        if (q1,q2) in taus:
            all_taus.append(taus[(q1,q2)])
        else:
            all_taus.append(gt)
        #if (q1,q2) in iswaps or (q2,q1) in iswaps:
        #    taus.append(c * np.pi/(2 * 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq))
        #    #print("1st loop:")
        #    #print("qubit freq 1",qubit_freqs[q1])
        #    #print("qubit freq 2",qubit_freqs[q2])
        #    #print("g",np.pi/(2*0.5*np.sqrt(qubit_freqs[q1]*qubit_freqs[q2])*Cqq))
        #    #print("c:",c)
        #    #print("tau from 1st loop:",c * np.pi/(2 * 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq))
        #else:
        #    taus.append(c * np.pi/(np.sqrt(2) * 0.5 * np.sqrt(qubit_freqs[q1] * (qubit_freqs[q2]+alphas[q2]))* Cqq))
    prob_swap = swap_channel(coupling, qubit_freqs, all_taus)
    prob_leak = leak_channel(coupling, qubit_freqs, alphas, all_taus)
    for (i, (q1,q2)) in enumerate(coupling):
        # cphase
        if not((q1,q2) in iswaps or (q2,q1) in iswaps):
            success *= 1 - prob_swap[i]
            swap_success *= 1 - prob_swap[i]
        ## sqrt iswap
        #elif ((q1,q2) in tau_special) or ((q2,q1) in tau_special):
        #    print("sqrt iswap 0.5-prob_swap ",abs(0.5 - prob_swap[i]))
        #    success *= 1 - abs(0.5 - prob_swap[i])
        #    swap_success *= 1 - abs(0.5 - prob_swap[i])
        ## iswap
        #else:
        #    #print("2nd loop:")
        #    print("iswap prob_swap ",prob_swap[i])
        #    #print("qubit freq 1",qubit_freqs[q1])
        #    #print("qubit freq 2",qubit_freqs[q2])
        #    #print("taus[i]:",taus[i])
        #    #print("g",np.pi/(2*0.5*np.sqrt(qubit_freqs[q1]*qubit_freqs[q2])*Cqq))
        #    #print("j", 0.5*np.sqrt(qubit_freqs[q1]*qubit_freqs[q2])*Cqq)
        #    success *= prob_swap[i]
        #    swap_success *= prob_swap[i]
        #    print("iswap success ",success)
        #    print("iswap swap_success ",swap_success)
        success *= 1 - prob_leak[i]
        leak_success *= 1 - prob_leak[i]
    return 1 - success, 1-swap_success, 1-leak_success


def compute_crosstalk_iswaps_gmon(iswaps, coupling, qubit_freqs, taus, gt, tilings, tiling_idx):
    # returns error rate of simultaneous iswaps
    if len(iswaps)==0:
        active = []
    else:
        active = tilings[tiling_idx]
    success = 1.0
    swap_success = 1.0
    leak_success = 1.0
    Cqq = CQQ
    # one complete swap: tau = pi/2g
    all_taus = []
    for (q1,q2) in coupling:
        if (q1,q2) in taus:
            all_taus.append(taus[(q1,q2)])
        else:
            all_taus.append(gt)
        #taus.append(c * np.pi/(2 * 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq))
    #taus = [np.pi/(2 * 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2])) * Cqq for (q1,q2) in coupling]
    #taus = [tau * tau_special[(q1,q2)] if (q1,q2) in tau_special or (q2,q1) in tau_special else tau for (i,tau) in taus]
    prob_swap = swap_channel_google(coupling, active, qubit_freqs, all_taus)
    alphas = [ALPHA for f in qubit_freqs] #TODO
    prob_leak = leak_channel_google(coupling, active, qubit_freqs, alphas, all_taus)
    #print("swap: ", prob_swap)
    #print("leak: ", prob_leak)
    for (i, (q1,q2)) in enumerate(coupling):
        if not((q1,q2) in iswaps or (q2,q1) in iswaps):
            success *= 1 - prob_swap[i]
            swap_success *= 1 - prob_swap[i]
        ## sqrt iswap
        #elif ((q1,q2) in tau_special) or ((q2,q1) in tau_special):
        #    print("sqrt iswap 0.5-prob_swap ",abs(0.5 - prob_swap[i]))
        #    success *= 1 - abs(0.5 - prob_swap[i])
        #    swap_success *= 1 - abs(0.5 - prob_swap[i])
        ## iswap
        #else:
        #    #print("2nd loop:")
        #    print("iswap prob_swap ",prob_swap[i])
        #    #print("qubit freq 1",qubit_freqs[q1])
        #    #print("qubit freq 2",qubit_freqs[q2])
        #    #print("taus[i]:",taus[i])
        #    #print("g",np.pi/(2*0.5*np.sqrt(qubit_freqs[q1]*qubit_freqs[q2])*Cqq))
        #    #print("j", 0.5*np.sqrt(qubit_freqs[q1]*qubit_freqs[q2])*Cqq)
        #    success *= prob_swap[i]
        #    swap_success *= prob_swap[i]
        #    print("iswap success ",success)
        #    print("iswap swap_success ",swap_success)

        success *= 1 - prob_leak[i]
        leak_success *= 1 - prob_leak[i]
    return 1 - success, 1-swap_success, 1-leak_success

def compute_crosstalk_cphase_gmon(iswaps, coupling, qubit_freqs, taus, gt, tilings, tiling_idx):
    # returns error rate of simultaneous cphases
    if len(cphases)==0:
        active = []
    else:
        active = tilings[tiling_idx]
    success = 1.0
    swap_success = 1.0
    leak_success = 1.0
    Cqq = CQQ
    all_taus = []
    alphas = [ALPHA for f in qubit_freqs] #TODO
    for (q1,q2) in coupling:
        if (q1, q2) in taus:
            all_taus.append(taus[(q1,q2)])
        else:
            all_taus.append(gt)
        #    c = 1.0
        #    f1 = qubit_freqs[q1]
        #    f2 = qubit_freqs[q2]+alphas[q2]
        #elif (q1,q2) in tau_special:
        #    c = tau_special[(q1,q2)]
        #    f1 = qubit_freqs[q1]
        #    f2 = qubit_freqs[q2]
        #elif (q2,q1) in tau_special:
        #    c = tau_special[(q2,q1)]
        #    f1 = qubit_freqs[q1]
        #    f2 = qubit_freqs[q2]
        #else:
        #    c = 1.0
        #    f1 = qubit_freqs[q1]
        #    f2 = qubit_freqs[q2]
        #taus.append(c * np.pi/(np.sqrt(2) * 0.5 * np.sqrt(f1 * f2) * Cqq))
    # one complete swap and back: tau = 2 * (pi/(2*sqrt(2)*g) = pi/sqrt(2)g
    #taus = [np.pi/(np.sqrt(2) * 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2])) * Cqq for (q1,q2) in coupling]
    #taus = [tau * tau_special[(q1,q2)] if (q1,q2) in tau_special or (q2,q1) in tau_special else tau for tau in taus]
    prob_swap = swap_channel_google(coupling, active, qubit_freqs, all_taus)
    prob_leak = leak_channel_google(coupling, active, qubit_freqs, alphas, all_taus)
    #print("swap: ", prob_swap)
    #print("leak: ", prob_leak)
    for (i, (q1,q2)) in enumerate(coupling):
        if not((q1,q2) in iswaps or (q2,q1) in iswaps):
            success *= 1 - prob_swap[i]
            swap_success *= 1 - prob_swap[i]
        success *= 1 - prob_leak[i]
        leak_success *= 1 - prob_leak[i]
    return 1 - success, 1-swap_success, 1-leak_success

def compute_crosstalk_flexible_gmon(cphases, iswaps, coupling, qubit_freqs, taus, gt, tilings, tiling_idx):
    if len(cphases)==0 and len(iswaps)==0:
        active = []
    else:
        active = tilings[tiling_idx]
    success = 1.0
    swap_success = 1.0
    leak_success = 1.0
    Cqq = CQQ
    all_taus = []
    alphas = [ALPHA for f in qubit_freqs]
    for (q1,q2) in coupling:
        if (q1,q2) in taus:
            all_taus.append(taus[(q1,q2)])
        else:
            all_taus.append(gt)
        #if (q1,q2) in iswaps or (q2,q1) in iswaps:
        #    taus.append(c * np.pi/(2 * 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq))
        #    #print("1st loop:")
        #    #print("qubit freq 1",qubit_freqs[q1])
        #    #print("qubit freq 2",qubit_freqs[q2])
        #    #print("g",np.pi/(2*0.5*np.sqrt(qubit_freqs[q1]*qubit_freqs[q2])*Cqq))
        #    #print("c:",c)
        #    #print("tau from 1st loop:",c * np.pi/(2 * 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq))
        #else:
        #    taus.append(c * np.pi/(np.sqrt(2) * 0.5 * np.sqrt(qubit_freqs[q1] * (qubit_freqs[q2]+alphas[q2]))* Cqq))
    prob_swap = swap_channel_google(coupling, active, qubit_freqs, all_taus)
    prob_leak = leak_channel_google(coupling, active, qubit_freqs, alphas, all_taus)
    for (i, (q1,q2)) in enumerate(coupling):
        # cphase
        if not((q1,q2) in iswaps or (q2,q1) in iswaps):
            success *= 1 - prob_swap[i]
            swap_success *= 1 - prob_swap[i]
        ## sqrt iswap
        #elif ((q1,q2) in tau_special) or ((q2,q1) in tau_special):
        #    print("sqrt iswap 0.5-prob_swap ",abs(0.5 - prob_swap[i]))
        #    success *= 1 - abs(0.5 - prob_swap[i])
        #    swap_success *= 1 - abs(0.5 - prob_swap[i])
        ## iswap
        #else:
        #    #print("2nd loop:")
        #    print("iswap prob_swap ",prob_swap[i])
        #    #print("qubit freq 1",qubit_freqs[q1])
        #    #print("qubit freq 2",qubit_freqs[q2])
        #    #print("taus[i]:",taus[i])
        #    #print("g",np.pi/(2*0.5*np.sqrt(qubit_freqs[q1]*qubit_freqs[q2])*Cqq))
        #    #print("j", 0.5*np.sqrt(qubit_freqs[q1]*qubit_freqs[q2])*Cqq)
        #    success *= prob_swap[i]
        #    swap_success *= prob_swap[i]
        #    print("iswap success ",success)
        #    print("iswap swap_success ",swap_success)
        success *= 1 - prob_leak[i]
        leak_success *= 1 - prob_leak[i]
    return 1 - success, 1-swap_success, 1-leak_success

def get_flux_noise(omega, delta):
    # Return the deviation of frequency due to gaussian noise in flux
    # omega is the target frequency, delta is the std dev.
    Ejs = EJS
    Ejl = EJL
    Ec = EC
    gamma1 = Ejl/Ejs
    d1 = (gamma1 - 1)/(gamma1 + 1)
    flux_est = np.arccos((omega-3.75)/1.25)/(2*np.pi)
    flux_new = flux_est + np.random.normal(0,delta)
    def _Em(phi, m):
        # Asymmetric transmon
        return -(Ejs + Ejl)/2 + np.sqrt(4*(Ejs + Ejl)*Ec*np.sqrt(np.cos(np.pi*phi)**2 + d1**2*np.sin(np.pi*phi)**2))*(m + 1/2) - Ec*(6*m**2 + 6*m + 3)/12
    omega_old = _Em(flux_est, 1) - _Em(flux_est, 0)
    omega_new = _Em(flux_new, 1) - _Em(flux_new, 0)
    if omega_new > 5.0:
        omega_new = 5.0
    elif omega_new < 2.5:
        omega_new = 2.5
    return omega_new - omega_old


def greedy_reschedule(coupling, edges):
    #print("=== coupling ===")
    #print(coupling)
    #print("--- edges ---")
    #print(edges)
    #print("===")
    newedges = []
    leftover = []
    # Priority based on order in the edges list
    ind = [False for _ in range(len(edges))]

    for (i, (u,v)) in enumerate(edges):
        conf = False
        for (j, (w,x)) in enumerate(edges):
            if (not i==j) and ind[j]:
                if ((u,w) in coupling) or ((u,x) in coupling) or ((v,w) in coupling) or ((v,x) in coupling):

                    ind[i] = False
                    conf = True
                    break
        if not conf:
            ind[i] = True

    for (i, e) in enumerate(edges):
        if (ind[i]):
            newedges.append(e)
        else:
            leftover.append(e)
    #if (len(leftover) > 0):
    #    print("=== newedges ===")
    #    print(newedges)
    #    print("--- leftover ---")
    #    print(leftover)
    #    print("===")
    return newedges, leftover, ind


def update_data(freqsdata, gatesdata, qubit_freqs, gates, num_q):
    freqsdata.append([qubit_freqs[i] for i in range(num_q)])
    gatesdata.append([g for g in gates])
    return

def write_circuit(freqsdata, gatesdata, outf):
    #print(freqsdata)
    freq_dset = outf.create_dataset("freqs", data=np.array(freqsdata))

    #print(gatesdata)
    dt = np.dtype([('name', np.string_, 32), ('qargs', np.int_, (2,))])
    for (t, gates) in enumerate(gatesdata):
        grp = outf.create_group("time%d" % t)
        npdata = np.array([(np.string_(name), qargs) for (name, qargs) in gates], dtype=dt)
        gate_dset = grp.create_dataset("gates", data=npdata, dtype=dt)
    print("\nCircuit data written to " + outf.filename + "\n")
    return

def iswap():
    return Operator([[1., 0., 0., 0.],[0., 0., -1.j, 0.],[0., -1.j, 0., 0.], [0., 0., 0., 1.]])

def sqrtiswap():
    return Operator([[1., 0., 0., 0.],[0., 1./np.sqrt(2), -1.j/np.sqrt(2), 0.],[0., -1.j/np.sqrt(2), 1./np.sqrt(2), 0.], [0., 0., 0., 1.]])

def decompose_layer(layer, decomp):
    new_qregs = layer.qregs
    new_cregs = layer.cregs
    if len(new_qregs) > 1:
        print("Multiple q registers.")
    if len(new_cregs) == 0:
        new_circuit = QuantumCircuit(new_qregs[0])
    elif len(new_cregs) == 1:
        new_circuit = QuantumCircuit(new_qregs[0], new_cregs[0])
    else:
        print("Multiple c registers.")
        new_circuit = QuantumCircuit(new_qregs[0])

    for inst, qargs, cargs in layer.data:
        #print(inst.name)
        if (inst.name == 'measure' or inst.name == 'barrier'):
            new_circuit.data.append((inst, qargs, cargs))
        elif (inst.name == 'cx' or inst.name == 'cnot' or inst.label == 'cnot'):
            if decomp == 'iswap':
                new_circuit.rz(-np.pi/2, qargs[0])
                new_circuit.rx(np.pi/2, qargs[1])
                new_circuit.rz(np.pi/2, qargs[1])
                new_circuit.unitary(iswap(), [qargs[0], qargs[1]], label='iswap')
                new_circuit.rx(np.pi/2, qargs[0])
                new_circuit.unitary(iswap(), [qargs[0], qargs[1]], label='iswap')
                new_circuit.rz(np.pi/2, qargs[1])
            elif decomp == 'cphase' or decomp == 'mixed':
                new_circuit.h(qargs[1])
                new_circuit.cz(qargs[0],qargs[1])
                new_circuit.h(qargs[1])
            else:
                print("Decomposition method %s not recognized. Try iswap or cphase." % decomp)
        elif (inst.name == 'swap' or inst.label == 'swap'):
            if decomp == 'iswap' or decomp == 'mixed':
                new_circuit.unitary(sqrtiswap(), [qargs[0], qargs[1]], label='sqrtiswap')
                new_circuit.rx(np.pi/2, qargs[0])
                new_circuit.rx(np.pi/2, qargs[1])
                new_circuit.unitary(sqrtiswap(), [qargs[0], qargs[1]], label='sqrtiswap')
                new_circuit.rx(-np.pi/2, qargs[1])
                new_circuit.rx(-np.pi/2, qargs[0])
                new_circuit.ry(np.pi/2, qargs[1])
                new_circuit.ry(np.pi/2, qargs[0])
                new_circuit.unitary(sqrtiswap(), [qargs[0], qargs[1]], label='sqrtiswap')
                new_circuit.ry(-np.pi/2, qargs[1])
                new_circuit.ry(-np.pi/2, qargs[0])
            elif decomp == 'cphase':
                new_circuit.h(qargs[1])
                new_circuit.cz(qargs[0],qargs[1])
                new_circuit.h(qargs[1])
                new_circuit.h(qargs[0])
                new_circuit.cz(qargs[1],qargs[0])
                new_circuit.h(qargs[0])
                new_circuit.h(qargs[1])
                new_circuit.cz(qargs[0],qargs[1])
                new_circuit.h(qargs[1])
            else:
                print("Decomposition method %s not recognized." % decomp)
        elif (inst.name == 'iswap' or inst.label == 'iswap'):
            new_circuit.unitary(iswap(), [qargs[0], qargs[1]], label='iswap')
        else:
            new_circuit.data.append((inst, qargs, cargs))
    return get_layer_circuits(new_circuit)

def decompose_layer_flexible(layer, G_crosstalk, verbose):
    new_qregs = layer.qregs
    new_cregs = layer.cregs
    if len(new_qregs) > 1:
        print("Multiple q registers.")
    if len(new_cregs) == 0:
        new_circuit = QuantumCircuit(new_qregs[0])
    elif len(new_cregs) == 1:
        new_circuit = QuantumCircuit(new_qregs[0], new_cregs[0])
    else:
        print("Multiple c registers.")
        new_circuit = QuantumCircuit(new_qregs[0])

    to_decompose = [] # all the two qubit gates (q1,q2) in this layer that are also in G_crosstalk.nodes()
    cnot_gates = []
    swap_gates = []
    for inst, qargs, cargs in layer.data:
        if (inst.name == 'measure' or inst.name == 'barrier'):
            continue
        if (inst.name == 'cx' or inst.name == 'cnot' or inst.label == 'cnot' or inst.name == 'swap' or inst.label == 'swap'):
            q1, q2 = qargs[0].index, qargs[1].index
            if (q1,q2) in G_crosstalk.nodes():
                to_decompose.append((q1,q2))
                if (inst.name == 'cx' or inst.name == 'cnot' or inst.label == 'cnot'):
                    cnot_gates.append((q1,q2))
                else:
                    swap_gates.append((q1,q2))
            elif (q2,q1) in G_crosstalk.nodes():
                to_decompose.append((q2,q1))
                if (inst.name == 'cx' or inst.name == 'cnot' or inst.label == 'cnot'):
                    cnot_gates.append((q1,q2))
                else:
                    swap_gates.append((q1,q2))
            else:
                print("Error: 2-qubit gate not a node in G_crosstalk")
    if verbose == 0:
        print("to_decompose:",to_decompose)
        print("CNOT:",cnot_gates)
        print("SWAP:",swap_gates)

    decomp_iswap = [] # the gates that will be decomposed with iswap
    decomp_cphase = [] # the gates that will be decomposed with cphase
    temp = []
    while to_decompose:
        if to_decompose[-1] in cnot_gates:
            decomp_cphase.append(to_decompose[-1])
            temp.append((to_decompose.pop(),1)) # 0-iswap, 1-cphase
        else:
            decomp_iswap.append(to_decompose[-1])
            temp.append((to_decompose.pop(),0)) # 0-iswap, 1-cphase
        while temp:
            v, d = temp.pop()
            for n in G_crosstalk[v]: # loop through v's neighbors in G_crosstalk
                # if we need to decompose v's neighbor
                if n in to_decompose:
                    to_decompose.remove(n)
                    # if v is decomposed with iswap, decompose n with cphase
                    if d == 0:
                        decomp_cphase.append(n)
                        temp.append((n,1))
                    # if v is decomposed with cphase, decompose n with iswap
                    else:
                        decomp_iswap.append(n)
                        temp.append((n,0))
    if verbose == 0:
        print("decomp_iswap:",decomp_iswap)
        print("decomp_cphase:",decomp_cphase)

    for inst, qargs, cargs in layer.data:
        if (inst.name == 'measure' or inst.name == 'barrier'):
            new_circuit.data.append((inst, qargs, cargs))
        elif (inst.name == 'cx' or inst.name == 'cnot' or inst.label == 'cnot'):
            q1, q2 = qargs[0].index, qargs[1].index
            if (q1,q2) in decomp_iswap or (q2,q1) in decomp_iswap:
                new_circuit.rz(-np.pi/2, qargs[0])
                new_circuit.rx(np.pi/2, qargs[1])
                new_circuit.rz(np.pi/2, qargs[1])
                new_circuit.unitary(iswap(), [qargs[0], qargs[1]], label='iswap')
                new_circuit.rx(np.pi/2, qargs[0])
                new_circuit.unitary(iswap(), [qargs[0], qargs[1]], label='iswap')
                new_circuit.rz(np.pi/2, qargs[1])
            elif (q1,q2) in decomp_cphase or (q2,q1) in decomp_cphase:
                new_circuit.h(qargs[1])
                new_circuit.cz(qargs[0],qargs[1])
                new_circuit.h(qargs[1])
        elif (inst.name == 'swap' or inst.label == 'swap'):
            q1, q2 = qargs[0].index, qargs[1].index
            if (q1,q2) in decomp_iswap or (q2,q1) in decomp_iswap:
                new_circuit.unitary(sqrtiswap(), [qargs[0], qargs[1]], label='sqrtiswap')
                new_circuit.rx(np.pi/2, qargs[0])
                new_circuit.rx(np.pi/2, qargs[1])
                new_circuit.unitary(sqrtiswap(), [qargs[0], qargs[1]], label='sqrtiswap')
                new_circuit.rx(-np.pi/2, qargs[1])
                new_circuit.rx(-np.pi/2, qargs[0])
                new_circuit.ry(np.pi/2, qargs[1])
                new_circuit.ry(np.pi/2, qargs[0])
                new_circuit.unitary(sqrtiswap(), [qargs[0], qargs[1]], label='sqrtiswap')
                new_circuit.ry(-np.pi/2, qargs[1])
                new_circuit.ry(-np.pi/2, qargs[0])
            elif (q1,q2) in decomp_cphase or (q2,q1) in decomp_cphase:
                new_circuit.h(qargs[1])
                new_circuit.cz(qargs[0],qargs[1])
                new_circuit.h(qargs[1])
                new_circuit.h(qargs[0])
                new_circuit.cz(qargs[1],qargs[0])
                new_circuit.h(qargs[0])
                new_circuit.h(qargs[1])
                new_circuit.cz(qargs[0],qargs[1])
                new_circuit.h(qargs[1])
        elif (inst.name == 'iswap' or inst.label == 'iswap'):
            new_circuit.unitary(iswap(), [qargs[0], qargs[1]], label='iswap')
        else:
            new_circuit.data.append((inst, qargs, cargs))
    return get_layer_circuits(new_circuit), len(decomp_iswap), len(decomp_cphase), len(cnot_gates), len(swap_gates)

def reschedule_layer(layers, coupling, verbose):
    if (len(layers) < 1):
        print("Too few layers for reschedule_layer")
        sys.exit()

    new_qregs = layers[0].qregs
    new_cregs = layers[0].cregs
    if len(new_qregs) > 1:
        print("Multiple q registers.")
    if len(new_cregs) == 0:
        new_circuit = QuantumCircuit(new_qregs[0])
    elif len(new_cregs) == 1:
        new_circuit = QuantumCircuit(new_qregs[0], new_cregs[0])
    else:
        print("Multiple c registers.")
        new_circuit = QuantumCircuit(new_qregs[0])

    num_layers = len(layers)
    #print("nl %d" % num_layers)
    idx = 1
    pending_g = []
    for g in layers[0].data:
        pending_g.append(g)
    while len(pending_g) > 0 or idx < num_layers:
        #print("idx %d" % idx)
        #print(coupling)
        active_q = [] # active qubits in this layer
        next_pd_g = [] # gates/instructions
        for inst, qargs, cargs in pending_g:
            if len(qargs) == 2:
                # two-qubit gates
                q1 = qargs[0].index
                q2 = qargs[1].index
                #print("Looking: (%d,%d)" % (q1, q2))
                #print("Active: ", active_q)
                conf = False
                if (q1 in active_q or q2 in active_q):
                    conf = True
                if not conf:
                    for (n1, n2) in coupling:
                        if (q1 == n1 and n2 in active_q) or (q2 == n1 and n2 in active_q) or (q1 == n2 and n1 in active_q) or (q2 == n2 and n1 in active_q):
                            conf = True
                            if verbose == 0:
                                print("Delay gate (%d, %d) due to crosstalk" % (q1,q2))
                            break
                if not conf:
                    active_q.append(q1)
                    active_q.append(q2)
                    new_circuit.data.append((inst, qargs, cargs))
                else:
                    next_pd_g.append((inst, qargs, cargs))


            elif len(qargs) == 1:
                # single-qubit gates
                q1 = qargs[0].index
                conf = False
                if (q1 in active_q):
                    conf = True
                if not conf:
                    active_q.append(q1)
                    new_circuit.data.append((inst, qargs, cargs))
                else:
                    next_pd_g.append((inst, qargs, cargs))
            else:
                # Multi-qubit gates; Barrier?
                numq = len(qargs)
                conf = False
                for qii in xrange(numq):
                    qi = qargs[qii].index
                    if (qi in active_q):
                        conf = True
                if not conf:
                    active_q.append(qi)
                    new_circuit.data.append((inst, qargs, cargs))
                else:
                    next_pd_g.append((inst, qargs, cargs))

        pending_g = next_pd_g
        if (len(pending_g) == 0 and idx < num_layers):
            for g in layers[idx].data:
                pending_g.append(g)
            idx += 1
        elif (idx < num_layers and not layers[idx].data[0][0].name == 'barrier'):
            for g in layers[idx].data:
                pending_g.append(g)
            idx += 1
        new_circuit.barrier(new_qregs[0])
    return get_layer_circuits(new_circuit)

def tiling_layer(layers, tilings, pattern_offset, verbose):
    if (len(layers) < 1):
        print("Too few layers for tiling_layer")
        sys.exit()

    new_qregs = layers[0].qregs
    new_cregs = layers[0].cregs
    if len(new_qregs) > 1:
        print("Multiple q registers.")
    if len(new_cregs) == 0:
        new_circuit = QuantumCircuit(new_qregs[0])
    elif len(new_cregs) == 1:
        new_circuit = QuantumCircuit(new_qregs[0], new_cregs[0])
    else:
        print("Multiple c registers.")
        new_circuit = QuantumCircuit(new_qregs[0])

    num_layers = len(layers)
    #print("nl %d" % num_layers)
    idx = 1
    pattern_idx = 0
    pending_g = []
    patterns = []
    for g in layers[0].data:
        pending_g.append(g)
    while len(pending_g) > 0 or idx < num_layers:
        #print("idx %d" % idx)
        #print(coupling)
        active_q = [] # active qubits in this layer
        next_pd_g = [] # gates/instructions
        two_q_layer = False
        two_q_tiling_delay = [] # 2q gates delayed due to tiling
        for inst, qargs, cargs in pending_g:
            if len(qargs) == 2:
                # two-qubit gates
                q1 = qargs[0].index
                q2 = qargs[1].index
                #print("Looking: (%d,%d)" % (q1, q2))
                #print("Active: ", active_q)
                conf = False
                if (q1 in active_q or q2 in active_q):
                    conf = True
                    next_pd_g.append((inst, qargs, cargs))
                    if (verbose == 0):
                        print("Delay gate (%d, %d)" % (q1,q2))
                #for (n1, n2) in coupling:
                #    if (q1 == n1 and n2 in active_q) or (q2 == n1 and n2 in active_q) or (q1 == n2 and n1 in active_q) or (q2 == n2 and n1 in active_q):
                #        conf = True
                #        if verbose==0:
                #            print("Delay gate (%d, %d)" % (q1,q2))
                #        break
                if not conf:
                    if (not((q1,q2) in tilings[(pattern_offset+pattern_idx) % 8] or (q2,q1) in tilings[(pattern_offset+pattern_idx) % 8])):
                        conf = True
                        two_q_tiling_delay.append((inst, qargs, cargs))
                        if verbose==0:
                            print("Delay gate (%d, %d) due to tiling" % (q1,q2))

                if not conf:
                    two_q_layer = True
                    active_q.append(q1)
                    active_q.append(q2)
                    new_circuit.data.append((inst, qargs, cargs))

            elif len(qargs) == 1:
                # single-qubit gates
                q1 = qargs[0].index
                conf = False
                if (q1 in active_q):
                    conf = True
                    if (verbose == 0):
                        print("Delay gate %d" % q1)
                if not conf:
                    active_q.append(q1)
                    new_circuit.data.append((inst, qargs, cargs))
                else:
                    next_pd_g.append((inst, qargs, cargs))
            else:
                # Multi-qubit gates; Barrier?
                numq = len(qargs)
                conf = False
                for qii in range(numq):
                    qi = qargs[qii].index
                    if (qi in active_q):
                        conf = True
                if not conf:
                    active_q.append(qi)
                    new_circuit.data.append((inst, qargs, cargs))
                else:
                    next_pd_g.append((inst, qargs, cargs))

        if two_q_layer == False and len(two_q_tiling_delay) > 0:
            two_q_possible = []
            for inst, qargs, cargs in two_q_tiling_delay:
                q1 = qargs[0].index
                q2 = qargs[1].index
                if (q1 in active_q or q2 in active_q):
                    next_pd_g.append((inst, qargs, cargs))
                else:
                    two_q_possible.append((inst, qargs, cargs))
            if len(two_q_possible) > 0:
                two_q_layer = True
                inst, qargs, cargs = two_q_possible[0]
                q1 = qargs[0].index
                q2 = qargs[1].index
                if (q1,q2) in tilings[pattern_offset % 8] or (q2,q1) in tilings[pattern_offset % 8]:
                    pattern_idx = 0
                elif (q1,q2) in tilings[(pattern_offset+1) % 8] or (q2,q1) in tilings[(pattern_offset+1) % 8]:
                    pattern_idx = 1
                elif (q1,q2) in tilings[(pattern_offset+2) % 8] or (q2,q1) in tilings[(pattern_offset+2) % 8]:
                    pattern_idx = 2
                elif (q1,q2) in tilings[(pattern_offset+3) % 8] or (q2,q1) in tilings[(pattern_offset+3) % 8]:
                    pattern_idx = 3
                elif (q1,q2) in tilings[(pattern_offset+4) % 8] or (q2,q1) in tilings[(pattern_offset+4) % 8]:
                    pattern_idx = 4
                elif (q1,q2) in tilings[(pattern_offset+5) % 8] or (q2,q1) in tilings[(pattern_offset+5) % 8]:
                    pattern_idx = 5
                elif (q1,q2) in tilings[(pattern_offset+6) % 8] or (q2,q1) in tilings[(pattern_offset+6) % 8]:
                    pattern_idx = 6
                elif (q1,q2) in tilings[(pattern_offset+7) % 8] or (q2,q1) in tilings[(pattern_offset+7) % 8]:
                    pattern_idx = 7
                for inst, qargs, cargs in two_q_possible:
                    q1 = qargs[0].index
                    q2 = qargs[1].index
                    if (q1,q2) in tilings[(pattern_offset+pattern_idx) % 8] or (q2,q1) in tilings[(pattern_offset+pattern_idx) % 8]:
                        if (q1 in active_q or q2 in active_q):
                            next_pd_g.append((inst, qargs, cargs))
                        else:
                            new_circuit.data.append((inst, qargs, cargs))
                            active_q.append(q1)
                            active_q.append(q2)
                    else:
                        next_pd_g.append((inst, qargs, cargs))

        pending_g = next_pd_g
        if (len(pending_g) == 0 and idx < num_layers):
            for g in layers[idx].data:
                pending_g.append(g)
            idx += 1
        elif (idx < num_layers and not layers[idx].data[0][0].name == 'barrier'):
            for g in layers[idx].data:
                pending_g.append(g)
            idx += 1
        patterns.append((pattern_idx+pattern_offset) % 8)
        new_circuit.barrier(new_qregs[0])
        patterns.append((pattern_idx+pattern_offset) % 8)
        if (two_q_layer):
            pattern_idx += 1
    return get_layer_circuits(new_circuit), (pattern_offset+pattern_idx) % 8, patterns




"""
Input: A dictionary with (key, value) = (edge, int color)
Return: A dictionary with (key, value) = (edge, int color), such that if a<b
then more nodes are colored with a than b
"""
def relabel_coloring(int_coloring):
    num_int = len(set(int_coloring.values()))
    arr = [] # [(color, popularity) for each color]
    for i in range(num_int):
        num = len([k for k,v in int_coloring.items() if v == i])
        arr.append((i,num))
    arr.sort(key=lambda x:x[1]) # sort arr by popularity
    new_coloring = {}
    for i in range(num_int):
        old_color = arr[num_int-i-1][0] # the new color is i
        # collect the list of nodes with color = old_color
        temp = [k for k,v in int_coloring.items() if v == old_color]
        # color them with the new color
        for k in temp:
            new_coloring[k] = i
    return new_coloring

def limit_colors(layers, lim_colors, G_crosstalk, verbose):
    # print("Limit colors:")
    if len(layers) < 1 and verbose == 0:
        print("Too few layers for limit_colors")
        sys.exit()
    new_qregs = layers[0].qregs
    new_cregs = layers[0].cregs
    if len(new_qregs) > 1 and verbose == 0:
        print("Multiple q registers.")
    if len(new_cregs) == 0:
        new_circuit = QuantumCircuit(new_qregs[0])
    elif len(new_cregs) == 1:
        new_circuit = QuantumCircuit(new_qregs[0], new_cregs[0])
    else:
        if verbose == 0:
            print("Multiple c registers.")
        new_circuit = QuantumCircuit(new_qregs[0])
    num_layers = len(layers)
    # print("number of layers:",num_layers)
    idx = 1
    pending_g = []
    for g in layers[0].data:
        pending_g.append(g)
    while len(pending_g) > 0 or idx < num_layers:
        # print("idx (of the next layer, start with 1):",idx)
        #print("idx %d" % idx)
        #print(coupling)
        next_pd_g = [] # gates/instructions
        edges = []
        for inst, qargs, cargs in pending_g:
            if len(qargs) == 2:
                # two-qubit gates
                q1 = qargs[0].index
                q2 = qargs[1].index
                edge = (q1, q2)
                edges.append(edge)
                edges.append((q2,q1))
                next_pd_g.append((inst, qargs, cargs))
            else: # 1-qubit gate
                new_circuit.data.append((inst, qargs, cargs))
        # put the gates in the next layer into pending_g
        pending_g = []
        if (idx < num_layers):
            for g in layers[idx].data:
                pending_g.append(g)
        # if there's at least 1 2-q gate
        if (len(next_pd_g) > 0):
            #print(" There are 2-q gates")
            # color the layer
            subgraph = nx.subgraph(G_crosstalk, edges)
            #print(G_crosstalk.nodes())
            #print(subgraph.nodes())
            int_coloring = nx.coloring.greedy_color(subgraph)
            #print("int_coloring:")
            #print(int_coloring)
            int_coloring = relabel_coloring(int_coloring)
            num_int = len(set(int_coloring.values())) # the number of colors
            # if the coloring uses more colors than allowed
            if num_int > lim_colors:
                #print(" Layer to be split, num_int:",num_int)
                # split the layer
                split = int(num_int/lim_colors) # the number of new layers
                if num_int%lim_colors != 0:
                    split += 1
                #print("number of sublayers:",split)
                # for all but the last of the new layers
                for i in range(split-1):
                    # append the gates that belong to the ith sublayer
                    for inst, qargs, cargs in next_pd_g:
                        q1 = qargs[0].index
                        q2 = qargs[1].index
                        try:
                            ic = int_coloring[(q1,q2)]
                        except:
                            ic = int_coloring[(q2,q1)]
                        if ic >= lim_colors*i and ic < lim_colors*(i+1):
                            new_circuit.data.append((inst, qargs, cargs))
                    new_circuit.barrier(new_qregs[0])
                    #print("Added sublayer:", i)
                # obtain the gates in the last sublayer
                last_layer = []
                for inst, qargs, cargs in next_pd_g:
                    q1 = qargs[0].index
                    q2 = qargs[1].index
                    try:
                        ic = int_coloring[(q1,q2)]
                    except:
                        ic = int_coloring[(q2,q1)]
                    if ic >= (split-1)*lim_colors:
                        last_layer.append((inst, qargs, cargs))
                #print("Last sublayer obtained, len:",len(last_layer))
                # for the last sublayer, we check if it can be combined with the next layer
                if idx >= num_layers:
                    #print("Already the last layer, just append")
                    for g in last_layer:
                        new_circuit.data.append(g)
                    new_circuit.barrier(new_qregs[0])
                else:
                    conf = False
                    #print("check if can be combined with the next layer")
                    active_next = [] # qubits used in the next layer
                    for inst, qargs, cargs in pending_g:
                        q1 = qargs[0].index
                        active_next.append(q1)
                        if len(qargs) == 2:
                            q2 = qargs[1].index
                            active_next.append(q2)
                    #print("active_next:",active_next)
                    for inst, qargs, cargs in last_layer:
                        q1 = qargs[0].index
                        q2 = qargs[1].index
                        if q1 in active_next or q2 in active_next:
                            conf = True
                            break
                    if conf: # the last sublayer cannot be combined with the next layer
                        #print("cannot be combined")
                        for g in last_layer:
                            new_circuit.data.append(g)
                        new_circuit.barrier(new_qregs[0])
                    else: # the last sublayer can be combined with the next layer
                        #print("can be combined")
                        for g in last_layer:
                            pending_g.append(g)
            else: # if the layer doesn't need to be split
                #print(" Layer doesn't need to be split")
                for inst, qargs, cargs in next_pd_g:
                    new_circuit.data.append((inst, qargs, cargs))
                new_circuit.barrier(new_qregs[0])
        else: # this layer only has 1-q gates
            #print(" There aren't 2-q gates")
            new_circuit.barrier(new_qregs[0])
        idx += 1
    return get_layer_circuits(new_circuit)

def get_iswap_time(edges, qubit_freqs, tau_special):
    gt = 0.0 # in ns
    Cqq = CQQ
    for (q1, q2) in edges:
        if (q1,q2) in tau_special:
            c = tau_special[(q1,q2)]
        elif (q2,q1) in tau_special:
            c = tau_special[(q2,q1)]
        else:
            c = 1.0
        # iswap time: pi/2g
        J = 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq
        t = c * np.pi/(2 * J)
        if t > gt: gt = t
    return gt

def get_cphase_time(edges, qubit_freqs, tau_special):
    gt = 0.0 # in ns
    Cqq = CQQ
    alphas = [ALPHA for f in qubit_freqs] #TODO
    for (q1, q2) in edges:
        if (q1,q2) in tau_special:
            c = tau_special[(q1,q2)]
        elif (q2,q1) in tau_special:
            c = tau_special[(q2,q1)]
        else:
            c = 1.0
        # cphase time: pi/sqrt(2)g
        J = 0.5 * np.sqrt(qubit_freqs[q1] * (qubit_freqs[q2]+alphas[q2])) * Cqq
        t = c * np.pi/(np.sqrt(2) * J)
        if t > gt: gt = t
    return gt

def get_max_time(gt, taus):
    tmax = gt
    for t in taus.values():
        if t > tmax: tmax = t
    return tmax

def success_rate_rand_coloring(device, circuit, scheduler, d, decomp):
    success_rate = 0.0
    avg_success = 0.0
    worst_success = 0.0
    depth_before = 0
    depth_after = 0
    total_time = 0.0
    max_colors = 0
    t_act = []
    t_2q = []
    return success_rate, avg_success, worst_success, depth_before, depth_after, total_time, max_colors, t_act, t_2q



def success_rate_full_coloring(device, circuit, scheduler, d, decomp, outfile, verbose):
    outf = h5py.File(outfile + ".hdf5", "w")
    # output data stores the information needed for hamiltonian simulation
    # Groups: freqs and gates (for each timestep)
    # Dataset format: [list of freqs for each timestep]
    #                 groups of [('gatename',list of qubits)]
    freqsdata = []
    gatesdata = []
    width = device.side_length
    height = device.side_length
    num_q = width * height
    omega_max = device.omega_max
    delta_int = device.delta_int
    delta_ext= device.delta_ext
    delta_park = device.delta_park

    G_connect = get_connectivity_graph(width, height)
    park_coloring = nx.coloring.greedy_color(G_connect)
    num_park = len(set(park_coloring.values()))
    G_crosstalk = get_aug_line_graph(width, height, d)
    if full_uni:
        int_coloring = {}
        for node in G_crosstalk.nodes:
            int_coloring[node] = 0
    else:
        int_coloring = nx.coloring.greedy_color(G_crosstalk)
    #print(int_coloring)
    num_int = len(set(int_coloring.values()))
    def _build_color_map():
        # negative colors for parking, non-negative colors for interaction.
        step_park = delta_park / num_park
        step_int = delta_int / num_int
        colormap = dict()
        for c in range(num_int):
            colormap[str(c)] = omega_max - c * step_int
        for c in range(num_park):
            colormap[str(-(c+1))]= omega_max - delta_int - delta_ext - c * step_park
        return colormap
    color_to_freq = _build_color_map()
    def _park_freq(c):
        omg = color_to_freq[str(-(c+1))]
        return omg + get_flux_noise(omg, sigma)#+np.random.normal(0,sigma)
    def _int_freq(c):
        omg = color_to_freq[str(c)]
        return omg + get_flux_noise(omg, sigma)#+np.random.normal(0,sigma)
    def _initial_frequency():
        freqs = dict()
        for q in range(width*height):
            freqs[q] = _park_freq(park_coloring[q])
        return freqs

    t_act = np.zeros(num_q)
    t_2q = np.zeros(num_q)
    active_list = [False for i in range(num_q)]
    success_rate = 1.0
    single_qb_err = 0.0015
    tot_success = 0.0
    tot_cnt = 0
    worst_success = 1.0
    max_colors = num_int # max number of colors used
    # Mapping
    coupling = get_nearest_neighbor_coupling_list(width, height)
    circ_mapped = get_map_circuit(circuit, coupling)

    Cqq = CQQ

    #circuit.draw(output='mpl')
    # Check scheduler
    if (scheduler == 'hybrid'):
        print("Hybrid scheduler to be implemented.")
        sys.exit(2)
    else:
        # scheduler == qiskit or greedy
        layers = get_layer_circuits(circ_mapped)
        qubit_freqs = _initial_frequency()
        alphas = [ALPHA for f in qubit_freqs] #TODO
        #leftover = []
        #left_gates= []
        num_layers = len(layers)
        idx = 0
        total_time = 0.0 # in ns
        total_tcz = 0.0 # total time spent on CZ gates
        # total number of gates decomposed with iswap and cphase, used only if d=flexible
        num_iswap = 0
        num_cphase = 0
        num_cnot = 0
        num_swap = 0
        #while (idx < num_layers or len(leftover) > 0):
        for idx in range(num_layers):
            all_gates = []
            layer_circuit = layers[idx]
            if verbose == 0:
                print(layer_circuit)
            print(idx, "-----------------")
            #print('leftover len: ' + str(len(leftover)))
            #print('left_gates len: ' + str(len(left_gates)))
            #decomp_layer = decompose_layer(layer_circuit, decomp)
            if decomp == 'flexible':
                if verbose == 0:
                    print("Decompose layer", idx)
                decomp_layer, temp_iswap, temp_cphase, temp_cnot, temp_swap = decompose_layer_flexible(layer_circuit, G_crosstalk, verbose)
                num_iswap += temp_iswap
                num_cphase += temp_cphase
                num_cnot += temp_cnot
                num_swap += temp_swap
            else:
                decomp_layer = decompose_layer(layer_circuit, decomp)
            if (scheduler == 'greedy'):
                resched_layer = reschedule_layer(decomp_layer, coupling, verbose)
            else:
                resched_layer = decomp_layer
            #print("Layers: %d %d" % (len(decomp_layer), len(resched_layer)))
            for layer in resched_layer:
                print(layer)
                edges = []
                edges_cphase = [] # if (q1,q2) is in it, then (q2,q1) is also in it
                edges_iswaps = [] # if (q1,q2) is in it, then (q2,q1) is also in it
                #edges = [e for e in leftover]
                #curr_gates = [e for e in left_gates]
                curr_gates = []
                taus = {} # For storing couplings that needed sqrt iswaps
                gt = 0.0
                layer_time = 0.0 # ns
                barrier = False
                single_qb_err = 0.0015
                single_qb_err_acc = 1.0
                for g, qargs, _ in layer.data:
                    #print(g, qargs)
                    #print(g.qasm())
                    if g.name == "barrier": barrier = True
                    if g.name == "measure": continue
                    if len(qargs) == 1:
                        all_gates.append((g.qasm(),(qargs[0].index, -1)))
                        active_list[qargs[0].index] = True
                        gt = GATETIMES[g.name]
                        if gt > layer_time: layer_time = gt
                        single_qb_err_acc *= 1 - single_qb_err
                    if len(qargs) == 2:
                        q1, q2 = qargs[0].index, qargs[1].index
                        active_list[q1] = True
                        active_list[q2] = True
                        edges.append((q1, q2))
                        #print((q1,q2))
                        curr_gates.append((g.qasm(),(q1, q2)))
                        try:
                            ic = int_coloring[(q1,q2)]
                        except:
                            ic = int_coloring[(q2,q1)]
                        f = _int_freq(ic)
                        #print('freq:',f)
                        if (g.name == 'unitary' and g.label == 'iswap'):
                            qubit_freqs[q1] = f
                            qubit_freqs[q2] = f
                            taus[(q1,q2)] = np.pi / (2 * 0.5 * np.sqrt(f*f) * Cqq)
                            edges_iswaps.append((q1,q2))
                        elif (g.name == 'unitary' and g.label == 'sqrtiswap'):
                            qubit_freqs[q1] = f
                            qubit_freqs[q2] = f
                            #tau_special[(q1,q2)] = 0.5
                            taus[(q1,q2)] = 0.5 * np.pi / (2 * 0.5 * np.sqrt(f*f) * Cqq)
                            edges_iswaps.append((q1,q2))
                        elif (g.name == 'cz' or g.label == 'cz'):
                            qubit_freqs[q1] = f
                            qubit_freqs[q2] = f - alphas[q2] # b/c match f1 with f2+alpha
                            taus[(q1,q2)] = np.pi / (np.sqrt(2) * 0.5 * np.sqrt(f*f) * Cqq) # f is interaction freq
                            edges_cphase.append((q1,q2))
                        else:
                            print("Gate %s(%s) not recognized. Supports iswap, sqrtiswap, cz." % (g.name, g.label))
                        t_2q[q1] += taus[(q1,q2)]
                        t_2q[q2] += taus[(q1,q2)]
                success_rate *= single_qb_err_acc
                #if (scheduler == 'greedy'):
                #    edges, leftover, ind = greedy_reschedule(coupling, edges)
                #    for i in range(len(ind)):
                #        if (ind[i]):
                #            all_gates.append(curr_gates[i])
                #        else:
                #            left_gates.append(curr_gates[i])
                #else:
                #    for i in range(len(curr_gates)):
                #        all_gates.append(curr_gates[i])
                for i in range(len(curr_gates)):
                    all_gates.append(curr_gates[i])


                #print("all_gates:")
                #print(all_gates)
                #print("edges:")
                #print(edges)
                #print("leftover:")
                #print(leftover)
                #print("qubit_freqs:")
                #print(qubit_freqs)
                #err = compute_crosstalk(edges, coupling, qubit_freqs)
                if barrier:
                    err, swap_err, leak_err = 0.0, 0.0, 0.0
                    gt = 0.0
                else:
                    if decomp == 'iswap':
                        #err, swap_err, leak_err = compute_crosstalk_iswaps(edges, coupling, qubit_freqs, tau_special)
                        #gt = get_iswap_time(edges, qubit_freqs, tau_special)
                        gt = get_max_time(gt, taus)
                        err, swap_err, leak_err = compute_crosstalk_iswaps(edges_iswaps, coupling, qubit_freqs, taus, gt)
                    elif decomp == 'cphase':
                        gt = get_max_time(gt, taus)
                        err, swap_err, leak_err = compute_crosstalk_cphase(edges_iswaps, coupling, qubit_freqs, taus, gt)
                        #err, swap_err, leak_err = compute_crosstalk_cphase(edges, coupling, qubit_freqs, tau_special)
                        #gt = get_cphase_time(edges, qubit_freqs, tau_special)
                    elif decomp == 'flexible' or decomp == 'mixed':
                        gt = get_max_time(gt, taus)
                        err, swap_err, leak_err = compute_crosstalk_flexible(edges_cphase, edges_iswaps, coupling, qubit_freqs, taus, gt)
                    else:
                        print("Decomposition method %s not recognized. Try iswap or cphase." % decomp)
                        gt = 0.0
                if gt > layer_time: layer_time = gt
                for qubit in range(num_q):
                    if active_list[qubit]:
                        t_act[qubit] += layer_time
                update_data(freqsdata, gatesdata, qubit_freqs, all_gates, num_q)
                success_rate *= 1 - err
                if verbose == 0:
                    print("Layer success: %12.10f (swap: %12.10f, leak: %12.10f)" % (1-err, swap_err, leak_err))
                    print("Layer time:", layer_time)
                # Reset the frequencies
                for (q1, q2) in edges:
                    qubit_freqs[q1] = _park_freq(park_coloring[q1])
                    qubit_freqs[q2] = _park_freq(park_coloring[q2])
                if not barrier:
                    tot_success += (1 - err)*single_qb_err_acc
                    tot_cnt += 1
                worst_success = min(worst_success, 1 - err)
                total_time += layer_time
            #idx += 1
        write_circuit(freqsdata, gatesdata, outf)
    avg_success = tot_success / tot_cnt
    if decomp=='flexible' and verbose == 0:
        print("Total number of CNOT gates:", num_cnot)
        print("Total number of SWAP gates:", num_swap)
        print("Total number of 2-qubit gates that are decomposed with iSWAP:", num_iswap)
        if num_iswap+num_cphase > 0:
            print("Proportion:",num_iswap/(num_iswap+num_cphase))
        print("Total number of 2-qubit gates that are decomposed with CPHASE:", num_cphase)
        if num_iswap+num_cphase > 0:
            print("Proportion:",num_cphase/(num_iswap+num_cphase))
    return success_rate, avg_success, worst_success, idx, tot_cnt, total_time, max_colors, t_act, t_2q


def success_rate_layer_coloring(device, circuit, scheduler, d, decomp, outfile, lim_colors, verbose):
    outf = h5py.File(outfile + ".hdf5", "w")
    # output data stores the information needed for hamiltonian simulation
    # Groups: freqs and gates (for each timestep)
    # Dataset format: [list of freqs for each timestep]
    #                 [('gatename',list of qubits)]
    freqsdata = []
    gatesdata = []
    width = device.side_length
    height = device.side_length
    num_q = width * height
    omega_max = device.omega_max
    delta_int = device.delta_int
    delta_ext= device.delta_ext
    delta_park = device.delta_park

    G_connect = get_connectivity_graph(width, height)
    park_coloring = nx.coloring.greedy_color(G_connect)
    num_park = len(set(park_coloring.values()))
    G_crosstalk = get_aug_line_graph(width, height, d)
    coupling = get_nearest_neighbor_coupling_list(width, height)
    Cqq = CQQ

    def _build_park_color_map():
        # negative colors for parking, non-negative colors for interaction.
        step_park = delta_park / num_park
        #step_int = delta_int / num_int
        colormap = dict()
        #for c in range(num_int):
        #    colormap[str(c)] = omega_max - c * step_int
        for c in range(num_park):
            colormap[str(-(c+1))]= omega_max - delta_int - delta_ext - c * step_park
        return colormap
    def _add_int_color_map(colormap, n_int):
        step_int = delta_int / n_int
        for c in range(n_int):
            colormap[str(c)] = omega_max - c * step_int

    color_to_freq = _build_park_color_map()
    def _park_freq(c):
        omg = color_to_freq[str(-(c+1))]#+np.random.normal(0,sigma)
        return omg + get_flux_noise(omg, sigma)#+np.random.normal(0,sigma)
    def _initial_frequency():
        freqs = dict()
        for q in range(width*height):
            freqs[q] = _park_freq(park_coloring[q])
        return freqs
    circ_mapped = get_map_circuit(circuit, coupling)
    #circuit.draw(output='mpl')
    t_act = np.zeros(num_q)
    t_2q = np.zeros(num_q)
    active_list = [False for i in range(num_q)]
    success_rate = 1.0
    tot_success = 0.0
    tot_cnt = 0
    max_colors = 0 # max number of colors used
    worst_success = 1.0
    # Check scheduler
    if (scheduler == 'hybrid'):
        print("Hybrid scheduler to be implemented.")
        sys.exit(2)
    else:
        #leftover = []
        #left_gates= []
        layers = get_layer_circuits(circ_mapped)
        qubit_freqs = _initial_frequency()
        num_layers = len(layers)
        idx = 0
        total_time = 0.0 # ns
        total_tcz = 0.0
        alphas = [ALPHA for f in qubit_freqs] #TODO

        if verbose == 0:
            print("Num of layers:", num_layers)

        # total number of gates decomposed with iswap and cphase, used only if d=flexible
        num_iswap = 0
        num_cphase = 0
        num_cnot = 0
        num_swap = 0

        #while (idx < num_layers or len(leftover) > 0):
        for idx in range(num_layers):
            all_gates = []
            layer_circuit = layers[idx]
            if verbose == 0:
                print(layer_circuit)
            print(idx, "-----------------")
            if decomp == 'flexible':
                if verbose == 0:
                    print("Decompose layer", idx)
                decomp_layer, temp_iswap, temp_cphase, temp_cnot, temp_swap = decompose_layer_flexible(layer_circuit, G_crosstalk, verbose)
                num_iswap += temp_iswap
                num_cphase += temp_cphase
                num_cnot += temp_cnot
                num_swap += temp_swap
            else:
                decomp_layer = decompose_layer(layer_circuit, decomp)
            if (scheduler == 'greedy'):
                resched_layer = reschedule_layer(decomp_layer, coupling, verbose)
            else:
                resched_layer = decomp_layer
            if lim_colors > 0:
                resched_layer = limit_colors(resched_layer, lim_colors, G_crosstalk, verbose)
            for layer in resched_layer:
                print(layer.qasm())
                # Pre-fill edges for constructing (undirected) xtalk graph
                #edges = [leftover[i//2] if i%2==0 else (leftover[i//2][1], leftover[i//2][0]) for i in range(2*len(leftover))]
                edges = []
                edges_cphase = [] # if (q1,q2) is in it, then (q2,q1) is also in it
                edges_iswaps = [] # if (q1,q2) is in it, then (q2,q1) is also in it
                #curr_gates = [e for e in left_gates]
                #leftover = []
                #left_gates = []
                taus = {} # For storing coupling times
                gt = 0.0
                layer_time = 0.0
                barrier = False
                for _, qargs, _ in layer.data:
                    if len(qargs) == 2:
                        q1, q2 = qargs[0].index, qargs[1].index
                        edge = (q1, q2)
                        edges.append(edge)
                        edges.append((q2,q1)) # because undirected graph
                #print("+Edges:")
                #print(edges)
                if (len(edges) > 0):
                    #idx += 1
                    #continue
                    subgraph = nx.subgraph(G_crosstalk, edges)
                    #print(G_crosstalk.nodes())
                    #print(subgraph.nodes())
                    int_coloring = nx.coloring.greedy_color(subgraph)
                    #print("+int_coloring:")
                    #print(int_coloring)
                    num_int = len(set(int_coloring.values()))
                    while lim_colors > 0 and num_int > lim_colors:
                        # need to recolor the layer, cannot use greedy_color because it is not seeded
                        # int_coloring = nx.coloring.greedy_color(subgraph,strategy='random_sequential')
                        int_coloring = {}
                        nodes = list(subgraph)
                        np.random.shuffle(nodes)
                        for u in nodes:
                            # Set to keep track of colors of neighbours
                            neighbour_colors = {int_coloring[v] for v in subgraph[u] if v in int_coloring}
                            # Find the first unused color.
                            temp_color = 0
                            while temp_color in neighbour_colors:
                                temp_color += 1
                            # Assign the new color to the current node.
                            int_coloring[u] = temp_color
                        num_int = len(set(int_coloring.values()))
                    if verbose == 0:
                        print("num_int: ", num_int)
                    int_coloring = relabel_coloring(int_coloring)
                    if num_int > max_colors: max_colors = num_int
                    #TODO?
                    if num_int == 0:
                        idx += 1
                        continue
                    _add_int_color_map(color_to_freq, num_int)
                def _int_freq(c):
                    omg =  color_to_freq[str(c)]#+np.random.normal(0,sigma)
                    return omg + get_flux_noise(omg, sigma)
                #print(layer)
                #print("-----------------")
                # Refill edges and curr_gates
                #edges = [e for e in leftover]
                edges = []
                #curr_gates = [e for e in left_gates]
                curr_gates = []
                single_qb_err = 0.0015
                single_qb_err_acc = 1.0
                for g, qargs, _ in layer.data:
                    if g.name == "barrier": barrier = True
                    if g.name == "measure": continue
                    #print(qargs)
                    #print(qargs[0].index)
                    if len(qargs) == 1: # single qubit gates
                        all_gates.append((g.qasm(),(qargs[0].index, -1)))
                        active_list[qargs[0].index] = True
                        gt = GATETIMES[g.name]
                        if gt > layer_time: layer_time = gt
                        single_qb_err_acc *= 1 - single_qb_err
                    elif len(qargs) == 2:
                        q1, q2 = qargs[0].index, qargs[1].index
                        active_list[q1] = True
                        active_list[q2] = True
                        edges.append((q1, q2))
                        curr_gates.append((g.qasm(),(q1, q2)))
                        try:
                            f = _int_freq(int_coloring[(q1, q2)])
                        except:
                            f = _int_freq(int_coloring[(q2, q1)])
                        #qubit_freqs[q1] = f
                        #qubit_freqs[q2] = f
                        if (g.name == 'unitary' and g.label == 'iswap'):
                            qubit_freqs[q1] = f
                            qubit_freqs[q2] = f
                            taus[(q1,q2)] = np.pi / (2 * 0.5 * np.sqrt(f*f) * Cqq)
                            edges_iswaps.append((q1,q2))
                        elif (g.name == 'unitary' and g.label == 'sqrtiswap'):
                            qubit_freqs[q1] = f
                            qubit_freqs[q2] = f
                            taus[(q1,q2)] = 0.5 * np.pi / (2 * 0.5 * np.sqrt(f*f) * Cqq)
                            edges_iswaps.append((q1,q2))
                        elif (g.name == 'cz' or g.label == 'cz'):
                            qubit_freqs[q1] = f
                            qubit_freqs[q2] = f - alphas[q2] # b/c match f1 with f2+alpha
                            taus[(q1,q2)] = np.pi / (np.sqrt(2) * 0.5 * np.sqrt(f*f) * Cqq) # f is interaction freq
                            edges_cphase.append((q1,q2))
                        else:
                            print("Gate %s(%s) not recognized. Supports iswap, sqrtiswap, cz." % (g.name, g.label))
                        t_2q[q1] += taus[(q1,q2)]
                        t_2q[q2] += taus[(q1,q2)]
                success_rate *= single_qb_err_acc
                #if (scheduler == 'greedy'):
                #    edges, leftover, ind = greedy_reschedule(coupling, edges)
                #    for i in range(len(ind)):
                #        if (ind[i]):
                #            all_gates.append(curr_gates[i])
                #        else:
                #            left_gates.append(curr_gates[i])
                #else:
                #    for i in range(len(curr_gates)):
                #        all_gates.append(curr_gates[i])
                for i in range(len(curr_gates)):
                    all_gates.append(curr_gates[i])

                #print("edges:")
                #print(edges)
                print("qubit_freqs:")
                print(qubit_freqs)
                #err = compute_crosstalk(edges, coupling, qubit_freqs)
                if barrier:
                    err, swap_err, leak_err = 0.0, 0.0, 0.0
                    gt = 0.0
                else:
                    if decomp == 'iswap':
                        gt = get_max_time(gt, taus)
                        err, swap_err, leak_err = compute_crosstalk_iswaps(edges_iswaps, coupling, qubit_freqs, taus, gt)
                        #gt = get_iswap_time(edges, qubit_freqs, taus)
                    elif decomp == 'cphase':
                        gt = get_max_time(gt, taus)
                        err, swap_err, leak_err = compute_crosstalk_cphase(edges_iswaps, coupling, qubit_freqs, taus, gt)
                        #gt = get_cphase_time(edges, qubit_freqs, taus)
                    elif decomp == 'flexible' or decomp == 'mixed':
                        gt = get_max_time(gt, taus)
                        err, swap_err, leak_err = compute_crosstalk_flexible(edges_cphase, edges_iswaps, coupling, qubit_freqs, taus, gt)
                        #gt_iswap = get_iswap_time(edges_iswaps, qubit_freqs, tau_special)
                        #gt_cphase = get_cphase_time(edges_cphase, qubit_freqs, tau_special)
                        #gt = max(gt_iswap, gt_cphase)
                    else:
                        print("Decomposition method %s not recognized. Try iswap, cphase, or flexible." % decomp)
                        gt = 0.0
                if gt > layer_time: layer_time = gt
                for qubit in range(num_q):
                    if active_list[qubit]:
                        t_act[qubit] += layer_time
                if verbose == 0:
                    print("Layer success: %12.10f (swap: %12.10f, leak: %12.10f)" % (1-err, swap_err, leak_err))
                    print("Layer time:", layer_time)
                update_data(freqsdata, gatesdata, qubit_freqs, all_gates, num_q)
                success_rate *= 1 - err
                # Reset the frequencies
                for (q1, q2) in edges:
                    qubit_freqs[q1] = _park_freq(park_coloring[q1])
                    qubit_freqs[q2] = _park_freq(park_coloring[q2])
                if not barrier:
                    tot_success += (1 - err)*single_qb_err_acc
                    tot_cnt += 1
                worst_success = min(worst_success, 1 - err)
                total_time += layer_time
            idx += 1
        write_circuit(freqsdata, gatesdata, outf)

    avg_success = tot_success / tot_cnt

    if decomp=='flexible' and verbose == 0:
        print("Total number of CNOT gates:", num_cnot)
        print("Total number of SWAP gates:", num_swap)
        print("Total number of 2-qubit gates that are decomposed with iSWAP:", num_iswap)
        if num_iswap+num_cphase > 0:
            print("Proportion:",num_iswap/(num_iswap+num_cphase))
        print("Total number of 2-qubit gates that are decomposed with CPHASE:", num_cphase)
        if num_iswap+num_cphase > 0:
            print("Proportion:",num_cphase/(num_iswap+num_cphase))

    return success_rate, avg_success, worst_success, idx, tot_cnt, total_time, max_colors, t_act, t_2q


def gen_tiling_pattern(device):
    # Generate coupler activation pattern, assume 2D grid of qubits
    width = device.side_length
    height = device.side_length
    num_q = width * height
    patternA = []
    patternB = []
    patternC = []
    patternD = []

    def _qubit_number(row, col):
        return row * width + col

    # horizontal edges
    for row in range(height):
        for col in range(width - 1):
            if (col+row)%2 == 1:
                patternC.append((_qubit_number(row, col), _qubit_number(row, col + 1)))
            else:
                patternD.append((_qubit_number(row, col), _qubit_number(row, col + 1)))
            #coupling_list.append((_qubit_number(row, col), _qubit_number(row, col + 1)))
            #if directed:
            #    coupling_list.append((_qubit_number(row, col + 1), _qubit_number(row, col)))

    # vertical edges
    for col in range(width):
        for row in range(height - 1):
            if (col+row)%2 == 1:
                patternB.append((_qubit_number(row, col), _qubit_number(row + 1, col)))
            else:
                patternA.append((_qubit_number(row, col), _qubit_number(row + 1, col)))
            #coupling_list.append((_qubit_number(row, col), _qubit_number(row + 1, col)))
            #if directed:
            #    coupling_list.append((_qubit_number(row + 1, col), _qubit_number(row, col)))

    return [patternA,patternB,patternC,patternD,patternC,patternD,patternA,patternB]




def success_rate_google_like(device, circuit, scheduler, d, decomp, outfile, verbose):
    # Full coloring with tunable transmon
    # Called when freq=google, also supporting scheduler=tiling
    outf = h5py.File(outfile + ".hdf5", "w")
    # output data stores the information needed for hamiltonian simulation
    # Groups: freqs and gates (for each timestep)
    # Dataset format: [list of freqs for each timestep]
    #                 groups of [('gatename',list of qubits)]
    freqsdata = []
    gatesdata = []
    width = device.side_length
    height = device.side_length
    num_q = width * height
    omega_max = device.omega_max
    delta_int = device.delta_int
    delta_ext= device.delta_ext
    delta_park = device.delta_park

    tilings = gen_tiling_pattern(device)
    G_connect = get_connectivity_graph(width, height)
    # park_coloring = nx.coloring.greedy_color(G_connect)
    # num_park = len(set(park_coloring.values()))
    G_crosstalk = get_aug_line_graph(width, height, d)
    int_freqs = {}
    if num_q == 4:
        int_freqs = {(0,1):6.619, (0,2):6.65, (1,3):6.657, (2,3):6.677}
        park_freqs = {0:6.605,1:6.638,2:6.694,3:6.681}
    elif num_q == 9:
        park_freqs = {0:6.605,1:6.638,2:6.565,3:6.694,4:6.681,5:6.601,6:6.643,7:6.621,8:6.646}
        int_freqs = {(0,1):6.619, (1,2):6.601, (3,4):6.677, (4,5):6.635, (6,7):6.631, (7,8):6.646, (0,3):6.65, (1,4):6.657, (2,5):6.585, (3,6):6.667, (4,7):6.645,(5,8):6.646}
    elif num_q == 16:
        park_freqs = {0:6.605,1:6.638,2:6.565,3:6.555,4:6.694,5:6.681,6:6.601,7:6.626,8:6.643,9:6.621,10:6.646,11:6.657,12:6.712,13:6.671,14:6.586,15:6.623}
        int_freqs = {(0,1):6.619,(1,2):6.601,(2,3):6.565,(4,5):6.677,(5,6):6.635,(6,7):6.595}
        int_freqs[(8,9)] = 6.631
        int_freqs[(9,10)] = 6.646
        int_freqs[(10,11)] = 6.646
        int_freqs[(12,13)] = 6.69
        int_freqs[(13,14)] = 6.631
        int_freqs[(14,15)] = 6.623
        int_freqs[(0,4)] = 6.65
        int_freqs[(1,5)] = 6.657
        int_freqs[(2,6)] = 6.585
        int_freqs[(3,7)] = 6.592
        int_freqs[(4,8)] = 6.667
        int_freqs[(5,9)] = 6.645
        int_freqs[(6,10)] = 6.646
        int_freqs[(7,11)] = 6.642
        int_freqs[(8,12)] = 6.68
        int_freqs[(9,13)] = 6.645
        int_freqs[(10,14)] = 6.646
        int_freqs[(11,15)] = 6.633
    num_int = len(set(int_freqs.values()))
    # num_park = len(set(park_freqs.values()))

    #print(int_coloring)
    # num_int = len(set(int_coloring.values()))
    # def _build_color_map():
        # negative colors for parking, non-negative colors for interaction.
    #    step_park = delta_park / num_park
    #    step_int = delta_int / num_int
    #    colormap = dict()
    #    for c in range(num_int):
    #        colormap[str(c)] = omega_max - c * step_int
    #    for c in range(num_park):
    #        colormap[str(-(c+1))]= omega_max - delta_int - delta_ext - c * step_park
    #    return colormap
    # color_to_freq = _build_color_map()
    # def _park_freq(c):
    #    omg = color_to_freq[str(-(c+1))]#+np.random.normal(0,sigma)
    #    return omg + get_flux_noise(omg, sigma)
    # def _int_freq(c):
    #    omg = color_to_freq[str(c)]#+np.random.normal(0,sigma)
    #    return omg + get_flux_noise(omg, sigma)
    def _initial_frequency():
        freqs = dict()
        for q in range(width*height):
            freqs[q] = park_freqs[q]
        return freqs

    t_act = np.zeros(num_q)
    t_2q = np.zeros(num_q)
    active_list = [False for i in range(num_q)]
    success_rate = 1.0
    single_qb_err = 0.0015
    tot_success = 0.0
    tot_cnt = 0
    worst_success = 1.0
    max_colors = num_int # max number of colors used
    # Mapping
    coupling = get_nearest_neighbor_coupling_list(width, height)
    circ_mapped = get_map_circuit(circuit, coupling)

    Cqq = CQQ

    #circuit.draw(output='mpl')
    # Check scheduler
    if (scheduler == 'hybrid'):
        print("Hybrid scheduler to be implemented.")
        sys.exit(2)
    else:
        # scheduler == qiskit or greedy
        layers = get_layer_circuits(circ_mapped)
        qubit_freqs = _initial_frequency()
        alphas = [ALPHA for f in qubit_freqs] #TODO
        #leftover = []
        #left_gates= []
        num_layers = len(layers)
        idx = 0
        pattern_offset = 0
        total_time = 0.0 # in ns
        # total number of gates decomposed with iswap and cphase, used only if d=flexible
        num_iswap = 0
        num_cphase = 0
        num_cnot = 0
        num_swap = 0
        #while (idx < num_layers or len(leftover) > 0):

        tiling_idx = 0
        for idx in range(num_layers):
            all_gates = []
            layer_circuit = layers[idx]
            if verbose == 0:
                print(layer_circuit)
            print(idx, "-----------------")
            #print('leftover len: ' + str(len(leftover)))
            #print('left_gates len: ' + str(len(left_gates)))
            #decomp_layer = decompose_layer(layer_circuit, decomp)
            if decomp == 'flexible':
                if verbose == 0:
                    print("Decompose layer", idx)
                decomp_layer, temp_iswap, temp_cphase, temp_cnot, temp_swap = decompose_layer_flexible(layer_circuit, G_crosstalk, verbose)
                num_iswap += temp_iswap
                num_cphase += temp_cphase
                num_cnot += temp_cnot
                num_swap += temp_swap
            else:
                decomp_layer = decompose_layer(layer_circuit, decomp)
            if (scheduler == 'greedy'):
                resched_layer = reschedule_layer(decomp_layer, coupling, verbose)
            elif (scheduler == 'tiling'):
                resched_layer, pattern_offset, patterns = tiling_layer(decomp_layer, tilings, pattern_offset, verbose)
                print(patterns)
            else:
                resched_layer = decomp_layer
            #print("Layers: %d %d" % (len(decomp_layer), len(resched_layer)))

            for (iid,layer) in enumerate(resched_layer):
                if (scheduler == 'tiling'):
                    tiling_idx = patterns[iid]
                    print("tiling_idx", tiling_idx)
                print(layer.qasm())
                edges = []
                edges_cphase = [] # if (q1,q2) is in it, then (q2,q1) is also in it
                edges_iswaps = [] # if (q1,q2) is in it, then (q2,q1) is also in it
                #edges = [e for e in leftover]
                #curr_gates = [e for e in left_gates]
                curr_gates = []
                taus = {} # For storing couplings that needed sqrt iswaps
                gt = 0.0
                layer_time = 0.0 # ns
                barrier = False
                single_qb_err = 0.0015
                for g, qargs, _ in layer.data:
                    #print(g, qargs)
                    #print(g.qasm())
                    if g.name == "barrier": barrier = True
                    if g.name == "measure": continue
                    if len(qargs) == 1:
                        all_gates.append((g.qasm(),(qargs[0].index, -1)))
                        active_list[qargs[0].index] = True
                        gt = GATETIMES[g.name]
                        if gt > layer_time: layer_time = gt
                        success_rate *= 1 - single_qb_err
                    if len(qargs) == 2:
                        q1, q2 = qargs[0].index, qargs[1].index
                        active_list[q1] = True
                        active_list[q2] = True
                        edges.append((q1, q2))
                        #print((q1,q2))
                        curr_gates.append((g.qasm(),(q1, q2)))
                        try:
                            f = int_freqs[(q1,q2)]
                        except:
                            f = int_freqs[(q2,q1)]
                        #print('freq:',f)
                        if (g.name == 'unitary' and g.label == 'iswap'):
                            qubit_freqs[q1] = f
                            qubit_freqs[q2] = f
                            taus[(q1,q2)] = np.pi / (2 * 0.5 * np.sqrt(f*f) * Cqq)
                            edges_iswaps.append((q1,q2))
                        elif (g.name == 'unitary' and g.label == 'sqrtiswap'):
                            qubit_freqs[q1] = f
                            qubit_freqs[q2] = f
                            #tau_special[(q1,q2)] = 0.5
                            taus[(q1,q2)] = 0.5 * np.pi / (2 * 0.5 * np.sqrt(f*f) * Cqq)
                            edges_iswaps.append((q1,q2))
                        elif (g.name == 'cz' or g.label == 'cz'):
                            qubit_freqs[q1] = f
                            qubit_freqs[q2] = f - alphas[q2] # b/c match f1 with f2+alpha
                            taus[(q1,q2)] = np.pi / (np.sqrt(2) * 0.5 * np.sqrt(f*f) * Cqq) # f is interaction freq
                            edges_cphase.append((q1,q2))
                        else:
                            print("Gate %s(%s) not recognized. Supports iswap, sqrtiswap, cz." % (g.name, g.label))
                        t_2q[q1] += taus[(q1,q2)]
                        t_2q[q2] += taus[(q1,q2)]

                #if (scheduler == 'greedy'):
                #    edges, leftover, ind = greedy_reschedule(coupling, edges)
                #    for i in range(len(ind)):
                #        if (ind[i]):
                #            all_gates.append(curr_gates[i])
                #        else:
                #            left_gates.append(curr_gates[i])
                #else:
                #    for i in range(len(curr_gates)):
                #        all_gates.append(curr_gates[i])
                for i in range(len(curr_gates)):
                    all_gates.append(curr_gates[i])



                #print("all_gates:")
                #print(all_gates)
                #print("edges:")
                #print(edges)
                #print("leftover:")
                #print(leftover)
                #print("qubit_freqs:")
                #print(qubit_freqs)
                #err = compute_crosstalk(edges, coupling, qubit_freqs)
                print("qubit_freqs:")
                print(qubit_freqs)
                if barrier:
                    err, swap_err, leak_err = 0.0, 0.0, 0.0
                    gt = 0.0
                else:
                    if decomp == 'iswap':
                        #err, swap_err, leak_err = compute_crosstalk_iswaps(edges, coupling, qubit_freqs, tau_special)
                        #gt = get_iswap_time(edges, qubit_freqs, tau_special)
                        gt = get_max_time(gt, taus)
                        if scheduler == 'tiling':
                            err, swap_err, leak_err = compute_crosstalk_iswaps_gmon(edges_iswaps, coupling, qubit_freqs, taus, gt, tilings, tiling_idx)
                        else:
                            err, swap_err, leak_err = compute_crosstalk_iswaps(edges_iswaps, coupling, qubit_freqs, taus, gt)
                    elif decomp == 'cphase':
                        gt = get_max_time(gt, taus)
                        if scheduler == 'tiling':
                            err, swap_err, leak_err = compute_crosstalk_cphase_gmon(edges_iswaps, coupling, qubit_freqs, taus, gt, tilings, tiling_idx)
                        else:
                            err, swap_err, leak_err = compute_crosstalk_cphase(edges_iswaps, coupling, qubit_freqs, taus, gt)
                        #err, swap_err, leak_err = compute_crosstalk_cphase(edges, coupling, qubit_freqs, tau_special)
                        #gt = get_cphase_time(edges, qubit_freqs, tau_special)
                    elif decomp == 'flexible' or decomp == 'mixed':
                        gt = get_max_time(gt, taus)
                        if scheduler == 'tiling':
                            err, swap_err, leak_err = compute_crosstalk_flexible_gmon(edges_cphase, edges_iswaps, coupling, qubit_freqs, taus, gt, tilings, tiling_idx)
                            print("Err,swap_err,leak_err:", err, swap_err, leak_err)
                        else:
                            err, swap_err, leak_err = compute_crosstalk_flexible(edges_cphase, edges_iswaps, coupling, qubit_freqs, taus, gt)
                    else:
                        print("Decomposition method %s not recognized. Try iswap or cphase." % decomp)
                        gt = 0.0
                if gt > layer_time: layer_time = gt
                for qubit in range(num_q):
                    if active_list[qubit]:
                        t_act[qubit] += layer_time
                update_data(freqsdata, gatesdata, qubit_freqs, all_gates, num_q)
                success_rate *= 1 - err
                if verbose == 0:
                    print("Layer success: %12.10f (swap: %12.10f, leak: %12.10f)" % (1-err, swap_err, leak_err))
                    print("Layer time:", layer_time)
                # Reset the frequencies
                for (q1, q2) in edges:
                    qubit_freqs[q1] = park_freqs[q1]
                    qubit_freqs[q2] = park_freqs[q2]
                if not barrier:
                    tot_success += 1 - err
                    tot_cnt += 1
                worst_success = min(worst_success, 1 - err)
                total_time += layer_time
            #idx += 1
        write_circuit(freqsdata, gatesdata, outf)
    avg_success = tot_success / tot_cnt
    if decomp=='flexible' and verbose == 0:
        print("Total number of CNOT gates:", num_cnot)
        print("Total number of SWAP gates:", num_swap)
        print("Total number of 2-qubit gates that are decomposed with iSWAP:", num_iswap)
        if num_iswap+num_cphase > 0:
            print("Proportion:",num_iswap/(num_iswap+num_cphase))
        print("Total number of 2-qubit gates that are decomposed with CPHASE:", num_cphase)
        if num_iswap+num_cphase > 0:
            print("Proportion:",num_cphase/(num_iswap+num_cphase))
    return success_rate, avg_success, worst_success, idx, tot_cnt, total_time, max_colors, t_act, t_2q

def _tests():
    """A function to run tests on the module"""
    theta = [np.random.random() for _ in range(8)]
    circuit = get_uccsd_circuit('LiH', theta)
    print(circuit)

if __name__ == "__main__":
    _tests()
