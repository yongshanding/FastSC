"""
util.py -  Utility functions
"""

"""
Input: A dictionary with (key, value) = (edge, int color)
Return: A dictionary with (key, value) = (edge, int color), such that if a<b
then more nodes are colored with a than b
"""
from fastsc.models import Qbit
from qiskit import QuantumCircuit 
from fastsc.util import get_layer_circuits
import numpy as np
import networkx as nx
from qiskit.quantum_info.operators import Operator

def iswap():
    return Operator([[1., 0., 0., 0.],[0., 0., -1.j, 0.],[0., -1.j, 0., 0.], [0., 0., 0., 1.]])

def sqrtiswap():
    return Operator([[1., 0., 0., 0.],[0., 1./np.sqrt(2), -1.j/np.sqrt(2), 0.],[0., -1.j/np.sqrt(2), 1./np.sqrt(2), 0.], [0., 0., 0., 1.]])

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

"""
Returns an array of Qbit objects
"""
def get_qubits(circ):
    arr = []
    for i in range(len(circ.qubits)):
        qq = Qbit(circ.qubits[i], i, None)
        arr.append(qq)
    return arr

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
    return get_layer_circuits(new_circuit)

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

def get_max_time(gt, taus):
    tmax = gt
    for t in taus.values():
        if t > tmax: tmax = t
    return tmax

