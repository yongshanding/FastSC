from fastsc.util import get_connectivity_graph, get_aug_line_graph, get_map_circuit, get_layer_circuits, get_nearest_neighbor_coupling_list, gen_tiling_pattern
import networkx as nx
import numpy as np
from ..models import Sycamore_device, IR, Qbit, Inst
from .util import relabel_coloring, get_qubits, decompose_layer, decompose_layer_flexible, reschedule_layer, get_max_time
from qiskit import QuantumCircuit

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


def google_like(device, circuit, scheduler, d, decomp, verbose, res_coup=0.0):
    freqsdata = []
    gatesdata = []
    width = device.side_length
    height = device.side_length
    num_q = width * height
    q_arr = get_qubits(circuit)
    ALPHA = device.alpha

    tilings = gen_tiling_pattern(device)
    G_connect = get_connectivity_graph(width, height)
    G_crosstalk = get_aug_line_graph(width, height, d)
    syc_device = Sycamore_device(device, num_q, res_coup)
    int_freqs = syc_device.int_freqs
    park_freqs = syc_device.park_freqs
    num_int = len(set(int_freqs.values()))

    def _initial_frequency():
        freqs = dict()
        for q in range(width*height):
            freqs[q] = park_freqs[q]
        return freqs

    t_act = np.zeros(num_q)
    t_2q = np.zeros(num_q)
    active_list = [False for i in range(num_q)]
    success_rate = 1.0
    # single_qb_err = 0.0015
    tot_success = 0.0
    tot_cnt = 0
    worst_success = 1.0
    max_colors = num_int # max number of colors used
    # Mapping
    coupling = device.coupling
    circ_mapped = get_map_circuit(circuit, coupling)

    Cqq = device.cqq
    park_freqs = _initial_frequency()
    alphas = [ALPHA for f in park_freqs]
    for i in range(num_q):
        q_arr[i].idle_freq = [park_freqs[i], park_freqs[i]+alphas[i]]
    ir = IR(qubits = q_arr, width = num_q, coupling = coupling, alpha = ALPHA)

    #circuit.draw(output='mpl')
    # Check scheduler
    if (scheduler == 'hybrid'):
        print("Hybrid scheduler to be implemented.")
        sys.exit(2)
    else:
        # scheduler == qiskit or greedy
        layers = get_layer_circuits(circ_mapped)
        qubit_freqs = park_freqs
        num_layers = len(layers)
        idx = 0
        pattern_offset = 0
        total_time = 0.0 # in ns
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
                decomp_layer= decompose_layer_flexible(layer_circuit, G_crosstalk, verbose)
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
                    if verbose == 0:
                        print("tiling_idx", tiling_idx)
                print(layer.qasm())
                insts = []
                edges = []
                #edges_cphase = [] # if (q1,q2) is in it, then (q2,q1) is also in it
                #edges_iswaps = [] # if (q1,q2) is in it, then (q2,q1) is also in it
                #edges = [e for e in leftover]
                #curr_gates = [e for e in left_gates]
                #curr_gates = []
                taus = {} # For storing couplings that needed sqrt iswaps
                gt = 0.0
                layer_time = 0.0 # ns
                barrier = False
                # single_qb_err = 0.0015
                for g, qargs, cargs in layer.data:
                    #print(g, qargs)
                    #print(g.qasm())
                    if g.name == "barrier": barrier = True
                    if g.name == "measure": continue
                    if len(qargs) == 1:
                        #all_gates.append((g.qasm(),(qargs[0].index, -1)))
                        active_list[qargs[0].index] = True
                        gt = device.gate_times[g.name]
                        if gt > layer_time: layer_time = gt
                        insts.append(Inst(g,qargs,cargs, None, gt))
                        # success_rate *= 1 - single_qb_err
                    if len(qargs) == 2:
                        q1, q2 = qargs[0].index, qargs[1].index
                        active_list[q1] = True
                        active_list[q2] = True
                        #edges.append((q1, q2))
                        #print((q1,q2))
                        #curr_gates.append((g.qasm(),(q1, q2)))
                        try:
                            f = int_freqs[(q1,q2)]
                        except:
                            f = int_freqs[(q2,q1)]
                        # f += get_flux_noise(device, f, sigma)
                        #print('freq:',f)
                        if (g.name == 'unitary' and g.label == 'iswap'):
                            f1 = f
                            f2 = f
                            taus[(q1,q2)] = np.pi / (2 * 0.5 * np.sqrt(f*f) * Cqq)
                        elif (g.name == 'unitary' and g.label == 'sqrtiswap'):
                            f1 = f
                            f2 = f
                            #tau_special[(q1,q2)] = 0.5
                            taus[(q1,q2)] = 0.5 * np.pi / (2 * 0.5 * np.sqrt(f*f) * Cqq)
                        elif (g.name == 'cz' or g.label == 'cz'):
                            f1 = f
                            f2 = f - alphas[q2] # b/c match f1 with f2+alpha
                            taus[(q1,q2)] = np.pi / (np.sqrt(2) * 0.5 * np.sqrt(f*f) * Cqq) # f is interaction freq
                        else:
                            print("Gate %s(%s) not recognized. Supports iswap, sqrtiswap, cz." % (g.name, g.label))
                        t_2q[q1] += taus[(q1,q2)]
                        t_2q[q2] += taus[(q1,q2)]
                        insts.append(Inst(g, qargs, cargs, [f1, f2], taus[(q1,q2)]))
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
                #for i in range(len(curr_gates)):
                #    all_gates.append(curr_gates[i])
                if not barrier:
                    coupling_factors = None
                    if (scheduler == 'tiling'):
                        tiling_active = tilings[tiling_idx]
                        coupling_factors = []
                        for (qa,qb) in coupling:
                            if (qa,qb) in tiling_active or (qb,qa) in tiling_active:
                                coupling_factors.append(1.0)
                            else:
                                # strength of residual coupling, 0~1
                                coupling_factors.append(syc_device.res_coupling) 
                    ir.append_layer_from_insts(insts, coupling_factors)
                    gt = get_max_time(gt, taus)
                    tot_cnt += 1
                    if gt > layer_time: layer_time = gt
                    total_time += layer_time
                    for qubit in range(num_q):
                        if active_list[qubit]:
                            t_act[qubit] += layer_time
            idx += 1
    ir.depth_before = idx
    ir.depth_after = tot_cnt
    ir.total_time = total_time
    ir.max_colors = max_colors
    ir.t_act = t_act
    ir.t_2q = t_2q
    return ir

