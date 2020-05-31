from fastsc.util import get_connectivity_graph, get_aug_line_graph, get_map_circuit, get_layer_circuits, get_nearest_neighbor_coupling_list
import networkx as nx
import numpy as np
from fastsc.models import IR, Qbit, Inst
from .util import relabel_coloring, get_qubits, decompose_layer, decompose_layer_flexible, reschedule_layer, limit_colors, get_max_time

def color_dynamic(device, circuit, scheduler, d, decomp, lim_colors, verbose):
    freqsdata = []
    gatesdata = []
    width = device.side_length
    height = device.side_length
    num_q = width * height
    omega_max = device.omega_max
    delta_int = device.delta_int
    delta_ext= device.delta_ext
    delta_park = device.delta_park
    ALPHA = device.alpha

    G_connect = get_connectivity_graph(width, height)
    park_coloring = nx.coloring.greedy_color(G_connect)
    num_park = len(set(park_coloring.values()))
    G_crosstalk = get_aug_line_graph(width, height, d)
    coupling = device.coupling
    Cqq = device.cqq

    q_arr = get_qubits(circuit)

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
        omg = color_to_freq[str(-(c+1))]
        return omg #+ get_flux_noise(device, omg, sigma)
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
    park_freqs = _initial_frequency()
    alphas = [ALPHA for f in park_freqs]
    for i in range(num_q):
        q_arr[i].idle_freq = [park_freqs[i], park_freqs[i]+alphas[i]]
    ir = IR(qubits = q_arr, width = num_q, coupling = coupling, alpha = ALPHA)

    # Check scheduler
    if (scheduler == 'hybrid'):
        print("Hybrid scheduler to be implemented.")
        sys.exit(2)
    else:
        layers = get_layer_circuits(circ_mapped)
        num_layers = len(layers)
        idx = 0
        total_time = 0.0 # ns
        total_tcz = 0.0
        if verbose == 0:
            print("Num of layers:", num_layers)


        #while (idx < num_layers or len(leftover) > 0):
        for idx in range(num_layers):
            # all_gates = []
            layer_circuit = layers[idx]
            if verbose == 0:
                print(layer_circuit)
            print(idx, "-----------------")
            if decomp == 'flexible':
                if verbose == 0:
                    print("Decompose layer", idx)
                decomp_layer = decompose_layer_flexible(layer_circuit, G_crosstalk, verbose)
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
                insts = []
                # Pre-fill edges for constructing (undirected) xtalk graph
                #edges = [leftover[i//2] if i%2==0 else (leftover[i//2][1], leftover[i//2][0]) for i in range(2*len(leftover))]
                edges = []
                #edges_cphase = [] # if (q1,q2) is in it, then (q2,q1) is also in it
                #edges_iswaps = [] # if (q1,q2) is in it, then (q2,q1) is also in it
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
                    if num_int == 0:
                        idx += 1
                        continue
                    _add_int_color_map(color_to_freq, num_int)
                def _int_freq(c):
                    omg =  color_to_freq[str(c)]#+np.random.normal(0,sigma)
                    return omg #+ get_flux_noise(device, omg, sigma)
                #print(layer)
                #print("-----------------")
                # Refill edges and curr_gates
                #edges = [e for e in leftover]
                # edges = []
                #curr_gates = [e for e in left_gates]
                # curr_gates = []
                # single_qb_err = 0.0015
                # single_qb_err_acc = 1.0
                for g, qargs, cargs in layer.data:
                    if g.name == "barrier": barrier = True
                    if g.name == "measure": continue
                    #print(qargs)
                    #print(qargs[0].index)
                    if len(qargs) == 1: # single qubit gates
                        # all_gates.append((g.qasm(),(qargs[0].index, -1)))
                        active_list[qargs[0].index] = True
                        gt = device.gate_times[g.name]
                        if gt > layer_time: layer_time = gt
                        insts.append(Inst(g,qargs,cargs, None, gt))
                        # single_qb_err_acc *= 1 - single_qb_err
                    elif len(qargs) == 2:
                        q1, q2 = qargs[0].index, qargs[1].index
                        active_list[q1] = True
                        active_list[q2] = True
                        #edges.append((q1, q2))
                        #curr_gates.append((g.qasm(),(q1, q2)))
                        try:
                            f = _int_freq(int_coloring[(q1, q2)])
                        except:
                            f = _int_freq(int_coloring[(q2, q1)])
                        if (g.name == 'unitary' and g.label == 'iswap'):
                            f1 = f
                            f2 = f
                            taus[(q1,q2)] = np.pi / (2 * 0.5 * np.sqrt(f*f) * Cqq)
                        elif (g.name == 'unitary' and g.label == 'sqrtiswap'):
                            f1 = f
                            f2 = f
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
                # success_rate *= single_qb_err_acc
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
                    ir.append_layer_from_insts(insts)
                    gt = get_max_time(gt, taus)
                    tot_cnt += 1
                    if gt > layer_time: layer_time = gt
                    total_time += layer_time
                    for qubit in range(num_q):
                        if active_list[qubit]:
                            t_act[qubit] += layer_time
            idx += 1
    ir.t_act = t_act
    ir.t_2q = t_2q
    ir.depth_before = idx
    ir.depth_after = tot_cnt
    ir.total_time = total_time
    ir.max_colors = max_colors
    return ir 
