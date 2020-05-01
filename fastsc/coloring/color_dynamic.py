from fastsc.util import get_connectivity_graph, get_aug_line_graph, get_map_circuit, get_layer_circuits, get_nearest_neighbor_coupling_list
import networkx as nx
from ..models import IR, Qbit, Inst

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
