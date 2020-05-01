from fastsc.util import get_connectivity_graph, get_aug_line_graph, get_map_circuit, get_layer_circuits, get_nearest_neighbor_coupling_list
import networkx as nx
from ..models import IR, Qbit, Inst

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


def google_like(device, circuit, scheduler, d, decomp, verbose):
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
        # write_circuit(freqsdata, gatesdata, outf)
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
