from fastsc.util import get_connectivity_graph, get_aug_line_graph, get_map_circuit, get_layer_circuits, get_nearest_neighbor_coupling_list
import networkx as nx
from ..models import IR, Qbit, Inst

def static_coloring(device, circuit, scheduler, d, decomp, verbose, uniform_freq):
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

    q_arr = get_qubits(circuit)

    if full_uni:
        int_coloring = {}
        for node in G_crosstalk.nodes:
            int_coloring[node] = 0
    else:
        int_coloring = nx.coloring.greedy_color(G_crosstalk)
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
        return omg #+ get_flux_noise(device, omg, sigma)#+np.random.normal(0,sigma)
    def _int_freq(c):
        omg = color_to_freq[str(c)]
        return omg #+ get_flux_noise(device, omg, sigma)#+np.random.normal(0,sigma)
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
    coupling = device.coupling
    circ_mapped = get_map_circuit(circuit, coupling)

    Cqq = CQQ
    park_freqs = _initial_frequency()
    alphas = [ALPHA for f in park_freqs]
    for i in range(num_q):
        qrr[i].idle_freq = [park_freqs[i], park_freqs[i]+alphas[i]]
    ir = IR(qubits = q_arr, width = num_q)

    #circuit.draw(output='mpl')
    # Check scheduler
    if (scheduler == 'hybrid'):
        print("Hybrid scheduler to be implemented.")
        sys.exit(2)
    else:
        # scheduler == qiskit or greedy
        layers = get_layer_circuits(circ_mapped)
        num_layers = len(layers)
        idx = 0
        # qubit_freqs = park_freqs
        total_time = 0.0 # in ns
        total_tcz = 0.0 # total time spent on CZ gates
        #while (idx < num_layers or len(leftover) > 0):
        for idx in range(num_layers):
            # all_gates = []
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
                decomp_layer = decompose_layer_flexible(layer_circuit, G_crosstalk, verbose)
            else:
                decomp_layer = decompose_layer(layer_circuit, decomp)
            if (scheduler == 'greedy'):
                resched_layer = reschedule_layer(decomp_layer, coupling, verbose)
            else:
                resched_layer = decomp_layer
            #print("Layers: %d %d" % (len(decomp_layer), len(resched_layer)))
            for layer in resched_layer:
                print(layer)
                insts = []
                #edges = []
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
                # single_qb_err_acc = 1.0
                for g, qargs, cargs in layer.data:
                    #print(g, qargs)
                    #print(g.qasm())
                    if g.name == "barrier": barrier = True
                    if g.name == "measure": continue
                    if len(qargs) == 1:
                        #all_gates.append((g.qasm(),(qargs[0].index, -1)))
                        active_list[qargs[0].index] = True
                        gt = GATETIMES[g.name]
                        if gt > layer_time: layer_time = gt
                        insts.append(g,qargs,cargs, None, gt)
                        # single_qb_err_acc *= 1 - single_qb_err
                    if len(qargs) == 2:
                        q1, q2 = qargs[0].index, qargs[1].index
                        active_list[q1] = True
                        active_list[q2] = True
                        #edges.append((q1, q2))
                        #print((q1,q2))
                        #curr_gates.append((g.qasm(),(q1, q2)))
                        try:
                            ic = int_coloring[(q1,q2)]
                        except:
                            ic = int_coloring[(q2,q1)]
                        f = _int_freq(ic)
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
                        insts.append(g, qargs, cargs, [f1, f2], taus[(q1,q2)])
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
                #print("qubit_freqs:")
                #print(qubit_freqs)
                if not barrier:
                    ir.append_layer_from_insts(insts)
                    gt = get_max_time(gt, taus)
                    tot_cnt += 1
                    if gt > layer_time: layer_time = gt
                    total_time += layer_time
                    for qubit in range(num_q):
                        if active_list[qubit]:
                            t_act[qubit] += layer_time
                # Reset the frequencies
                #for (q1, q2) in edges:
                #    qubit_freqs[q1] = _park_freq(park_coloring[q1])
                #    qubit_freqs[q2] = _park_freq(park_coloring[q2])
                    # tot_success += (1 - err)*single_qb_err_acc
                # worst_success = min(worst_success, 1 - err)
            idx += 1
    return ir, idx, tot_cnt, total_time, max_colors, t_act, t_2q
