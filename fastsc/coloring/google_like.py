from fastsc.util import get_connectivity_graph, get_aug_line_graph, get_map_circuit, get_layer_circuits, get_nearest_neighbor_coupling_list
import networkx as nx

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
    G_crosstalk = get_aug_line_graph(width, height, d)
    int_freqs = {}
    # TODO: get the parking and int frequencies from device.py
    num_int = len(set(int_freqs.values()))
    # num_park = len(set(park_freqs.values()))


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

    return
