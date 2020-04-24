from fastsc.util import get_connectivity_graph, get_aug_line_graph, get_map_circuit, get_layer_circuits, get_nearest_neighbor_coupling_list
import networkx as nx


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

    return 
