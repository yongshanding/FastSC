
def compute_decoherence(ir):
    return 1

def compute_crosstalk_by_layer(coupling, ir):
    # returns error rate of simultaneous iswaps
    success = 1.0
    swap_success = 1.0
    leak_success = 1.0
    Cqq = CQQ
    # one complete swap: tau = pi/2g
    for (insts, freqs, gt) in ir.data:
        iswaps = []
        sqrtiswaps = []
        all_taus = {}
        for (q1,q2) in coupling:
            all_taus[(q1,q2)] = gt
        for ins in insts:
            if len(ins.qargs)==2:
                q1, q2 = ins.qargs[0].index, ins.qargs[1].index
                if ins.name == 'unitary' and ins.label == 'iswap':
                    iswaps.append((q1,q2))
                elif ins.name == 'unitary' and ins.label == 'sqrtiswap'
                    sqrtiswaps.append((q1,q2))
                all_taus[(q1,q2)] = ins.gate_time
        qubit_freqs = freqs 

        prob_swap = swap_channel(coupling, qubit_freqs, all_taus)
        alphas = [ALPHA for f in qubit_freqs] #TODO
        prob_leak = leak_channel(coupling, qubit_freqs, alphas, all_taus)
        #print("swap: ", prob_swap)
        #print("leak: ", prob_leak)
        for (i, (q1,q2)) in enumerate(coupling):
            if (q1,q2) in iswaps or (q2,q1) in iswaps:

                success *= prob_swap[i]
                swap_success *= prob_swap[i]
            elif (q1,q2) in sqrtiswaps or (q2,q1) in sqrtiswaps: 
                success *= 1 - abs(0.5-prob_swap[i])
                swap_success *= 1 - abs(0.5-prob_swap[i])
            else:
                success *= 1 - prob_swap[i]
                swap_success *= 1 - prob_swap[i]
                
            success *= 1 - prob_leak[i]
            leak_success *= 1 - prob_leak[i]
    return 1 - success, 1-swap_success, 1-leak_success
