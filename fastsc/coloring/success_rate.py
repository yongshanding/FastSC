from ..models import swap_channel, leak_channel, get_flux_noise
import numpy as np
import math, random

def compute_decoherence(device, t_act, t_2q):
    # calculate decoherence
    decoh = 1.
    random.seed(60615)
    np.random.seed(60615)
    for i in range(device.qubits):
        # relaxation times for qubit i
        T1_mean = np.random.uniform(20000.,35000.)
        T1 = np.random.normal(T1_mean,10000.)
        while T1 <= 0:
            T1 = np.random.normal(T1_mean,10000.)
        T2 = T1*min(2.,np.random.normal(1.8,0.2))
        T1_tilde = T1*np.random.normal(0.7,0.1)
        T2_tilde = T2*np.random.normal(0.5,0.1)
        t1 = t_act[i]-t_2q[i]
        t2 = t_2q[i] # time spent on 2q gates
        decoh *= math.exp(-1.*t1/T1-t1/T2)
        decoh *= math.exp(-1.*t2/T1_tilde-t2/T2_tilde)
    return decoh

# TODO: add single qubit errors: single_qb_err = 0.0015
# TODO: provide residual_factors for swap_channel and leak_channel
def compute_crosstalk_by_layer(device, ir):
    # returns error rate of simultaneous iswaps
    success = 1.0
    swap_success = 1.0
    leak_success = 1.0
    Cqq = device.cqq
    ALPHA = device.alpha
    coupling = device.coupling
    # one complete swap: tau = pi/2g
    for (insts, freqs, gt) in ir.data:
        # Iterate over layer
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
                elif ins.name == 'unitary' and ins.label == 'sqrtiswap':
                    sqrtiswaps.append((q1,q2))
                all_taus[(q1,q2)] = ins.gate_time
        qubit_freqs = [f + get_flux_noise(f, device) for f in freqs]

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
