import numpy as np


def swap_channel(coupling, residual_factors, qubit_freqs, taus, Cqq=0.019):
    # return prob of rabi oscillation between 01 and 10 (i.e. iswap)
    # for each pair of coupled qubits
    # taus is list of hold durations: max swap happens at 2pi/g
    res = []
    if residual_factors==None:
        residual_factors = [1.0]*len(coupling)
    #print(taus)
    for (i, (q1,q2)) in enumerate(coupling):
        J = 0.5 * np.sqrt(qubit_freqs[q1] * qubit_freqs[q2]) * Cqq
        delta_omega = np.absolute(qubit_freqs[q1] - qubit_freqs[q2])
        if delta_omega <= 1e-5:
            residual_coupling = J
        else:
            residual_coupling = J**2 / (delta_omega * 4 * np.pi)
        tau = taus[(q1,q2)]
        res.append(residual_factors[i]*np.sin(residual_coupling * tau)**2)
    return res

def leak_channel(coupling, residual_factors, qubit_01freqs, qubit_12freqs, taus, Cqq = 0.019):
    # return prob of rabi oscillation between 11 and (20+02)/sqrt(2) (i.e. leakage)
    # for each pair of coupled qubits
    # taus is list of hold durations
    res = []
    if residual_factors==None:
        residual_factors = [1.0]*len(coupling)
    for (i, (q1,q2)) in enumerate(coupling):
        # 02
        f1 = qubit_01freqs[q1]
        f2 = qubit_12freqs[q2]
        J = 0.5 * np.sqrt(f1 * f2) * Cqq
        delta_omega = np.absolute(f1 - f2)
        if delta_omega <= 1e-5:
            residual_coupling = np.sqrt(2) * J
        else:
            residual_coupling = (np.sqrt(2)*J)**2 / (delta_omega * 4 * np.pi)
        tau = taus[(q1,q2)]
        epsilon_leak02 = residual_factors[i]*np.sin(residual_coupling * tau)**2
        #print("delta_omega: ", delta_omega, "epsilon_leak: ", epsilon_leak02)
        res.append(epsilon_leak02)
    return res

