import qiskit

def get_bv_circuit(nQubits, hiddenString):
    qr = qiskit.QuantumRegister(nQubits, name='q')
    # for recording the measurement on qr
    cr = qiskit.ClassicalRegister(nQubits-1)

    bvCircuit = qiskit.QuantumCircuit(qr, cr)

    # Apply Hadamard gates to the first
    # (nQubits - 1) before querying the oracle
    for i in range(nQubits-1):
        bvCircuit.h(qr[i])

    # Apply 1 and Hadamard gate to the last qubit
    # for storing the oracle's answer
    bvCircuit.x(qr[nQubits-1])
    bvCircuit.h(qr[nQubits-1])

    # Apply the inner-product oracle
    hiddenString = hiddenString[::-1]
    for i in range(len(hiddenString)):
        if hiddenString[i] == "1":
            bvCircuit.cx(qr[i], qr[nQubits-1])
    hiddenString = hiddenString[::-1]

    # Apply Hadamard gates after querying the oracle
    for i in range(nQubits-1):
        bvCircuit.h(qr[i])

    # Measurement
    for i in range(nQubits-1):
        bvCircuit.measure(qr[i], cr[i])

    return bvCircuit
