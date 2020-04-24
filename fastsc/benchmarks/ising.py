from qiskit import QuantumCircuit, QuantumRegister

def get_ising_circuit(num_qubits, layers=5, x_theta=1.0, z_theta=1.0):
    """Linear ising model, a la https://arxiv.org/pdf/1511.03316.pdf."""
    qc = QuantumCircuit(num_qubits)
    for i in range(num_qubits):
        qc.h(i)
    for m in range(layers):
        for i in range(0, num_qubits - 1, 2):
            qc.cx(i, i + 1)
            qc.rz(z_theta*(1-m/layers), i + 1)
            qc.cx(i, i + 1)

        for i in range(1, num_qubits - 1, 2):
            qc.cx(i, i + 1)
            qc.rz(z_theta*(1-m/layers), i + 1)
            qc.cx(i, i + 1)

        for i in range(num_qubits):
            qc.rx(x_theta*m/layers, i)
    return qc
