
from qiskit import QuantumCircuit
import numpy as np

def get_qgan_circuit(num_qubits, layers=1):
    qc = QuantumCircuit(num_qubits)
    for _ in range(layers):
        for i in range(num_qubits):
            qc.ry(np.random.random() * 2 * np.pi, i)
        for i in range(0, num_qubits, 2):
            qc.cz(i, (i + 1) % num_qubits)
        for i in range(1, num_qubits - 1, 2):
            qc.cz(i, (i + 1) % num_qubits)
    return qc
