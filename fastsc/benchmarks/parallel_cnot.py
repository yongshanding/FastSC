
from qiskit import QuantumCircuit
import numpy as np

def get_parallel_cnot(nQubits):

    cnot_circ = QuantumCircuit(nQubits)
    side = int(np.sqrt(nQubits))
    a = int(nQubits / 2 - 1)
    cnot_circ.cnot(a, a+1)
    cnot_circ.cnot(a+side, a+side+1)

    return cnot_circ
