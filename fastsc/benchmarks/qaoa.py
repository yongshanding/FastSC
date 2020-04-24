from qiskit import QuantumCircuit
import numpy as np
import networkx as nx

def get_qaoa_circuit(num_qubits, probability, p):  # maxcut on an N-node ErdosRenyi graph
    G = nx.erdos_renyi_graph(num_qubits, probability)
    parameters = [np.random.random() for _ in range(2 * p)] # random values for parameters

    qc = QuantumCircuit(num_qubits)
    for i in range(num_qubits):
        qc.h(i)
    for _ in range(p):
        edges = list(G.edges())[:]
        beta, gamma = parameters.pop(0), parameters.pop(0)

        remaining_edges = []
        occupied_qubits = []
        while len(edges) != 0:
            for edge in edges:  # Greedily parallelize edges
                (i, j) = edge
                if i in occupied_qubits or j in occupied_qubits:
                    remaining_edges.append(edge)
                else:
                    qc.cx(i, j)
                    qc.rz(gamma, j)
                    qc.cx(i, j)
                    occupied_qubits.extend([i, j])
            edges = remaining_edges
            remaining_edges, occupied_qubits = [], []
        for i in range(num_qubits):
            qc.rx(beta, i)

    return qc
