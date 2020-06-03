
from qiskit import QuantumCircuit
from qiskit.quantum_info.operators import Operator
import numpy as np
import networkx as nx

def iswap(params):
    return Operator([[1., 0., 0., 0.],[0., 0., -1.j, 0.],[0., -1.j, 0., 0.], [0., 0., 0., 1.]])

def sqrtx():
    """Returns the single qubit gate Sqrt(X)"""
    return Operator([[(1.+1.j)/2,(1.-1.j)/2],[(1.-1.j)/2,(1.+1.j)/2]])

def sqrty():
    """Returns the single qubit gate Sqrt(Y)"""
    return Operator([[(1.+1.j)/2,(-1-1.j)/2],[(1.+1.j)/2,(1.+1.j)/2]])

def sqrtw():
    """Returns the single qubit gate Sqrt(W)"""
    return Operator([[(1.+1.j)/2,-1.j/np.sqrt(2)],[1./np.sqrt(2),(1.+1.j)/2]])

def get_xeb_circuit(num_qubits, depth):
    """Cross-entropy benchmark circuit based on Google supremacy paper"""
    assert int(np.sqrt(num_qubits)) ** 2 == num_qubits, 'only doing square grids for now'
    qc = QuantumCircuit(num_qubits)
    dim = int(np.sqrt(num_qubits))
    G = nx.convert_node_labels_to_integers(nx.grid_2d_graph(dim, dim))
    # partitioning the couplings into 4 sets
    setA = []
    setB = []
    setC = []
    setD = []
    for (u,v) in G.edges:
        assert u < v, 'assuming that every nx edge (u,v) satisfies u < v'
        # if the edge is horizontal
        # row: int(u/dim), col: u%dim
        if v-u == 1:
            if (int(u/dim)+u%dim)%2 == 0:
                setD.append((u,v))
            else:
                setC.append((u,v))
        else: # if the edge is vertical
            if (int(u/dim)+u%dim)%2 == 0:
                setA.append((u,v))
            else:
                setB.append((u,v))

    # testing
    print("Size of ABCD:",len(setA), len(setB), len(setC), len(setD))

    # add gates to the circuits
    # cycle 1 - assign each qubit one of the 3 single-qubit gates, uniformly at random
    cur_rand = np.random.random(num_qubits)
    for i in range(num_qubits):
        if cur_rand[i] < 1/3:
            qc.unitary(sqrtx(),[i],label='sqrtx')
            cur_rand[i] = 0 # labels: 0 - x, 1 - y, 2 - w
        elif cur_rand[i] < 2/3:
            qc.unitary(sqrty(),[i],label='sqrty')
            cur_rand[i] = 1
        else:
            qc.unitary(sqrtw(),[i],label='sqrtw')
            cur_rand[i] = 2
    # 2-qubit gates in cycle 1
    for (u,v) in setA:
        qc.cz(u,v)
    for i in range(1,depth):
        # make sure that gates do not repeat sequentially
        prev_rand = cur_rand
        cur_rand = np.random.random(num_qubits)
        # add single qubit gates in cycle i
        for j in range(num_qubits):
            if prev_rand[j] == 0: # if the last gate is sqrtx
                if cur_rand[j] < 0.5:
                    qc.unitary(sqrty(),[j],label='sqrty')
                    cur_rand[j] = 1
                else:
                    qc.unitary(sqrtw(),[j],label='sqrtw')
                    cur_rand[j] = 2
            elif prev_rand[j] == 1: # if the last gate is sqrty
                if cur_rand[j] < 0.5:
                    qc.unitary(sqrtx(),[j],label='sqrtx')
                    cur_rand[j] = 0
                else:
                    qc.unitary(sqrtw(),[j],label='sqrtw')
                    cur_rand[j] = 2
            else: # if the last gate is sqrtw
                if cur_rand[j] < 0.5:
                    qc.unitary(sqrty(),[j],label='sqrty')
                    cur_rand[j] = 1
                else:
                    qc.unitary(sqrtx(),[j],label='sqrtx')
                    cur_rand[j] = 0
        # add 2-qubit gates, repeating the sequence "ABCDCDAB"
        if i%8 == 0 or i%8 == 6: # A
            for (u,v) in setA:
                qc.cz(u,v)
        elif i%8 == 1 or i%8 == 7: # B
            for (u,v) in setB:
                qc.cz(u,v)
        elif i%8 == 2 or i%8 == 4: # C
            for (u,v) in setC:
                qc.cz(u,v)
        else: # D
            for (u,v) in setD:
                qc.cz(u,v)
    return qc

def get_xeb_iswap_barriers_circuit(num_qubits, depth):
    """Cross-entropy benchmark circuit based on Google supremacy paper, with iSWAPs instead of CZs"""
    assert int(np.sqrt(num_qubits)) ** 2 == num_qubits, 'only doing square grids for now'
    qc = QuantumCircuit(num_qubits)
    dim = int(np.sqrt(num_qubits))
    G = nx.convert_node_labels_to_integers(nx.grid_2d_graph(dim, dim))
    # partitioning the couplings into 4 sets
    setA = []
    setB = []
    setC = []
    setD = []
    for (u,v) in G.edges:
        assert u < v, 'assuming that every nx edge (u,v) satisfies u < v'
        # if the edge is horizontal
        # row: int(u/dim), col: u%dim
        if v-u == 1:
            if (int(u/dim)+u%dim)%2 == 0:
                setD.append((u,v))
            else:
                setC.append((u,v))
        else: # if the edge is vertical
            if (int(u/dim)+u%dim)%2 == 0:
                setA.append((u,v))
            else:
                setB.append((u,v))

    # testing
    print("Size of ABCD:",len(setA), len(setB), len(setC), len(setD))

    # add gates to the circuits
    # cycle 1 - assign each qubit one of the 3 single-qubit gates, uniformly at random
    cur_rand = np.random.random(num_qubits)
    for i in range(num_qubits):
        if cur_rand[i] < 1/3:
            qc.unitary(sqrtx(),[i],label='sqrtx')
            cur_rand[i] = 0 # labels: 0 - x, 1 - y, 2 - w
        elif cur_rand[i] < 2/3:
            qc.unitary(sqrty(),[i],label='sqrty')
            cur_rand[i] = 1
        else:
            qc.unitary(sqrtw(),[i],label='sqrtw')
            cur_rand[i] = 2
    qc.barrier()
    # 2-qubit gates in cycle 1
    for (u,v) in setA:
        qc.unitary(iswap([u, v]), [u, v], label='iswap')
    qc.barrier()
    for i in range(1,depth):
        # make sure that gates do not repeat sequentially
        prev_rand = cur_rand
        cur_rand = np.random.random(num_qubits)
        # add single qubit gates in cycle i
        for j in range(num_qubits):
            if prev_rand[j] == 0: # if the last gate is sqrtx
                if cur_rand[j] < 0.5:
                    qc.unitary(sqrty(),[j],label='sqrty')
                    cur_rand[j] = 1
                else:
                    qc.unitary(sqrtw(),[j],label='sqrtw')
                    cur_rand[j] = 2
            elif prev_rand[j] == 1: # if the last gate is sqrty
                if cur_rand[j] < 0.5:
                    qc.unitary(sqrtx(),[j],label='sqrtx')
                    cur_rand[j] = 0
                else:
                    qc.unitary(sqrtw(),[j],label='sqrtw')
                    cur_rand[j] = 2
            else: # if the last gate is sqrtw
                if cur_rand[j] < 0.5:
                    qc.unitary(sqrty(),[j],label='sqrty')
                    cur_rand[j] = 1
                else:
                    qc.unitary(sqrtx(),[j],label='sqrtx')
                    cur_rand[j] = 0
        qc.barrier()
        # add 2-qubit gates, repeating the sequence "ABCDCDAB"
        if i%8 == 0 or i%8 == 6: # A
            for (u,v) in setA:
                qc.unitary(iswap([u, v]), [u, v], label='iswap')
        elif i%8 == 1 or i%8 == 7: # B
            for (u,v) in setB:
                qc.unitary(iswap([u, v]), [u, v], label='iswap')
        elif i%8 == 2 or i%8 == 4: # C
            for (u,v) in setC:
                qc.unitary(iswap([u, v]), [u, v], label='iswap')
        else: # D
            for (u,v) in setD:
                qc.unitary(iswap([u, v]), [u, v], label='iswap')
        qc.barrier()
    return qc

def get_xeb_iswap_circuit(num_qubits, depth):
    """Cross-entropy benchmark circuit based on Google supremacy paper, with iSWAPs instead of CZs"""
    assert int(np.sqrt(num_qubits)) ** 2 == num_qubits, 'only doing square grids for now'
    qc = QuantumCircuit(num_qubits)
    dim = int(np.sqrt(num_qubits))
    G = nx.convert_node_labels_to_integers(nx.grid_2d_graph(dim, dim))
    # partitioning the couplings into 4 sets
    setA = []
    setB = []
    setC = []
    setD = []
    for (u,v) in G.edges:
        assert u < v, 'assuming that every nx edge (u,v) satisfies u < v'
        # if the edge is horizontal
        # row: int(u/dim), col: u%dim
        if v-u == 1:
            if (int(u/dim)+u%dim)%2 == 0:
                setD.append((u,v))
            else:
                setC.append((u,v))
        else: # if the edge is vertical
            if (int(u/dim)+u%dim)%2 == 0:
                setA.append((u,v))
            else:
                setB.append((u,v))

    # testing
    print("Size of ABCD:",len(setA), len(setB), len(setC), len(setD))

    # add gates to the circuits
    # cycle 1 - assign each qubit one of the 3 single-qubit gates, uniformly at random
    cur_rand = np.random.random(num_qubits)
    for i in range(num_qubits):
        if cur_rand[i] < 1/3:
            qc.unitary(sqrtx(),[i],label='sqrtx')
            cur_rand[i] = 0 # labels: 0 - x, 1 - y, 2 - w
        elif cur_rand[i] < 2/3:
            qc.unitary(sqrty(),[i],label='sqrty')
            cur_rand[i] = 1
        else:
            qc.unitary(sqrtw(),[i],label='sqrtw')
            cur_rand[i] = 2
    # 2-qubit gates in cycle 1
    for (u,v) in setA:
        qc.unitary(iswap([u, v]), [u, v], label='iswap')
    for i in range(1,depth):
        # make sure that gates do not repeat sequentially
        prev_rand = cur_rand
        cur_rand = np.random.random(num_qubits)
        # add single qubit gates in cycle i
        for j in range(num_qubits):
            if prev_rand[j] == 0: # if the last gate is sqrtx
                if cur_rand[j] < 0.5:
                    qc.unitary(sqrty(),[j],label='sqrty')
                    cur_rand[j] = 1
                else:
                    qc.unitary(sqrtw(),[j],label='sqrtw')
                    cur_rand[j] = 2
            elif prev_rand[j] == 1: # if the last gate is sqrty
                if cur_rand[j] < 0.5:
                    qc.unitary(sqrtx(),[j],label='sqrtx')
                    cur_rand[j] = 0
                else:
                    qc.unitary(sqrtw(),[j],label='sqrtw')
                    cur_rand[j] = 2
            else: # if the last gate is sqrtw
                if cur_rand[j] < 0.5:
                    qc.unitary(sqrty(),[j],label='sqrty')
                    cur_rand[j] = 1
                else:
                    qc.unitary(sqrtx(),[j],label='sqrtx')
                    cur_rand[j] = 0
        # add 2-qubit gates, repeating the sequence "ABCDCDAB"
        if i%8 == 0 or i%8 == 6: # A
            for (u,v) in setA:
                qc.unitary(iswap([u, v]), [u, v], label='iswap')
        elif i%8 == 1 or i%8 == 7: # B
            for (u,v) in setB:
                qc.unitary(iswap([u, v]), [u, v], label='iswap')
        elif i%8 == 2 or i%8 == 4: # C
            for (u,v) in setC:
                qc.unitary(iswap([u, v]), [u, v], label='iswap')
        else: # D
            for (u,v) in setD:
                qc.unitary(iswap([u, v]), [u, v], label='iswap')
    return qc


