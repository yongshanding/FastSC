
from qiskit import QuantumCircuit
import numpy as np
from qiskit.quantum_info.operators import Operator

def sqrtx():
    """Returns the single qubit gate Sqrt(X)"""
    return Operator([[(1.+1.j)/2,(1.-1.j)/2],[(1.-1.j)/2,(1.+1.j)/2]])

def sqrty():
    """Returns the single qubit gate Sqrt(Y)"""
    return Operator([[(1.+1.j)/2,(-1-1.j)/2],[(1.+1.j)/2,(1.+1.j)/2]])

def sqrtw():
    """Returns the single qubit gate Sqrt(W)"""
    return Operator([[(1.+1.j)/2,-1.j/np.sqrt(2)],[1./np.sqrt(2),(1.+1.j)/2]])

def get_parallel_cnot(numQ, dep):
    assert numQ == 16, 'only doing numQ = 16'
    qc = QuantumCircuit(num_qubits)

    control1 = 5
    target1 = 6
    control2 = 9
    target2 = 10

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
    cnot_circ.cnot(control1, target1)
    cnot_circ.cnot(control2, target2)

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
        cnot_circ.cnot(control1, target1)
        cnot_circ.cnot(control2, target2)
    return qc

def get_parallel_cnot_barriers(numQ, dep):
    assert numQ == 16, 'only doing numQ = 16'
    qc = QuantumCircuit(num_qubits)

    control1 = 5
    target1 = 6
    control2 = 9
    target2 = 10

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
    cnot_circ.cnot(control1, target1)
    cnot_circ.cnot(control2, target2)
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
        cnot_circ.cnot(control1, target1)
        cnot_circ.cnot(control2, target2)
        qc.barrier()
    return qc
