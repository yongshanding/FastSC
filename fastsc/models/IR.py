"IR.py - IR (Intermediate Representation) of quantum circuit"

from qiskit import QuantumCircuit#, Qubit, Instruction
import copy

class Qbit(object):
    def __init__(self, q, omega01, omega12):
        self.register = q.register
        self.index = q.index
        self.idle_freq = [omega01, omega12]

class Inst(object):
    def __init__(self, ins, qargs, cargs, int_freq, gt):
        self.ins = ins
        self.name = ins.name
        self.label = ins.label
        self.qargs = qargs
        self.cargs = cargs
        self.int_freq = int_freq
        self.gate_time = gt
    def qasm(self):
        return self.ins.qasm()


class IR(object):
    def __init__(self, data=[], depth=0, qubits=[], width=0, coupling=[], alpha=-0.2):
        self.data = data
        self.depth = depth
        self.qubits = qubits
        self.width = width
        self.t_act = None
        self.t_2q = None
        self.num_1qg = 0
        self.depth_before = 0
        self.depth_after = 0
        self.total_time = None
        self.max_colors = None
        self.coupling = coupling
        self.alpha = alpha

    def append_layer(self, layer):
        self.data.append(layer)
        self.depth += 1

    def append_layer_from_insts(self, insts, coupling_factors=None):
        active = []
        gt = 0
        freqs = [[self.qubits[i].idle_freq[0], self.qubits[i].idle_freq[1]] for i in range(self.width)]
        if coupling_factors == None:
            coupling_factors = [1.0] * len(self.coupling)

        for inst in insts:
            int_freq = inst.int_freq
            for i in range(len(inst.qargs)):
                if inst.qargs[i] in active:
                    print("Warning: Qubit " + str(q.index) + " cannot be used twice in one time step.")
                else:
                   active.append(inst.qargs[i])
                if int_freq != None:
                    # two-qubit gate
                    freqs[inst.qargs[i].index][0] = int_freq[i] # omega01
                    freqs[inst.qargs[i].index][1] = int_freq[i] + self.alpha # omega12
            if int_freq == None:
                # single-qubit gate
                self.num_1qg += 1
            if inst.gate_time >= gt:
                gt = inst.gate_time
        self.data.append((insts, freqs, gt, coupling_factors))

