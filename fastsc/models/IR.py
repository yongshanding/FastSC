"IR.py - IR (Intermediate Representation) of quantum circuit"

from qiskit import QuantumCircuit#, Qubit, Instruction

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
    def __init__(self, data=[], depth=0, qubits=[], width=0):
        self.data = data
        self.depth = depth
        self.qubits = qubits
        self.width = width

    def append_layer(self, layer):
        self.data.append(layer)
        self.depth += 1

    def append_layer_from_insts(self, insts, coupling_factors=[]):
        active = []
        gt = 0
        freqs = [qubits[i].idle_freq for i in range(self.width)]

        for ins in insts:
            int_freq = ins.int_freq
            for i in range(len(inst.qargs)):
                if inst.qargs[i] in active:
                    print("Warning: Qubit " + str(q.index) + " cannot be used twice in one time step.")
                else:
                   active.append(inst.qargs[i])
                if int_freq is not None:
                    freqs[inst.qargs[i].index] = int_freq[i]
            if ins.gate_time >= gt:
                gt = ins.gate_time
        self.data.append((insts, freqs, gt, coupling_factors))
