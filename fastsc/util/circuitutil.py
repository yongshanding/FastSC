"""
circuitutil.py - A module for extending Qiskit circuit functionality.
"""

import numpy as np
import pickle
#from qiskit import Aer, BasicAer, QuantumCircuit, QuantumRegister, execute
#from qiskit.extensions.standard import *
#from qiskit.mapper import CouplingMap, swap_mapper
#from qiskit.converters import circuit_to_dag, dag_to_circuit
#from qiskit.transpiler import PassManager, transpile
#from qiskit.transpiler.passes import (BasicSwap, CXCancellation, HCancellation)

#from qiskit import Aer, BasicAer, 
from qiskit import QuantumCircuit, QuantumRegister, execute
from qiskit.circuit.library.standard_gates import *
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.transpiler import PassManager, CouplingMap
from qiskit.compiler import transpile
from qiskit.transpiler.passes import (BasicSwap, CXCancellation)

import re, math
import networkx as nx


### CONSTANTS ###

# NOTICE: GATE_TO_PULSE_TIME is kept here for dependency reasons, but all
# future references to this dict and any other experimental constants
# should be kept in fqc/data/data.py.
# See Gate_Times.ipnyb for determination of these pulse times
GATE_TO_PULSE_TIME = {'h': 1.4, 'cx': 3.8, 'rz': 0.4, 'rx': 2.5, 'x': 2.5, 'swap': 7.4, 'id': 0.0}

#unitary_backend = BasicAer.get_backend('unitary_simulator')
#state_backend = Aer.get_backend('statevector_simulator')

### FUNCTIONS ###

#def get_unitary(circuit):
#    """Given a qiskit circuit, produce a unitary matrix to represent it.
#    Args:
#    circuit :: qiskit.QuantumCircuit - an arbitrary quantum circuit
#    Returns:
#    matrix :: np.matrix - the unitary representing the circuit
#    """
#    job = execute(circuit, unitary_backend)
#    unitary = job.result().get_unitary(circuit, decimals=10)
#    return np.matrix(unitary)

def get_map_circuit(circuit, coupling_list=None):

    #merge_rotation_gates(circuit)

    coupling_map = None if coupling_list is None else CouplingMap(couplinglist=coupling_list)
    bs = BasicSwap(coupling_map=coupling_map)
    pass_manager = PassManager(bs)
    # Some CNOT identities are interleaved between others,
    # for this reason a second pass is required. More passes
    # may be required for other circuits.
    pass_manager.append(CXCancellation())
    #if coupling_map is not None:
    #    pass_manager.append(BasicSwap(coupling_map))

    optimized_circuit = transpile(circuit, #backend=state_backend,
                                  #coupling_map=coupling_map,
                                  pass_manager=pass_manager)

    return optimized_circuit

def get_layer_circuits(circuit):
    """Returns a list of circuits representing each simultaneous layer of circuit."""
    dagcircuit = circuit_to_dag(circuit)
    layer_circuits = []
    for layer in dagcircuit.layers():
        layer_circuits.append(dag_to_circuit(layer['graph']))
    return layer_circuits


#######################################
### Connectivity Graph Construction ###
#######################################

def get_connectivity_graph(qubits, topology='grid', param=None):
    attempt = 0
    while attempt < 10:
        if topology == 'grid':
            # assume square grid
            side = int(np.sqrt(qubits))
            G = nx.grid_2d_graph(side, side)
        elif topology == 'erdosrenyi':
            if param == None:
                print("Erdos Renyi graph needs parameter p.")
            G = nx.fast_gnp_random_graph(qubits, param)
        elif topology == 'turan':
            if param == None:
                print("Turan graph needs parameter r.")
            G = nx.turan_graph(qubits, param)
        elif topology == 'regular':
            if param == None:
                print("d-regular graph needs parameter d.")
            G = nx.random_regular_graph(param, qubits)
        elif topology == 'cycle':
            G = nx.cycle_graph(qubits)
        elif topology == 'wheel':
            G = nx.wheel_graph(qubits)
        elif topology == 'complete':
            G = nx.complete_graph(qubits)
        elif topology == 'hexagonal':
            # assume square hexagonal grid, node = 2(m+1)**2-2
            side = int(np.sqrt((qubits+2)/2))-1
            G = nx.hexagonal_lattice_graph(side, side)
        elif topology == 'path':
            G = nx.path_graph(qubits)
        elif topology == 'ibm_falcon':
            # https://www.ibm.com/blogs/research/2020/07/qv32-performance/
            # 27 qubits
            G = nx.empty_graph(27)
            G.name = "ibm_falcon"
            G.add_edges_from([(0,1),(1,2),(2,3),(3,4),(3,5),(5,6),
                              (6,7),(7,8),(8,9),(8,10),(10,11),
                              (1,12),(6,13),(11,14),
                              (12,15),(15,16),(16,17),(17,18),(17,19),
                              (19,20),(13,20),(20,21),(21,22),(22,23),
                              (22,24),(24,25),(14,25),(25,26)])
        elif topology == 'ibm_penguin':
            # 20 qubits
            G = nx.empty_graph(20)
            G.name = "ibm_penguin"
            G.add_edges_from([(0,1),(1,2),(2,3),(3,4),
                              (0,5),(5,6),(6,7),(7,8),(8,9),(4,9),
                              (5,10),(10,11),(11,12),(7,12),(12,13),(13,14),(9,14),
                              (10,15),(15,16),(16,17),(17,18),(18,19),(14,19)])
        elif topology == '1express': # path with express channels
            G = nx.convert_node_labels_to_integers(nx.path_graph(qubits))
            G.add_edges_from([(s,s+param) for s in range(0,qubits-param,param//2)])
        elif topology == '2express': # grid with express channels
            side = int(np.sqrt(qubits))
            G = nx.convert_node_labels_to_integers(nx.grid_2d_graph(side, side))
            G.add_edges_from([(s,s+param) for x in range(side) for s in range(x*side,x*side+side-param,param//2)]) # rows
            G.add_edges_from([(s,s+param*side) for y in range(side) for s in range(y,y+side*(side-param),param//2*side)]) # cols
        else:
            print("Topology %s not recognized; use empty graph instead." % topology)
            G = nx.empty_graph(qubits)
        if nx.is_connected(G) or nx.is_empty(G):
            break
        else:
            attempt += 1

    return nx.convert_node_labels_to_integers(G)



def get_grid_coupling_list(width, height, directed=True):
    """Returns a coupling list for nearest neighbor (rectilinear grid) architecture.
    Qubits are numbered in row-major order with 0 at the top left and
    (width*height - 1) at the bottom right.
    If directed is True, the coupling list includes both  [a, b] and [b, a] for each edge.
    """
    coupling_list = []

    def _qubit_number(row, col):
        return row * width + col

    # horizontal edges
    for row in range(height):
        for col in range(width - 1):
            coupling_list.append((_qubit_number(row, col), _qubit_number(row, col + 1)))
            if directed:
                coupling_list.append((_qubit_number(row, col + 1), _qubit_number(row, col)))

    # vertical edges
    for col in range(width):
        for row in range(height - 1):
            coupling_list.append((_qubit_number(row, col), _qubit_number(row + 1, col)))
            if directed:
                coupling_list.append((_qubit_number(row + 1, col), _qubit_number(row, col)))

    return coupling_list


####################################
### Crosstalk Graph Construction ###
####################################

def get_crosstalk_graph(g_conn, topology='grid', d=1):
    if topology == 'grid':
        nq = len(g_conn)
        side_length = int(np.sqrt(nq))
        return get_aug_line_graph(side_length, side_length, d)
    else:
        if d != 1: print("Distance > 1 not yet supported for non-grid.")
        out_graph = nx.line_graph(g_conn)
        augmenting = []
        vertices = list(out_graph.nodes())
        marked = {}
        for from_node in vertices:
            marked[from_node] = True
            for n in out_graph.neighbors(from_node):
                for to_node in out_graph.neighbors(n):
                    if not(to_node in marked): # undirected and prevent selfloop
                        augmenting.append((from_node, to_node))
        out_graph.add_edges_from(augmenting)
        return out_graph
        

def get_aug_line_graph(width, height, d):
    
    def _coordinate_distance(point_a, point_b):
        x_a, y_a = int(point_a) // height, int(point_a) % height
        x_b, y_b = int(point_b) // height, int(point_b) % height
        return math.fabs(x_a - x_b) + math.fabs(y_a - y_b)

    out_graph = nx.line_graph(get_connectivity_graph(width*height, 'grid'))
    augmenting_depth_one = []
    
    vertices = list(out_graph.nodes())
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            point1, point2 = vertices[i]
            point3, point4 = vertices[j]
            if _coordinate_distance(point1, point3) <= d or _coordinate_distance(point1, point4) <= d or _coordinate_distance(point2, point3) <= d or _coordinate_distance(point2, point4) <= d:
                augmenting_depth_one.append((vertices[i], vertices[j]))
    
    out_graph.add_edges_from(augmenting_depth_one)
    return out_graph



def gen_tiling_pattern(device):
    # Generate coupler activation pattern, assume 2D grid of qubits
    width = device.side_length
    height = device.side_length
    num_q = width * height
    patternA = []
    patternB = []
    patternC = []
    patternD = []

    def _qubit_number(row, col):
        return row * width + col

    # horizontal edges
    for row in range(height):
        for col in range(width - 1):
            if (col+row)%2 == 1:
                patternC.append((_qubit_number(row, col), _qubit_number(row, col + 1)))
            else:
                patternD.append((_qubit_number(row, col), _qubit_number(row, col + 1)))
            #coupling_list.append((_qubit_number(row, col), _qubit_number(row, col + 1)))
            #if directed:
            #    coupling_list.append((_qubit_number(row, col + 1), _qubit_number(row, col)))

    # vertical edges
    for col in range(width):
        for row in range(height - 1):
            if (col+row)%2 == 1:
                patternB.append((_qubit_number(row, col), _qubit_number(row + 1, col)))
            else:
                patternA.append((_qubit_number(row, col), _qubit_number(row + 1, col)))
            #coupling_list.append((_qubit_number(row, col), _qubit_number(row + 1, col)))
            #if directed:
            #    coupling_list.append((_qubit_number(row + 1, col), _qubit_number(row, col)))

    return [patternA,patternB,patternC,patternD,patternC,patternD,patternA,patternB]


