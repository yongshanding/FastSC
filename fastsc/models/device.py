"""
device.py - A module for defining attributes in the device..
"""

#from fqc.util import get_unitary
import numpy as np
import scipy.linalg as la
from scipy.special import factorial
import os,sys,inspect
# import h5py
import matplotlib.pyplot as plt
import random as rd
import time
import re, math
from datetime import datetime
from ..util import get_grid_coupling_list, get_connectivity_graph, get_crosstalk_graph


GATETIMES = {'unitary': 55,'rz': 30, 'z':30, 'u1': 30, 's': 30, 't': 30, 'rx': 30, 'x': 30, 'u2': 30, 'ry': 30, 'y': 30, 'u3': 30, 'h': 30, 'measure': 0.0, 'barrier': 0.0} # all in ns

class Device(object):
    """
    Fields:
    """
    def __init__(self, topology, qubits, omega_max, delta_int, delta_ext, delta_park, cqq=0.019, alpha=-0.2, EJS=8, EJL=20, EC=0.5, flux_sigma=0.0, error_1q_gate=0.001, d=1):
        if (topology == None):
            self.topology = 'grid'
        else:
            self.topology = topology
        self.qubits = qubits
        self.side_length = int(np.sqrt(qubits)) # only for square grid
        self.omega_max = omega_max
        self.delta_int = delta_int
        self.delta_ext = delta_ext
        self.delta_park = delta_park
        self.omega_min = omega_max - delta_int - delta_ext - delta_park
        self.cqq = cqq
        self.alpha = alpha
        self.ejs = EJS
        self.ejl = EJL
        self.ec = EC
        self.flux_sigma = flux_sigma
        self.g_connect = None
        self.coupling = None
        self.g_xtalk = None
        self.d_xtalk = d
        self.gate_times = GATETIMES
        self.error_1q_gate = error_1q_gate

    def build_graph(self, param=None):
        if self.topology=='grid':
            self.g_connect = get_connectivity_graph(self.qubits, self.topology)
            self.coupling = get_grid_coupling_list(self.side_length, self.side_length)
            self.g_xtalk = get_crosstalk_graph(self.g_connect, self.topology, self.d_xtalk)
        elif self.topology == 'erdosrenyi':
            self.g_connect = get_connectivity_graph(self.qubits, self.topology, param)
            self.coupling = self.g_connect.edges()
            self.g_xtalk = get_crosstalk_graph(self.g_connect, self.topology, self.d_xtalk)
        else:
            print("Topology not recognized.")
            self.coupling = None

class Sycamore_device(object):
    def __init__(self, device, size, res_coupling=0.0):
        if size != device.qubits:
            print("Warning: device size inconsistent. Device: " + str(device.qubits) + ", Sycamore_device: " + str(size))
        if size not in [4,9,16,25]:
            raise Exception("Wrong device size")
        if size == 4:
            self.int_freqs = {(0,1):6.619, (0,2):6.65, (1,3):6.657, (2,3):6.677}
            self.park_freqs = {0:6.605,1:6.638,2:6.694,3:6.681}
        elif size == 9:
            self.park_freqs = {0:6.605,1:6.638,2:6.565,3:6.694,4:6.681,5:6.601,6:6.643,7:6.621,8:6.646}
            self.int_freqs = {(0,1):6.619, (1,2):6.601, (3,4):6.677, (4,5):6.635, (6,7):6.631, (7,8):6.646, (0,3):6.65, (1,4):6.657, (2,5):6.585, (3,6):6.667, (4,7):6.645,(5,8):6.646}
        elif size == 16:
            self.park_freqs = {0:6.605,1:6.638,2:6.565,3:6.555,4:6.694,5:6.681,6:6.601,7:6.626,8:6.643,9:6.621,10:6.646,11:6.657,12:6.712,13:6.671,14:6.586,15:6.623}
            self.int_freqs = {(0,1):6.619,(1,2):6.601,(2,3):6.565,(4,5):6.677,(5,6):6.635,(6,7):6.595}
            self.int_freqs[(8,9)] = 6.631
            self.int_freqs[(9,10)] = 6.646
            self.int_freqs[(10,11)] = 6.646
            self.int_freqs[(12,13)] = 6.69
            self.int_freqs[(13,14)] = 6.631
            self.int_freqs[(14,15)] = 6.623
            self.int_freqs[(0,4)] = 6.65
            self.int_freqs[(1,5)] = 6.657
            self.int_freqs[(2,6)] = 6.585
            self.int_freqs[(3,7)] = 6.592
            self.int_freqs[(4,8)] = 6.667
            self.int_freqs[(5,9)] = 6.645
            self.int_freqs[(6,10)] = 6.646
            self.int_freqs[(7,11)] = 6.642
            self.int_freqs[(8,12)] = 6.68
            self.int_freqs[(9,13)] = 6.645
            self.int_freqs[(10,14)] = 6.646
            self.int_freqs[(11,15)] = 6.633
        else: # size == 25
            self.park_freqs = {0:6.612,1:6.571,2:6.605,3:6.638,4:6.565,5:6.687,6:6.661,7:6.694,8:6.681,9:6.601,10:6.634,11:6.628,12:6.643,13:6.621,14:6.646,15:6.707,16:6.665,17:6.712,18:6.671,19:6.586,20:6.775,21:6.734,22:6.766,23:6.729,24:6.594}
            self.int_freqs = {(0,1):6.592,(1,2):6.589,(2,3):6.619,(3,4):6.601}
            self.int_freqs[(5,6)] = 6.675
            self.int_freqs[(6,7)] = 6.678
            self.int_freqs[(7,8)] = 6.677
            self.int_freqs[(8,9)] = 6.635
            self.int_freqs[(10,11)] = 6.628
            self.int_freqs[(11,12)] = 6.632
            self.int_freqs[(12,13)] = 6.631
            self.int_freqs[(13,14)] = 6.646
            self.int_freqs[(15,16)] = 6.693
            self.int_freqs[(16,17)] = 6.690
            self.int_freqs[(17,18)] = 6.690
            self.int_freqs[(18,19)] = 6.631
            self.int_freqs[(20,21)] = 6.753
            self.int_freqs[(21,22)] = 6.766
            self.int_freqs[(22,23)] = 6.743
            self.int_freqs[(23,24)] = 6.594

            self.int_freqs[(0,5)] = 6.648
            self.int_freqs[(1,6)] = 6.617
            self.int_freqs[(2,7)] = 6.650
            self.int_freqs[(3,8)] = 6.657
            self.int_freqs[(4,9)] = 6.585
            self.int_freqs[(5,10)] = 6.660
            self.int_freqs[(6,11)] = 6.640
            self.int_freqs[(7,12)] = 6.667
            self.int_freqs[(8,13)] = 6.645
            self.int_freqs[(9,14)] = 6.646
            self.int_freqs[(10,15)] = 6.670
            self.int_freqs[(11,16)] = 6.646
            self.int_freqs[(12,17)] = 6.680
            self.int_freqs[(13,18)] = 6.645
            self.int_freqs[(14,19)] = 6.646
            self.int_freqs[(15,20)] = 6.741
            self.int_freqs[(16,21)] = 6.704
            self.int_freqs[(17,22)] = 6.712
            self.int_freqs[(18,23)] = 6.703
            self.int_freqs[(19,24)] = 6.594

        self.res_coupling = res_coupling
        omega_max = max(self.int_freqs.values())
        omega_min = min(self.int_freqs.values())
        if omega_max > device.omega_max:
            print("Warning: max freq inconsistent. Device: " + str(device.omega_max) + ", Sycamore_device: " + str(omega_max))
            self.omega_max = omega_max
        if omega_min < device.omega_min:
            print("Warning: min freq inconsistent. Device: " + str(device.omega_min) + ", Sycamore_device: " + str(omega_min))
            self.omega_min = omega_min

    def get_park_freq(self,q):
        if q not in self.park_freqs:
            raise Exception("Wrong qubit number")
        return self.park_freqs[q]

    def get_itr_freq(self,q1,q2):
        if (q1,q2) in self.int_freqs:
            return self.int_freqs[(q1,q2)]
        elif (q2,q1) in self.int_freqs:
            return self.int_freqs[(q1,q2)]
        else:
            raise Exception("Coupling not found")
