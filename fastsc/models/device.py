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
from ..util import get_nearest_neighbor_coupling_list

class Device(object):
    """
    Fields:
    """
    def __init__(self, side_length, omega_max, delta_int, delta_ext, delta_park, cqq=0.012, alpha=-0.2, EJS=8, EJL=20, EC=0.3, flux_sigma=0.0, coupling=None):
        self.qubits = side_length*side_length
        self.side_length = side_length
        self.omega_max = omega_max
        self.delta_int = delta_int
        self.delta_ext = delta_ext
        self.delta_park = delta_park
        self.cqq = cqq
        self.alpha = alpha
        self.ejs = EJS
        self.ejl = EJL
        self.ec = 0.3
        self.flux_sigma = flux_sigma
        if coupling==None:
            self.coupling = get_nearest_neighbor_coupling_list(side_length, side_length)
        else:
            self.coupling = coupling

class Sycamore_device(object):
    def __init__(self, size):
        if size not in [4,9,16]:
            raise Exception("Wrong device size")
        if size == 4:
            self.int_freqs = {(0,1):6.619, (0,2):6.65, (1,3):6.657, (2,3):6.677}
            self.park_freqs = {0:6.605,1:6.638,2:6.694,3:6.681}
        elif size == 9:
            self.park_freqs = {0:6.605,1:6.638,2:6.565,3:6.694,4:6.681,5:6.601,6:6.643,7:6.621,8:6.646}
            self.int_freqs = {(0,1):6.619, (1,2):6.601, (3,4):6.677, (4,5):6.635, (6,7):6.631, (7,8):6.646, (0,3):6.65, (1,4):6.657, (2,5):6.585, (3,6):6.667, (4,7):6.645,(5,8):6.646}
        else:
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
