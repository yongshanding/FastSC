"""
device.py - A module for defining attributes in the device..
"""

#from fqc.util import get_unitary
import numpy as np
import scipy.linalg as la
from scipy.special import factorial
import os,sys,inspect
import h5py
import matplotlib.pyplot as plt
import random as rd
import time
import re, math
from datetime import datetime


class Device(object):
    """
    Fields:
    """
    def __init__(self, side_length, omega_max, delta_int, delta_ext, delta_park):
        self.side_length = side_length
        self.omega_max = omega_max
        self.delta_int = delta_int
        self.delta_ext = delta_ext
        self.delta_park = delta_park
