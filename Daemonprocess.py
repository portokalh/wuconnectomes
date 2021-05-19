#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  18 10:38:04 2020

@author: Jacques Stout
Multiprocessing class
Used to have multiprocessing within a multiprocessing function
"""
import multiprocessing

from types import ModuleType, FunctionType



class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.

# Custom objects know their class.
# Function objects seem to know way too much, including modules.
# Exclude modules as well.
BLACKLIST = type, ModuleType, FunctionType

class MyPool(multiprocessing.pool.Pool):
    #This whole class allows us to run sub processes within subprocesses (subprocesses cannot exceed max count though)
    Process = NoDaemonProcess
