EMCEE_VERSION = '2.1.0'
import numpy as np
import matplotlib.pyplot as plt
import emcee
from emcee.utils import MPIPool
import os,sys
import subprocess
import tempfile as tmpf
import itertools
import operator
import job_array_commons as jac
import job_array_mesa as jam
from glob import glob
import read_mesa


def write_input_line(switches,arguments):
    
    line = ''
    for jj,switch in enumerate(switches):
        if switches[jj]   == 'logsdir':
            line += '--%s %s '%(switch,arguments[jj])
        elif switches[jj] == 'zamsdir':
            line += '--%s %s '%(switch,arguments[jj])
        elif switches[jj] == 'zams_in_dir':
            line += '--%s %s '%(switch,arguments[jj])
        elif switches[jj] == 'history_file':
            line += '--%s %s '%(switch,arguments[jj])
        elif switches[jj] == 'pulse_file':
            line += '--%s %s '%(switch,arguments[jj])
        else:
            line += '--%s %sd0 '%(switch,arguments[jj])
    return line


def make_arguments(adjpars,adjkeys,switches,tmpf1,tmpf2, logsPath=None,zamsPath=None,zams_inc=None):
    ## this function takes in adjpars, and the sorted adjkeys and switches and creates
    ## the list of arguments to be given, in combination with switches to ./rn
    ## IMPORTANT: This assumes a fixed LOGS dir using temporary history files

    combs = []

    for ii,key in enumerate(adjkeys):

        combs.append(adjpars[key]['value'])
    
    zams_name    = ''
    zams_suffix  = '_ZAMS.mod'

    
    if logsPath is not None:
        logs_name = logsPath  + logs_name
    else:
        logs_name = './LOGS/' + logs_name

    if zamsPath is not None:
        zams_name = zamsPath + zams_name
    else:
        zams_name = './ZAMS_models/'+zams_name

    arguments = np.hstack([logs_name,combs,zams_name+zams_suffix])

    if os.path.exists(logs_name):
        print 'Directory: %s already exists'%(logs_name)
    else:
        os.mkdir(logs_name)

    return arguments
