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
from glob import glob
import read_mesa


def lnlike(in_line):
    ## chi2 minimization
    ## [ FUTURE ] Include noise model

    print 'running chi2'
    try:
        #in_line = mr.write_input_line(switches,arguments)
        run_com = ['/STER/colej/mesa/mesa-r8118/star/kic493a_ms/rn' , in_line ]
    
        output,other = subprocess.Popen(run_com,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
    
        subout = output.split('\n')
        for line in subout:
            print line
    
        if ( ('stop because have dropped below central lower limit for h1') and ('termination code: xa_central_lower_limit') ) in subout:
        
            history_file = in_line.split(' ')[1]
            data = read_mesa.read_mesa_ascii(history_file+'/history.data')
            teff = 10**data[1]['log_Teff'][-1]
            logg = data[1]['log_g'][-1]
            dP   = data[1]['Asymptotic_dP'][-1]*86400.
            
            chisq = (teff-15100.)**2/150.**2 + (logg-3.95)**2/0.1**2 + (dP-6160.)**2/320.**2
            
            print 'Track evaluated! lnL=',-0.5*chisq
            return -0.5*chisq
        
        else:
            print 'Inf returned'
            return -np.inf
    except:
        print 'Except encountered'
        return -np.inf


def get_adjpars_grid(file,delimeter=' '):
    '''
    This function reads in a generic file with a given delimeter of a format Parameter Name - String Format - Format Multiplier -
    Lower Bound - Upper Bound - Step Size - Order Index
    
    '''
    adjpars = {}
    for line in open(file,'r'):
        if line.strip():
            vals = line.split(delimeter)
            par = vals[0]
            print par
            vals = vals[1:]
            adjpars[par.strip()]={'formatter':vals[0],'format_multiplier':float(vals[1]),'limits':(float(vals[2]),float(vals[3])),'step':float(vals[4]),'value': 0. ,'index':int(vals[5])}

    return adjpars

def get_adjpars_sample(file,delimeter=' '):
    from mesa import labels as ml

    '''
    This function reads in a generic file with a given delimeter of a format Parameter Name - prior type -
    Lower Bound - Upper Bound - Order Index

    '''
    adjpars = {}
    for line in open(file,'r'):
        if line.strip():
            vals = line.split(delimeter)
            par = vals[0]
            #print par
            vals = vals[1:]
            adjpars[par.strip()]={'prior_type':vals[0],'limits':(float(vals[1]),float(vals[2])),'value': 0. ,'index':int(vals[3]),'replace_string':'-'}

    strings = ml.get_labels(adjpars)
    for ii,key in enumerate(adjpars):
        adjpars[key]['replace_string']=strings[ii]

    return adjpars

def update_adjpars(theta,adjpars):
    ## assumes theta is ordered to agree with adjpars

    for ii,key in enumerate(adjpars):
        adjpars[key]['value'] = theta[ii]

    return adjpars

def get_sorted_keys_and_switches(adjpars,include_logs=False,include_zams_in=False,include_zams_out=False):
    keys  = adjpars.keys()
    inds = [ adjpars[key]['index'] for key in keys ]
    keys_sort = np.array(sorted(zip(keys,inds),key=lambda x: x[1]))
    
    switches = np.hstack(keys_sort[:,0])
    if include_logs is True:
        switches = np.hstack(['logsdir',switches])
    if include_zams_in is True:
        switches = np.hstack([switches,'zams_in_dir'])
    if include_zams_out is True:
        switches = np.hstack([switches,'zamsdir'])

    return  keys_sort[:,0],switches


def write_input_line(switches,arguments):
    
    line = ''
    for jj,switch in enumerate(switches):
        if switches[jj]   == 'logsdir':
            line += '--%s %s '%(switch,arguments[jj])
        elif switches[jj] == 'zamsdir':
            line += '--%s %s '%(switch,arguments[jj])
        elif switches[jj] == 'zams_in_dir':
            line += '--%s %s '%(switch,arguments[jj])
        elif 'history_file_' in switches[jj]:
            line += '--%s %s '%(switch,arguments[jj])
        else:
            line += '--%s %sd0 '%(switch,arguments[jj])
    return line


#def get_arguments(adjpars,adjkeys,switches):
    

def make_log_dir_single(adjpars,adjkeys,switches,logsPath=None,zamsPath=None,zams_inc=None):
    parts = []
    combs = []
    #keys  = adjpars.keys()
    #print keys
    #inds = [ adjpars[key]['index'] for key in keys ]

    #sortkeys = np.array(sorted(zip(keys,inds),key=lambda x: x[1]))

    for ii,key in enumerate(adjkeys):
        
        #print key
        parts.append(('{}'.format(adjpars[key]['formatter'])).format(adjpars[key]['format_multiplier']*adjpars[key]['value'] ))

        combs.append(adjpars[key]['value'])
        
    #arguments = []
    stack = np.hstack(parts)
    #print stack
    logs_name    = 'LOGS_'
    zams_name    = ''
    #zpath        = '/STER/colej/mesa/mesa-r8118/star/kic493a_pre/ZAMS_models/'
    zams_suffix  = '_ZAMS.mod'
    for jj,dir in enumerate(stack):
        logs_name  += dir
        if ( (zams_inc != None) and (np.any([ inc == sortkeys[jj]  for inc in zams_inc])==True) ):
            zams_name += dir
        elif ( (zams_inc != None) and (np.any([ inc == sortkeys[jj] for inc in zams_inc])==False) ):
            continue
        elif zams_inc == None:
            zams_name += dir
        
    logs_name = logs_name.strip('_')
    zams_name = zams_name.strip('_')
        

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

def lnprior(adjpars):
    # 

    if adjpars['overshoot_f_above_burn_h_core']['value'] < adjpars['overshoot_f0_above_burn_h_core']['value']:
	print '\t fov < f0,ov'
	return -np.inf
    else:
        if not np.all( [ adjpars[key]['limits'][0] < theta[ii] < adjpars[key]['limits'][1] for ii,key in enumerate(adjkeys) ] ):
    	    print '\tPriors Out Of Bounds!!!'
	    return -np.inf
	else:
            print '\tPriors Satisfied'
            #lnp = 0.
            return 0.


def init_parameters(chain_file,adjpars,adjkeys=None):
    """
    Initializes a new chain with passed parameters and priors.
    """
    #p0s = np.array([[p[0] + (p[1]-p[0])*np.random.rand() for p in priors] for i in xrange(nwalkers)])
    # Write the header:
    f_out = open(chain_file, 'w')
    f_out.write('# emcee version:   %s\n' % EMCEE_VERSION)
    #f_out.write('# PHOEBE version: %s\n' % PHOEBE_VERSION)
    f_out.write('# \n')
    f_out.write('# Number of parameters being adjusted: %d\n' % len(adjpars))
    f_out.write('# \n')
    f_out.write('#     Parameter:   Lower limit:   Upper limit:\n')
    if adjkeys is not None:
        for key in adjkeys:
            #print key
            #print adjpars[key]['limits'][0],adjpars[key]['limits'][1]
            f_out.write('#%15s %14.5f %14.5f\n' % (key, adjpars[key]['limits'][0], adjpars[key]['limits'][1] ) )
        f_out.close()
    else:
        for key in adjpars:
            f_out.write('#%15s %14.5f %14.5f\n' % (key, adjpars[key]['limits'][0], adjpars[key]['limits'][1] ) )
        f_out.close()
