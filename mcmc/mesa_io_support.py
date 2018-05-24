EMCEE_VERSION = '2.1.0'
MESA_VERSION  = 'r10108'
import numpy as np
import collections
import subprocess

import entropy.mesa.gyre_support as gy

def get_adjpars_mcmc(file):

        adjpars = collections.OrderedDict()
        for line in open(file,'r'):
                vals = line.split(' ')
                par  = vals[0].strip()
                adjpars[par]={ 'value':0., 'priors':(float(vals[1]),float(vals[2])),
                               'type':vals[3], 'index':vals[4].strip()
                             }
        #adjpars['central_h1_sigma']={'value':0.09,'priors':(0.0005,0.1),
                                     #'type':'flat','index':'primary'}
        return adjpars

def update_adjpars(theta,adjpars):
        ## assumes theta is ordered to agree with adjkeys

        for ii,key in enumerate(adjpars):
                #print key,theta[ii]
                adjpars[key]['value'] = theta[ii]

        return adjpars

def write_mesa_inlist(inlist, adjpars, history_file,
                      new_Y, new_Z,
                      pulse_file='', save_pulse_data_logical='.true.',
                      central_h1=0.0001, x_logic_ctrl_1 = '.true.',
                      max_age=1e36, nsteps_adjust_before_max_age=0,
                      base_inlist='/home/cole/python/entropy/mesa/inlist_MAMSIE_BASE_FE_NET'):

        print 'WRITING INLIST'
        with open(base_inlist, 'r') as f:
                lines = f.readlines()
        replacements = {
                'SAVE_PULSE_DATA_LOGICAL': '{}'.format(save_pulse_data_logical),
                'SAVE_PULSE_DATA_FILENAME': '{}'.format("'"+pulse_file+"'"),
                'STAR_HISTORY_NAME': '{}'.format("'"+history_file+"'"),
                'NEW_Y'  : '{:6.5f}'.format(new_Y),
                'NEW_Z'  : '{:6.5f}'.format(new_Z),
                'INITIAL_MASS'  : '{:8.6f}'.format(adjpars['initial_mass']['value']),
                'MIXING_LENGTH_ALPHA'  : '{:5.4f}'.format(adjpars['mixing_length_alpha']['value']),
                #'STEP_OVERSHOOT_F_ABOVE_BURN_H_CORE'  : '{:5.4f}'.format(adjpars['step_overshoot_f_above_burn_h_core']['value']),
                'OVERSHOOT_F_ABOVE_BURN_H_CORE'  : '{:5.4f}'.format(adjpars['overshoot_f_above_burn_h_core']['value']),
                'MIN_D_MIX'  : '{:8.4f}'.format(adjpars['min_D_mix']['value']),
                'X_LOGICAL_CTRL_1': '{}'.format(x_logic_ctrl_1),
                'CENTRAL_H1' : '{:6.5f}'.format(central_h1),
                'MAX_AGE' : '{:f}'.format(max_age),
                'NSTEPS_ADJUST': '{:1.0f}'.format(nsteps_adjust_before_max_age)
                        }

        new_lines = []
        for line in lines:
                new_line = line
                for key in replacements:
                        if (replacements[key] != ''):
                                new_line = new_line.replace(key, replacements[key])
                new_lines.append(new_line)

        with open(inlist, 'w') as f:
                f.writelines(new_lines)
        return



def run_mesa(run_com,history_file):


        mesa_proc = subprocess.Popen(run_com,bufsize=-1,stdout= subprocess.PIPE, stderr= subprocess.PIPE)

        subout = mesa_proc.stdout.readlines()

        if mesa_proc.returncode is not None:
                print '\t\t\t\tProcess Return code: ',mesa_proc.returncode
                try:
                        mesa_proc.kill()
                        print '\t\tSuccessfully terminated MESA Process'
                except:
                        print 'couldnt kill process'
                return None
        else:
                try:
                        mesa_proc.kill()
                        print '\t\tSuccessfully terminated MESA Process'
                except:
                        print 'Mesa process already terminated'
                print '\t\t\tMESA DONE'

                stop_crits = [ 'termination code: max_age', 'termination code: max_model_number', 'termination code: xa_central_lower_limit']

                if np.any([ crit in out for out in subout for crit in stop_crits ]) :

                        data_star         = np.genfromtxt(history_file,skip_header=5,names=True)

                        teff_star         = 10**data_star['log_Teff'][-1]
                        logg_star         = data_star['log_g'][-1]
                        mass_star         = data_star['star_mass'][-1]
                        logL_star         = data_star['log_L'][-1]
                        radius_star       = 10**data_star['log_R'][-1]
                        xc_star           = data_star['center_h1'][-1]
                        age_star          = data_star['star_age'][-1]
                        omega_crit        = data_star['Omega_crit'][-1]
                        asymptotic_dp     = data_star['Asymptotic_dP'][-1]
                        mass_cc           = data_star['core_mass_custom'][-1]

                        outputs = collections.OrderedDict()
                        outputs['Teff']          = {'value':teff_star}
                        outputs['logg']          = {'value':logg_star}
                        outputs['logL']          = {'value':logL_star}
                        outputs['radius']        = {'value':radius_star}
                        outputs['Xc']            = {'value':xc_star}
                        outputs['age']           = {'value':age_star}
                        outputs['omega_crit']    = {'value':omega_crit}
                        outputs['asymptotic_dp'] = {'value':asymptotic_dp}
                        outputs['mass_cc']       = {'value':mass_cc}

                        return outputs
                else:
                        return None





def init_parameters(chain_file,primary_adjpars,secondary_adjpars):
        """
        Initializes a new chain with passed parameters and priors.
        """
        # Write the header:
        f_out = open(chain_file, 'w')
        f_out.write('# emcee version:   %s\n' % EMCEE_VERSION)
        f_out.write('# \n')
        f_out.write('# Number of parameters being adjusted: %d\n' % (len(primary_adjpars)+len(secondary_adjpars)))
        f_out.write('# \n')
        f_out.write('#     Parameter:   Lower limit:   Upper limit:\n')
        for key in primary_adjpars:
                f_out.write('#  %31s %14.5f %14.5f\n' % (key, primary_adjpars[key]['priors'][0]
                            , primary_adjpars[key]['priors'][1] ) )
        for key in secondary_adjpars:
                f_out.write('#  %31s %14.5f %14.5f\n' % (key, secondary_adjpars[key]['priors'][0]
                            , secondary_adjpars[key]['priors'][1] ) )
        f_out.close()

def init_bookkeeper(filename,header_keys):

        fout = open(filename,'w')
        fout.write('# MCMC Book-keeping file\n')
        fout.write('# emcee version:   %s\n' % EMCEE_VERSION)
        fout.write('# MESA version:   %s\n' % MESA_VERSION)
        header = ''
        for key in header_keys:
                header += key+'\t'
        header += '\n'
        fout.write(header)
        fout.close()

def write_bookkeeper(filename,adjpars,outpars,lnlklhd,frange=10):
        with open(filename,'a') as ftemp:
                bookline = ''
                for key in adjpars:
                        bookline += '%f\t'%adjpars[key]['value']
                        if outpars is not None:
                                for key in outpars:
                                        bookline += '%f\t'%outpars[key]['value']
                                bookline += '%f'%lnlklhd
                        else:
                                for ll in range(frange):
                                        bookline += '-inf\t'
                bookline += '\n'
                ftemp.write(bookline)
        print 'Book-Keeper Updated'
