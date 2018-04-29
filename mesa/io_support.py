EMCEE_VERSION = '2.1.0'
MESA_VERSION  = 'r10108'
import numpy as np
import collections
import mesa.gyre_support as gy
import subprocess


def write_mesa_inlist_from_relax(inlist, adjpars, history_file,
                                 new_Y, new_Z,
                                 pulse_file='', save_pulse_data_logical='.true.',
                                 central_h1=0.0001, x_logic_ctrl_1 = '.true.',
                                 max_age=1e36, nsteps_adjust_before_max_age=0,
                                 base_inlist='/home/cole/python/mesa/inlist_MAMSIE_BASE'):

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

def write_mesa_inlist_from_pms(  inlist, adjpars, history_file,
                                 inital_Y, initial_Z, model_file='',
                                 max_age=1e36, nsteps_adjust_before_max_age=0,
                                 base_inlist='/home/cole/python/mesa/inlist_MAMSIE_BASE_PMS'):

        print 'WRITING INLIST'
        with open(base_inlist, 'r') as f:
                lines = f.readlines()
        replacements = {
                'SAVE_MODEL_FILENAME': '{}'.format("'"+model_file+"'"),
                'STAR_HISTORY_NAME': '{}'.format("'"+history_file+"'"),
                'INITIAL_Y'  : '{:6.5f}'.format(initial_Y),
                'INITIAL_Z'  : '{:6.5f}'.format(initial_Z),
                'INITIAL_MASS'  : '{:8.6f}'.format(adjpars['initial_mass']['value']),
                'MIXING_LENGTH_ALPHA'  : '{:5.4f}'.format(adjpars['mixing_length_alpha']['value']),
                'OVERSHOOT_F_ABOVE_BURN_H_CORE'  : '{:5.4f}'.format(adjpars['overshoot_f_above_burn_h_core']['value']),
                'MIN_D_MIX'  : '{:8.4f}'.format(adjpars['min_D_mix']['value']),
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

def write_mesa_inlist_load_pms(  inlist, adjpars, history_file,
                                 inital_Y, initial_Z, pms_model_file='',
                                 pulse_file='', save_pulse_data_logical='.false.',
                                 profile_file='', save_profile_data_logical='.false.',
                                 central_h1=0.0001, x_logic_ctrl_1 = '.true.',
                                 max_age=1e36, nsteps_adjust_before_max_age=0,
                                 base_inlist='/home/cole/python/mesa/inlist_MAMSIE_BASE_LOAD_PMS'):

        print 'WRITING INLIST'
        with open(base_inlist, 'r') as f:
                lines = f.readlines()
        replacements = {
                'SAVE_PULSE_DATA_LOGICAL': '{}'.format(save_pulse_data_logical),
                'SAVE_PULSE_DATA_FILENAME': '{}'.format("'"+pulse_file+"'"),
                'SAVE_PROFILE_DATA_LOGICAL': '{}'.format(save_pulse_data_logical),
                'SAVE_PROFILE_DATA_FILENAME': '{}'.format("'"+pulse_file+"'"),
                'PMS_MODEL_FILENAME': '{}'.format("'"+pms_model_file+"'"),
                'STAR_HISTORY_NAME': '{}'.format("'"+history_file+"'"),
                'INITIAL_Y'  : '{:6.5f}'.format(initial_Y),
                'INITIAL_Z'  : '{:6.5f}'.format(initial_Z),
                'INITIAL_MASS'  : '{:8.6f}'.format(adjpars['initial_mass']['value']),
                'MIXING_LENGTH_ALPHA'  : '{:5.4f}'.format(adjpars['mixing_length_alpha']['value']),
                'OVERSHOOT_F_ABOVE_BURN_H_CORE'  : '{:5.4f}'.format(adjpars['overshoot_f_above_burn_h_core']['value']),
                'X_LOGICAL_CTRL_1': '{}'.format(x_logic_ctrl_1),
                'CENTRAL_H1' : '{:6.5f}'.format(central_h1),
                'MIN_D_MIX'  : '{:8.4f}'.format(adjpars['min_D_mix']['value']),
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
        




