def get_labels(adjpars):

	labels = {
                  'initial_mass': 'INITIAL_MASS',
                  'initial_x':    'INITIAL_X'   ,
                  'initial_y':    'INITIAL_Y'   ,
                  'initial_z':    'INITIAL_Z'   ,
                  'overshoot_f0_above_burn_h_core': 'OVERSHOOT_F0',
                  'overshoot_f_above_burn_h_core': 'EXP_OVERSHOOT_F',
                  'step_overshoot_f_above_burn_h_core': 'STP_OVERSHOOT_F',
                  'min_D_mix':    'MIN_D_MIX',
                  'central_h1':   'CENTRAL_H1',
                  'omega_uni':    'OMEGA_UNI'
                 }


	return_labels = [labels[key] for key in adjpars]

	return return_labels
