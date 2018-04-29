import numpy as np


def parse_header(chain_file):
	header = ''
	f_in = open(chain_file, 'r')
	for l in f_in.readlines():
		if l[0] == '#':
			header += l[2:]
		else:
			break
	f_in.close()
	header = header.strip()

	adjpars, priors = [], []
	for l in header.split('\n'):
		try:
			c = l.split()
			par = c[0].strip()
			llim, ulim = float(c[1]), float(c[2])
			adjpars.append(par)
			priors.append((llim, ulim))
		except:
			continue

	return adjpars, priors

def get_labels(adjpars,passbands = ['0','1','2','3'],perr0=['0','1']):

	labels = {
		  'phoebe_hjd0':'$\\mathrm{HJD_{0}\\,[d]}$',
		  #'phoebe_hjd0':'$\\mathrm{BJD_{0}\\,[d]}$',
		  'phoebe_period':'$\\mathrm{P_{orb}\\,[d]}$',
		  'phoebe_sma':'$\\mathrm{a\\,[R_{\\odot}]}$',
		  'phoebe_rm':'$\\mathrm{q\\,\\frac{M_2}{M_1}}$',
		  'phoebe_vga':'$\\mathrm{\\gamma_0}$',
		  'phoebe_incl':'$i \\,\\mathrm{[\\deg]}$',
		  #'phoebe_perr0':'$\\omega_0-%s \\, \\mathrm{[rad]}$',
		  'phoebe_perr0':'$\\omega_0 \\, \\mathrm{[rad]}$',
		  'phoebe_dperdt':'$\\frac{d\\omega}{dt} \\, \\mathrm{[rad/sec]}$',
		  'phoebe_ecc':'$e$',
		  'phoebe_f1':'$\\frac{P_{rot,1}}{P_{orb}}$',
		  'phoebe_f2':'$\\frac{P_{rot,2}}{P_{orb}}$',
		  'phoebe_teff1':'$\\mathrm{T_{eff,1}\\,[K]}$',
		  'phoebe_teff2':'$\\mathrm{T_{eff,2}\\,[K]}$',
		  'phoebe_pot1':'$\\mathrm{\\Omega_{1}}$',
		  'phoebe_pot2':'$\\mathrm{\\Omega_{2}}$',
		  'phoebe_alb1':'$\\mathrm{A_1}$',
		  'phoebe_alb2':'$\\mathrm{A_2}$',
		  'phoebe_hla':'$\\mathrm{L_{1,%s}}$',
		  'phoebe_cla':'$\\mathrm{L_{2,%s}}$',
		  'phoebe_el3':'$\\mathrm{L_{3,%s}}$',
		  'phoebe_spots_colatitude':'$\\mathrm{Co-lat\\,[rad]}$',
		  'phoebe_spots_longitude':'$\\mathrm{Longitude\\,[rad]}$',
		  'phoebe_spots_radius':'$\\mathrm{R_{spot}\\,[rad]}$',
		  'phoebe_spots_tempfactor':'$\\mathrm{\\frac{T_{spot}}{T_{eff}}}$'
		 }

	hla_counter = 0
	cla_counter = 0
	el3_counter = 0

	return_labels = []
	for par in adjpars:
		print par
		if 'phoebe_hla' in par:
			hla_ind = int(par.split('-')[1])
			return_labels.append(labels['phoebe_hla']%passbands[hla_ind])
		elif 'phoebe_cla' in par:
			cla_ind = int(par.split('-')[1])
			return_labels.append(labels['phoebe_cla']%passbands[cla_ind])
		elif 'phoebe_el3' in par:
			el3_ind = int(par.split('-')[1])
			return_labels.append(labels['phoebe_el3']%passbands[el3_ind])
		#elif 'phoebe_perr0' in par:
		#	perr0_ind = int(par.split('-')[1])
		#	return_labels.append(labels['phoebe_perr0']%perr0[perr0_ind])
		else:
			return_labels.append(labels[par])

	return return_labels
