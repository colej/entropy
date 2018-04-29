import numpy as np
import keppy as kp


def run_binning(x,y,yerr=None,nbins=100,alias=False):
	# Binning function -- takes into account aliasing and error
	# propogation on errorbins
	bwidth = 1./nbins
	if alias==True:
		phStart,phStop = -0.6,0.6
	else:
		phStart,phStop = -0.5,0.5

	bins      = np.arange(phStart,phStop+bwidth,bwidth)
	bin_means = ( np.histogram(x,bins,weights=y)[0] / np.histogram(x,bins)[0] )
	if yerr is not None:
		bin_errs = ( np.histogram(x,bins,weights=yerr)[0] / np.histogram(x,bins)[0] )
	else:
		bin_errs = None

	return bwidth,bins,bin_means,bin_errs

def oversample_binning(x,y,yerr=None,phStart=-0.5,phStop=0.5,nbins=100):
	# Binning function -- takes into account aliasing and error
	# propogation on errorbins
	bwidth = (abs(phStart - phStop))/nbins

	print 'BWIDTH: ',bwidth
	bins      = np.arange(phStart,phStop+bwidth,bwidth)
	bin_means = ( np.histogram(x,bins,weights=y)[0] / np.histogram(x,bins)[0] )
	if yerr is not None:
		bin_errs = ( np.histogram(x,bins,weights=yerr)[0] / np.histogram(x,bins)[0] )
	else:
		bin_errs = None

	return bwidth,bins,bin_means,bin_errs


