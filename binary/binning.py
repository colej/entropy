import numpy as np


def run_aliasing(x,y,z):
    ax,ay,az = [],[],[]
    for ii in range(len(x)):
        if x[ii] < -0.4:
                ax.append(x[ii]+1)
                ay.append(y[ii])
                az.append(z[ii])
        if x[ii] > 0.4:
                ax.append(x[ii]-1)
                ay.append(y[ii])
                az.append(z[ii])
    return np.append(x,ax),np.append(y,ay),np.append(z,az)

def run_binning(x,y,yerr=None,nbins=100,alias=False):
    # Binning function -- takes into account aliasing and error
    # propogation on errorbins
    bwidth = 1./nbins
    phStart,phStop = -0.5,0.5

    bins = np.arange(phStart,phStop+bwidth,bwidth)

    bin_means = ( np.histogram(x,bins,weights=y)[0] / np.histogram(x,bins)[0] )

    if yerr is not None:
        bin_errs = ( np.histogram(x,bins,weights=yerr)[0] / np.histogram(x,bins)[0] )
        if alias==True:
            bins,bin_means,bin_errs = run_aliasing(bins,bin_means,bin_errs)
    else:
        bin_errs = None
        if alias==True:
            bins,bin_means,_ = run_aliasing(bins,bin_means,np.ones_like(bin_means))

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
