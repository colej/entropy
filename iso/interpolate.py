import numpy as np
from scipy.interpolate import interp1d

from entropy.mesa.mesa_isochrone_routines import hdf5_io_support as h5io

def load_track(file,skip_header=0,type='ascii'):

    if type=='ascii':
        track = np.genfromtxt(file,skip_header=skip_header,names=True)
    elif type=='hdf5':
        track = h5io.read_hdf5_history(file)

    return track

def interpolate_on_age(track,delta_t,keys=["M","logL","logTeff","logR","logroeff","logM_loss"],kind='slinear'):
    ## @input: track   --> evolutionary track with given
    ## @input: delta_t --> time step to interpolate evolutionary track to
    ## @input: keys    --> fields of the track that we want to interpolate along
    ## @input: kind    --> type of interoplation to use. Options are: nearest; linear; zero;
    ##                      slinear; quadratic; cubic

    ## Sanity Check
    if not kind in ['nearest', 'linear', 'zero', 'slinear', 'quadratic','cubic']:
        raise ValueError, "kind must be either 'nearest', 'linear', 'zero', 'slinear', 'quadratic','cubic'"

    ## Determine initial mass of track
    Mini = track['M'][0]

    ## Build array of ages sampled equidistantly from the beginning to end of
    ## the track with time step delta_t
    age    = track['age']
    i_ages = np.arange(min(age),max(age),delta_t)

    ## Build the interpolation function for each field of interest
    ## and use to populate new interpolated track
    i_track = { key: [] for key in keys }
    for key in keys:
        ifunc = interp1d(age, track[key], kind=kind)(i_ages)
        i_track[key] = ifunc

    return i_track
