import numpy as np
from scipy.interpolate import interp1d

def load_track(file,skip_header=0,type='ascii'):
    ## load starevol tracks

    if type=='ascii':
        track = np.genfromtxt(file,skip_header=skip_header,names=True)
        npms  = np.where( track["Ph"] != 1 ) ## --> skips the PMS phase
        track = track[npms]

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


def weighted_mass(obs,tracks,eval_keys=["logL","logTeff"]):
    ## @input: obs --> dictionary containing 'value', 'upper', and 'lower' per key
    ## @input: tracks --> tracks that will be considered in the weighted averaging
    ## @input: eval_keys --> keys for the values that will be evaluated
    ##
    ## @output: mass --> weighted gaussian average of the mass for all points considered
    ## @output: imass --> weighted gaussian average of the initial mass for all points considered
    ##
    ## It is important to note that the normalization factors 1/(2*sigma_L * sigma_Teff) cancel
    ## out. As such, we do not use them.


    ## Calculate normalization factor rho_(T,L) --> Eq (3) from Escorza et al. 2017
    cost  = 0. # --> cost for tracks considering masses per point
    icost = 0. # --> cost for tracks considering the initial mass of the tracks
    for track in tracks:

        for ii in range(len(track['M'])):
            track_cost = 1.
            for key in eval_keys:

                sigma = ( obs[key]['lower'] * ( track[key][ii] <= obs[key]['value'] )
                        + obs[key]['upper'] * ( track[key][ii] >= obs[key]['value'] ) )

                track_cost *= np.exp( -0.5* (obs[key]['value'] - track[key][ii])**2
                                           / sigma**2 )

            cost  += track_cost * track['M'][ii]
            icost += track_cost * track['M'][0]

    inverse_rho  = 1. / cost
    inverse_irho = 1. / icost

    ## Calculate weighted mass --> Eq (2) from Escorza et al. 2017
    wMass  = 0.
    wiMass = 0.
    for track in tracks:
        for ii in range(len(track['M'])):
            track_cost = 1.
            for key in eval_keys:
                sigma = ( obs[key]['lower'] * ( track[key][ii] <= obs[key]['value'] )
                        + obs[key]['upper'] * ( track[key][ii] >= obs[key]['value'] ) )

                track_cost *= np.exp( -0.5* (obs[key]['value'] - track[key][ii])**2
                                           / sigma**2 )

            track_cost *= (track['M'][ii]**2)*inverse_rho
            wMass  += track_cost * (track['M'][ii]**2)*inverse_rho
            wiMass += track_cost * (track['M'][0]**2)*inverse_irho

    return wMass,wiMass
