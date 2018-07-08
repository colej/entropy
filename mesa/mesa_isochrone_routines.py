import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, PchipInterpolator
from glob import glob

from entropy.mesa import hdf5_io_support as h5io
from entropy.general.smoothing import smooth

def distance_function(xj,wj,npoints):
    ## xj --> array of quantities used for calculating distance
    ## wj --> weights per quantity used for calculating distance
    ## npoints --> number of points per track
    ## calculates the formula: D_(i+1) = D_i + sqrt( sum j->N w_j * ( x_(j,i+1) - x_(j,i) )^2 )
    dist = [0.]
    dcurr = 0.
    for ii in range(1,npoints):
        sqdiff = 0
        for jj in range(len(wj)):
            sqdiff += wj[jj]*(xj[jj][ii]-xj[jj][ii-1])**2
        dcurr = dist[ii-1] + np.sqrt(sqdiff)
        dist.append( dcurr )
    return dist


def generate_primary_eeps(track):

    m0 = track["star_mass"][0]
    tags = []
    ## Pre-MS EEP
    pms_     = np.where( track['log_center_T'] > 5.0 )[0][0]
    ## ZAMS EEP
    if m0 <= 1.5:
        zams_     = np.where( ((track['log_Lnuc']/track['log_L'] > 0.999)  &
                               (track['center_h1'][0]-track['center_h1'] < 0.0015)) )[0][0]
    elif ( (m0>1.5) & (m0<=3.0) ):
        zams_     = np.where( ((track['log_Lnuc']/track['log_L'] > 0.9999)  &
                              (track['center_h1'][0]-track['center_h1'] < 0.006) &
                              (np.max(track['log_g'])-track['log_g'] < 0.01 )) )[0][0]
    else:
        zams_     = np.where( ((track['log_Lnuc']/track['log_L'] > 0.9999)  &
                              (track['center_h1'][0]-track['center_h1'] < 0.0015) &
                              (np.max(track['log_g'])-track['log_g'] < 0.01 )) )[0][0]
    # zams_     = np.where( ((track['log_Lnuc']/track['log_L'] > 0.999)  &
    #                       (track['center_h1'][0]-track['center_h1'] < 0.0015)) )[0][0]
    ## IAMS EEP
    iams_     = np.where( track['center_h1'] < 0.3 )[0][0]

    ## TAMS EEP
    tams_      = np.where( track['center_h1'] < 1e-12 )[0][0]
    Yc_at_tams = track['center_he4'][tams_]

    try:
        ## RGB-Tip EEP
        RGBTipL_   = np.where( ((track['log_L']==min(track['log_L'])) &
                               (Yc_at_tams-track['center_he4'] < 0.01)) )[0][0]
        RGBTipT_   = np.where( (track['log_Teff']==min(track['log_Teff'])) and
                               (Yc_at_tams-track['center_he4'] < 0.01) )[0][0]
        rgb_       = min(RGBTipT_,RGBTipL_)


        tags       = [True,True,True,True]

    except:
        rgb_       = -1


        tags       = [True,True,True,False]

    ## make eeps
    ## If genfromtxt
    #eep1 = { key:track[key][pms_ :zams_] for key in track.dtype.names } ## PMS  to ZAMS
    #eep2 = { key:track[key][zams_:iams_] for key in track.dtype.names } ## ZAMS to IAMS
    #eep3 = { key:track[key][iams_:tams_] for key in track.dtype.names } ## IAMS to TAMS
    #eep4 = { key:track[key][tams_:rgb_ ] for key in track.dtype.names } ## TAMS to RGBTip

    eep1 = { key:track[key][pms_ :zams_] for key in track.keys() } ## PMS  to ZAMS
    eep2 = { key:track[key][zams_:iams_] for key in track.keys() } ## ZAMS to IAMS
    eep3 = { key:track[key][iams_:tams_] for key in track.keys() } ## IAMS to TAMS
    eep4 = { key:track[key][tams_:rgb_ ] for key in track.keys() } ## TAMS to RGBTip


    return eep1,eep2,eep3,eep4,tags


def generate_secondary_eeps(primary_eeps,npoints=150):

    xkeys    = ['star_age','log_L','log_Teff','log_center_T','log_center_Rho']
    weights  = [0.05,0.125,2.,1.,1.]

    logAges  = np.array([ np.log10(eep['star_age']) for eep in primary_eeps ])
    logLs    = np.array([ eep['log_L'] for eep in primary_eeps ])
    logTs    = np.array([ eep['log_Teff'] for eep in primary_eeps ])
    logTcs   = np.array([ eep['log_center_T'] for eep in primary_eeps ])
    logRHOcs = np.array([ eep['log_center_Rho'] for eep in primary_eeps ])

    n_eeps = len(primary_eeps)
    xjs  = np.array([ [logAges[ii], logLs[ii], logTs[ii], logTcs[ii], logRHOcs[ii]] for ii in range(n_eeps)])
    lens = np.array( [ len(logAges[ii]) for ii in range( n_eeps ) ] )


    distances_per_eep = []
    for jj,xj in enumerate(xjs):
        d = distance_function(xj,weights,lens[jj])
        distances_per_eep.append(d)

    interpolated_primary_eeps = []

    for ii,eep in enumerate(primary_eeps):

        # interpolated_primary_eeps.append( { key: interp1d( distances_per_eep[ii], eep[key], kind='cubic') for jj,key in enumerate(eep) } )
        interpolated_primary_eeps.append( { key: PchipInterpolator( distances_per_eep[ii], eep[key], extrapolate=False) for jj,key in enumerate(eep) } )

    deltas = []
    for jj,distance in enumerate(distances_per_eep):
        deltas.append([])
        for ii in range(npoints):
            deltas[jj].append( distance[-1]/float(npoints) * ii )

    secondary_eeps = []

    for ii,eep in enumerate(primary_eeps):
        secondary_eeps.append( { key: interpolated_primary_eeps[ii][key]( deltas[ii] ) for key in eep } )

    return secondary_eeps


def construct_eep(track,npoints=100,keys=None):

    ## Find Primary EEPs
    primary_eeps = generate_primary_eeps(track)
    if keys is not None:
        keys = np.array(keys)
    else:
        keys = np.array( [ key for key in primary_eeps[0].keys() ] )

    ## Find Secondary EEPs
    secondary_eeps = generate_secondary_eeps(primary_eeps[:-1], npoints)

    # plt.plot(track['log_Teff'],track['log_g'],'k-')
    # for ii,eep in enumerate(primary_eeps[:-1]):
    #     plt.plot(eep['log_Teff'][0],eep['log_g'][0],'ro',ms=9)
    #     plt.plot(secondary_eeps[ii]['log_Teff'][1:],secondary_eeps[ii]['log_g'][1:],'bo',ms=3)
    #
    # plt.xlim(plt.xlim()[::-1])
    # plt.ylim(plt.ylim()[::-1])
    # plt.show()

    ## Construct single track of eeps
    eep_track    = {}
    for key in keys:
        eep_track[key] = []

    for ii,eep in enumerate(primary_eeps[:-1]):
        for key in keys:
            eep_track[key].extend( np.hstack( [ eep[key][0],secondary_eeps[ii][key][1:],eep[key][-1] ] ) )

    return eep_track


def write_eep(eep,savename):
    ## Function to write eep to file
    with open(savename,'w') as fout:
        keys = eep.keys()
        header = ''
        for key in keys:
            header += '%s\t'%key
        header +='\n'
        fout.write(header)
        npoints = len(eep[keys[0]])
        for ii in range(npoints):
            line = ''
            for key in keys:
                line += '%f\t'%eep[key][ii]
            line += '\n'
            fout.write(line)
    return


def interpolate_mass_tracks(tracks, masses, i_masses):
    keys = tracks[0].keys()
    npoints = len(tracks[0][keys[0]])
    interpolated_tracks = np.zeros( (len( i_masses ),len(keys),npoints) )

    for kk in range(npoints):
        for jj,key in enumerate(keys):
            quantity_at_npoint = [ track[key][kk] for track in tracks ]
            zipped = zip(masses,quantity_at_npoint)
            zipped.sort(key=lambda x:x[0])
            sorted_masses,sorted_quantity = zip(*zipped)
            # i_func = interp1d( masses, quantity_at_npoint, kind='cubic')
            i_func = PchipInterpolator( sorted_masses, sorted_quantity, extrapolate=False)
            i_val  = i_func(i_masses)
            for ii,val in enumerate(i_val):
                interpolated_tracks[ii,jj,kk] = val

    return interpolated_tracks,keys


def construct_isochrones(tracks,i_ages,savename):

    masses  = np.array( [ track['star_mass'][0] for track in tracks] )
    zipp = zip(masses,tracks)
    zipp.sort(key=lambda x:x[0])
    masses,tracks = zip(*zipp)

    npoints = len( tracks[0]['star_age'] )
    # keys    = ['star_age','star_mass','log_L','log_R','log_Teff','log_g',
    #            'log_cntr_T','log_cntr_Rho','log_cntr_P',
    #            'center_h1','center_he4',
    #            'surface_h1','surface_he3','surface_he4','surface_c12','surface_n14',
    #            'Omega_crit','Asymptotic_dP','core_mass_custom','k2_stellar_harmonic'
    #            ]

    #-> Quantities that we will include in our isochrones
    keys    = ['star_age','star_mass','log_L','log_R','log_Teff','log_g',
               'log_cntr_T','log_cntr_Rho','log_cntr_P',
               'center_h1','center_he4','surface_c12','mass_conv_core',
               'Omega_crit','Asymptotic_dP','core_mass_custom','k2_stellar_harmonic'
               ]

    #-> create a dictionary where each age is a dictionary for all quantities listed above
    isochrones = { 'age-%i'%cc: { key: [] for key in keys } for cc,i_age in enumerate(i_ages) }

    #-> Loop through all ages and create the isochrone for that age
    for cc,i_age in enumerate(i_ages):

        print 'ISOCHRONE FOR: ',i_age/1e6,' Myrs'
        mass0s = []
        #-> Loop over each track at the same EEP to generate an isochrone
        for n_eep in range(npoints):

            eep_subsample = { key: [] for key in keys }

            for ii,track in enumerate(tracks):
                #plt.plot(track['log_Teff'],track['log_g'],'k-')
                for key in keys:
                    eep_subsample[key].append(track[key][n_eep])

            eep_masses = np.array( eep_subsample['star_mass'] )
            eep_ages   = np.array( eep_subsample['star_age' ] )


            if ( ( min(eep_ages) < i_age < max(eep_ages) ) & (len(eep_masses) > 3) ) :

                zipped = zip(eep_ages,eep_masses)
                zipped.sort(key=lambda x:x[0])
                sorted_eep_ages,sorted_eep_masses = zip(*zipped)
                # i_age_func     = interp1d( eep_ages, eep_masses, kind='slinear' )
                i_age_func = PchipInterpolator(sorted_eep_ages, sorted_eep_masses, extrapolate=False)
                mass0      = i_age_func(i_age)
                mass0s.append(mass0)

                for key in keys:
                    cval   = np.array(eep_subsample[key])
                    zipped = zip(eep_masses, cval)
                    zipped.sort(key=lambda x:x[0])
                    sorted_eep_masses, sorted_cval = zip(*zipped)
                    i_func = PchipInterpolator(sorted_eep_masses, sorted_cval, extrapolate=False)
                    isochrones['age-%s'%cc][key].append( i_func( mass0 ) )

        #         plt.figure(1)
        #         plt.plot( eep_masses, np.log10(eep_ages), 'ko-')
        #         plt.axvline(i_age_func(i_age),color='red',alpha=0.5)
        #         plt.axhline(np.log10(i_age),color='red',alpha=0.5)
        #
        # plt.show()

                # plt.figure(2)
                # plt.plot(isochrones['age-%i'%cc]['log_Teff'],isochrones['age-%i'%cc]['log_g'],'k-',alpha=0.7)

    for cc,i_age in enumerate(i_ages):
        for key in keys:
            isochrones['age-%i'%cc][key] = np.hstack(isochrones['age-%i'%cc][key])

    # plt.figure(2)
    # plt.xlim(plt.xlim()[::-1])
    # plt.ylim(plt.ylim()[::-1])
    #
    # plt.show()

    with open(savename,'w') as fout:
        header = '# AGE[Myr]  %s'%' '.join( [ '%s'%key for key in keys ] )
        fout.write(header+'\n')
        for cc,i_age in enumerate(i_ages):
            for n in range(len(isochrones['age-%i'%cc]['star_mass'])):
                print n,'/',len(isochrones['age-%i'%cc]['star_mass'])
                fout.write( '%f %s \n'%(i_age,' '.join( [ '%.8f'%isochrones['age-%i'%cc][key][n] for key in keys] )) )




            #valid_eeps.append(np.hstack(valid_masses))

    #for mass in masses:
        #plt.axvline(mass,color='red',alpha=0.5)
    #for i_age in i_ages:
        #plt.axhline(np.log10(i_age),color='blue',alpha=0.5)

    #plt.show()

    return
