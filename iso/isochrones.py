import numpy as np
from scipy.interpolate import PchipInterpolator
import h5py

from entropy.mesa import hdf5_io_support as h5io


def construct_isochrones(eep_tracks,i_ages,keys,savename):


    masses  = np.array( [ eep_track['star_mass'][0] for eep_track in eep_tracks] )
    zipp = zip(masses,tracks)
    zipp.sort(key=lambda x:x[0])
    masses,tracks = zip(*zipp)

    npoints = len( tracks[0]['star_age'] )

    # #-> Quantities that we will include in our isochrones
    # keys    = ['star_age','star_mass','log_L','log_R','log_Teff','log_g',
    #            'log_cntr_T','log_cntr_Rho','log_cntr_P',
    #            'center_h1','center_he4','surface_c12','mass_conv_core',
    #            'Omega_crit','Asymptotic_dP','core_mass_custom','k2_stellar_harmonic'
    #            ]

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
