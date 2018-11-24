import numpy as np
from scipy.interpolate import PchipInterpolator
import h5py

from entropy.mesa import hdf5_io_support as h5io
from entropy.evol import primary_eeps as eepe
from entropy.evol import secondary_eeps as eese


def construct_eep(track,npoints=100,keys=None):

    if keys is not None:
        keys = np.array(keys)
    else:
        keys = np.array( [ key for key in track.keys() ] )

    ## Find Primary EEPs
    primary_eeps = list(eepe.generate_primary_eeps(track))
    pms_to_zams_peep, zams_to_mams_peep, mams_to_tams_peep, tams_to_rgbtip_peep = primary_eeps

    ## Find Secondary EEPs
    pms_to_zams_seep = eese.generate_secondary_eeps(pms_to_zams_peep,npoints=npoints)
    zams_to_mams_seep = eese.generate_secondary_eeps(zams_to_mams_peep,npoints=npoints)
    mams_to_tams_seep = eese.generate_secondary_eeps(mams_to_tams_peep,npoints=npoints)
    tams_to_rgbtip_seep = eese.generate_secondary_eeps(tams_to_rgbtip_peep,npoints=npoints)


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
        eep_track[key] = np.hstack([ pms_to_zams_seep[key], zams_to_mams_seep[key],
                                     mams_to_tams_seep[key], tams_to_rgbtip_seep[key] ])

    yield eep_track


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
