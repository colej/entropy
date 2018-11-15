import numpy as np


def time_to_ph(time, period=None, t0=None):
    '''
    converts time to phase from input ephemeris
    DOES NOT ACCOUNT FOR BARYCENTRIC OR HELIOCENTRIC CORRECTION

    hjd_to_ph(time, period, t0)

    input: time (float or array)
    input: period (float)
    input: t0 (float)
    output: phase (float or array)
    '''

    if t0 > 2400000:
        t0 -= 2400000.
    if time[0] > 2400000:
        time = np.array( [t-2400000. for t in time] )
    ph = np.array( [ -0.5+( ( t-t0-0.5*period ) % period ) / period for t in time ] )
    return ph


def ph_to_time (phase, period, t0, norb):
    '''
    converts phase to time from an input ephemeris and number of orbit norb
    time = ph_to_time(phase, period, t0, norb)
    input: phase (float or array)
    input: period (float)
    input: t0 (float)
    input: n (int)
    output: time (float or array)
    '''
    n = int(n)

    time = np.array( [ ph*period + t0 + norb*period for ph in phase ] )
    #~ (phase+0.5)*period+period/2.0+hjd0+n*period
    return time


def calc_sup_conj_phase(omega,eps):
    '''
    routine to calculate the phase of superior conjuction
    and phase of periastron

    '''

    ecc_anomoly    = 2.*np.arctan( np.tan( 0.25*np.pi - 0.5*omega )*((1.+eps)/(1.-eps))**-0.5 )
    mean_anomoly   = ecc_anomoly-eps*np.sin(ecc_anomoly)

    phase_sup_conj = 0.5*(mean_anomoly+omega)/np.pi - 0.25
    phase_periast  = phase_sup_conj - 0.5*mean_anomoly/np.pi

    return phase_sup_conj,phase_periast
