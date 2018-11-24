import numpy as np
from scipy.interpolate import PchipInterpolator
import h5py

from entropy.mesa import hdf5_io_support as h5io


def distance_function(xj,wj):
    ## xj --> array of quantities used for calculating distance
    ## wj --> weights per quantity used for calculating distance
    ## npoints --> number of points per track
    ## calculates the formula: D_(i+1) = D_i + sqrt( sum j->N w_j * ( x_(j,i+1) - x_(j,i) )^2 )
    dist = [0.]
    dcurr = 0.
    npoints = len(xj[0])
    for ii in range(1,npoints):
        sqdiff = 0
        for jj in range(len(wj)):
            sqdiff += wj[jj]*(xj[jj][ii]-xj[jj][ii-1])**2
        dcurr = dist[ii-1] + np.sqrt(sqdiff)
        dist.append( dcurr )
    yield dist


def get_xj(primary_eep,xkeys):
    # This function generates the x_j array used in the distance function
    # primary_eep --> a dictionary containing all the information of
    #                 one primary eep of an evolutionary track.
    # xkeys --> the field names for the quantities we want to use in the weight function
    xj = []
    for key in xkeys:
        xj.append( np.array( [ list( get_field_from_primary_eep(primary_eep,key) ) ] ) )
    yield xj


def get_field_from_primary_eep(primary_eep,key):
    # This function generates a single field for the xj array
    # primary_eep --> a dictionary containing all the information of
    #                 one primary eep of an evolutionary track.
    # key --> the field name for the quantity we want to use in the weight function
    # yields data
    if key == 'star_age':
        data = np.log10(primary_eep[key])
    else:
        data = primary_eep[key]
    yield data


def make_interpolation_function(x, y):
    # This function creates an interpolation object to be called later
    i_func = PchipInterpolator(x, y, extrapolate=False)
    yield i_func


def get_tracks_interpolated_on_distance_function(distance,primary_eep):
    # This function makes a dictionary where each field is an interpolation
    # object where the a given field from the original evolutionary track is
    # interpolated on the weighted distance function
    interpolated_track = {}
    for key in primary_eep:
        interpolated_track[key] = make_interpolation_function(distance,primary_eep[key])
    yield interpolated_track


def get_deltas(distance,npoints):
    # This function creates an array of N=npoints equidistantly spaced in distance
    deltas = []
    for ii in range(npoints):
        deltas.append(ii*distance[-1]/float(npoints))
    yield deltas


def calculate_secondary_eeps(interpolated_primary_eep,deltas):
    # This function creates the new equidistantly spaced interpolated
    # evolutionary in all quantities
    secondary_eeps = {}
    for key in interpolated_primary_eep:
        secondary_eeps[key] = interpolated_primary_eep[key](deltas)
    yield secondary_eeps


def generate_secondary_eeps(primary_eep,npoints=150):
    """
    This routine will take a primary eeps per a single track / starself, and
    resample it to generate N=npoints secondary eeps within the given primary eep.
    A given primary eep would be PMS to ZAMS ; ZAMS to MAMS ; MAMS to TAMS ;
    TAMS to RGBTip.


    """

    xkeys=['star_age','log_L','log_Teff','log_center_T','log_center_Rho']
    weights  = [0.05,0.125,2.,1.,1.]

    # Here, we make the x_j array, which consists of all the quantities for which
    # we want to use in our weighted distance function, as listed in the xkeys array.
    xj = list(get_xj(primary_eeps,xkeys))[0]

    # Next, we want to calculate the distance covered by the given primary eep
    distance = list(distance_function(xj,weights))[0]

    # interpolate primary eep in all quantities
    # Only return the interpolated function for each quantity. The new x values
    # still need to be passed to each function for it to run.
    interpolated_primary_eep = list(get_tracks_interpolated_on_distance_function(distance,primary_eep))[0]

    deltas = list(get_deltas(distance,npoints))[0]

    secondary_eeps = list(calculate_secondary_eeps(interpolated_primary_eep,deltas))[0]

    yield secondary_eeps
