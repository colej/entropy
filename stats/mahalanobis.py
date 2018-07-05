import numpy as np
import itertools


def MahalanobisDistance(YObs,EpsObs,Thetas,YTheo,YOther):
    #-> Original Author: May Gade Pedersen
    #-> Updated: Cole Johnston 2018-06-27

    # N: number of observed values
    # p: number of varied parameters in the grid (e.g. Mini, Xini, Xc etc.)
    # q: number of grid points
    # YObs: observed values, e.g. frequencies (a vector of length N)
    # EpsObs: errors on the observed values (a vector of length N)
    # Thetas: Matrix of parameters in the grid (dimensions q x p)
    # YTheo: corresponding theoretical values of the observed ones
    # (a vector of length N)

    # Determine number of grid points
    q = np.shape(Thetas)[0]

    # Convert to matrix format
    YObsMat = np.matrix(YObs).T
    YTheoMat = np.matrix(YTheo).T
    YOtherMat = np.matrix(YOther).T

    # Calculate the average on the theoretical values (e.g. frequencies)
    # over the entire grid. Returns a vector of dimension N x 1
    Yav = YTheoMat.mean(1)

    # Calculate the variance-covriance matrix
    N = len(YObs)
    V = np.zeros((N,N))

    for i in range(q):
        difference = np.subtract(YTheoMat[:,i],Yav)
        V += np.matmul(difference,difference.T)
    V = V/(q-1.)

    # Include observational errors in the variance-covariance matrix
    V = V + np.diag(EpsObs**2.)

    # Calculate Mahalanobis distances
    MD = np.zeros(q)
    Vinv = np.linalg.pinv(V,rcond=1e-12)
    for i in range(q):
        diff = (YTheoMat[:,i]-YObsMat)
        MD[i] = np.matmul(np.matmul(diff.T,Vinv),diff)[0][0]

    # return the results ordered from lowest to highest Mahalanobis
    # distance
    idx2 = np.argsort(MD)
    return np.concatenate((np.matrix(MD[idx2]).T,Thetas[idx2,:],YTheoMat[:,idx2].T,YOtherMat[:,idx2].T),axis=1)



def MahalanobisDistanceEB(YObsArr,EpsObsArr,Thetas,YTheo,YOther):
    #-> Original Author: May Gade Pedersen
    #-> Updated: Cole Johnston 2018-06-27

    # N: number of observed values
    # p: number of varied parameters in the grid (e.g. Mini, Xini, Xc etc.)
    # q: number of grid points
    # YObsArr: array of observed values per star, e.g. frequencies (a vector of length N)
    # EpsObsArr: array of errors on the observed values per star (a vector of length N)
    # Thetas: Matrix of parameters in the grid (dimensions q x p)
    # YTheo: corresponding theoretical values of the observed ones
    # (a vector of length N)

    # Determine number of grid points
    q = np.shape(Thetas)[0]

    YObsMatArr = [ np.matrix(YObs).T for YObs in YObsArr ]
    # Convert to matrix format
    YObsMat = np.matrix(YObs).T
    YTheoMat = np.matrix(YTheo).T
    YOtherMat = np.matrix(YOther).T

    # Calculate the average on the theoretical values (e.g. frequencies)
    # over the entire grid. Returns a vector of dimension N x 1
    Yav = YTheoMat.mean(1)

    # Calculate the variance-covriance matrix
    N = len(YObs)
    V = np.zeros((N,N))

    for i in range(q):
        difference = np.subtract(YTheoMat[:,i],Yav)
        V += np.matmul(difference,difference.T)
    V = V/(q-1.)


    # Include observational errors in the variance-covariance matrix
    VInvArr = [ np.linalg.pinv(V + np.diag(EpsObs**2.)) for EpsObs in EpsObsArr ]

    # Calculate Mahalanobis distances
    MDArr = [ np.zeros(q) for xx in YObsArr ]
    MDTot = np.zeros(q)
    # Vinv = np.linalg.pinv(V,rcond=1e-12)
    for i in range(q):
        diffArr = [ (YTheoMat[:,i]-YObsMat) for YObsMat in YObsMatArr ]
        for jj,MD in enumerate(MDArr):
            MD[i] = np.matmul(np.matmul(diffArr[jj].T,VInvArr[jj]),diffArr[jj])[0][0]



    # return the results ordered from lowest to highest Mahalanobis
    # distance
    idx2 = np.argsort(MD)

    return np.concatenate((np.matrix(MD[idx2]).T,Thetas[idx2,:],YTheoMat[:,idx2].T,YOtherMat[:,idx2].T),axis=1)
