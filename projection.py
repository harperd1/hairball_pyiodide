# from scipy import optimize
import numpy
np = numpy
import pandas as pd

def distance(r1, r2):
    """Euclidean distance between points.
        
        Parameters
        ----------
            r1 : list
                Coordinates of point 1.
            r2 : list
                Coordinates of point 2.
            
        Returns
        -------
            dist : float
                Euclidean distance between points 1 and 2.

    """
    dx = r1[0] - r2[0]
    dy = r1[1] - r2[1]
    d = (dx ** 2 + dy ** 2)**.5
    return d

# def DistErr(x, *args):
#     """Computes distance error function for scipy optimization. 
#     Copied from E3 in pp. 311 of ref. [1]
        
#     Parameters
#     ----------
#         x : np.array
#             1D array of coordinates to be optimized.
#         *args : dict
#             Other parameters (refer to scipy.optimize docs)
        
#     Returns
#     -------
#         E : np.array
#             Objective function

#     """
#     E = 0
#     distance_matrix,natoms = args
#     for i in range(natoms-1):
#         for j in range(i+1, natoms):
#             ri = [x[2*i], x[2*i+1]] #coordinates of atom1
#             rj = [x[2*j], x[2*j+1]] #coordinates of atom2
#             dij = distance(ri, rj)
#             intended_distance = distance_matrix[i][j]
#             E += np.abs(intended_distance-dij)**2
#     return np.asarray(E)

# def optimize():
#     x = np.reshape(initial_guess, 2*natoms)
#     if maxiter:
#         res1 = optimize.fmin_cg(DistErr, x, args=(distance_matrix, natoms), maxiter = maxiter)
#     else:
#         res1 = optimize.fmin_cg(DistErr, x, args=(distance_matrix, natoms))
#     X = np.reshape(res1, (natoms, 2))


def buildDistanceMatrix(connectivityMatrix):
    #Build a distance matrix such that the distance between each node is equal
    ### To the length of the shortest path between them
    ### Not based on any particular research - just seems like a reasonable approach

    def findPathLenth(idx1,idx2,pathLength=1,connectionsChecked=[],connectivityMatrix=connectivityMatrix,):

        directConnection = connectivityMatrix[idx1,idx2]
        if directConnection > 0.1:
            return pathLength

        else:
            intermediateConnections = [i for i,ii in enumerate(connectivityMatrix[idx1,:]) if ii>.1]
            intermediateConnections = [i for i in intermediateConnections if i not in connectionsChecked]
            if len(intermediateConnections) == 0:
                return 1000
            pathLength = min([findPathLenth(i,idx2,pathLength=pathLength+1,connectionsChecked=connectionsChecked+[idx1]) for i in intermediateConnections])
            return pathLength

 
    nNodes = connectivityMatrix.shape[0]
    distanceMatrix = np.zeros((nNodes,nNodes))

    for i in range(nNodes):
        for ii in range(nNodes):
            if i==ii:
                continue
            distanceMatrix[i,ii] = findPathLenth(i,ii)

    return distanceMatrix

def __GetCMDists(distanceMatrix):
    # Find the distance of each node to the "center of mass" of the matrix
    ### The goal is to get a represenation of the points that doesn't directly reference their positions
    ### This enables some clever math tricks and we'll use it later

    ### This is an implementation of equation 48 from the referenced paper, but I don't really get how it works
    ### I'm taking it on faith that the equation works.
    ### I have added a step so it returns the center mass distance at the end instead of the squared center mass distance
    ###### I only made that change because it's more intuitive to me

    natoms = distanceMatrix.shape[0]
    squaredCenterMassDistance = np.zeros(natoms)
    centerMassDistance = np.zeros(natoms)
    for i in range(natoms):
        for j in range(natoms):
            squaredCenterMassDistance[i] += distanceMatrix[i][j]**2/natoms
        for j in range(natoms):
            for k in range(j, natoms):
                squaredCenterMassDistance[i] -= (distanceMatrix[j][k])**2/natoms**2
        centerMassDistance[i] = squaredCenterMassDistance[i]**.5
    return centerMassDistance

def __buildMatrixA(distanceMatrix, centerMassDistance):
    # Build a matrix "A" as it appears in equation 51 from the referenced paper
    ### This euqation is explained in the text between equations 50 and 51
    ### Since we are not using any weights in our equation, this is also the same as matrix B in equation 54
    ### The paper proves that the best-fit for minimizing the "Strain" of fitting a distance matrix in
    ###### low-dimensions can be found based on the eigenvectors of this matrix
    ############### "Therefore, the global minimum Y subject to this condition 
    ############### is obtained by taking the three largest nonnegative eigenvalues of B, 
    ############### and scaling the corresponding eigenvectora by their square roots"

    natoms = distanceMatrix.shape[0]
    A = np.zeros((natoms, natoms))
    for i in range(natoms):
        for j in range(natoms):
            A[i][j] = (centerMassDistance[i]**2 + 
                       centerMassDistance[j]**2 - 
                       distanceMatrix[i][j]**2)/2
    return A

def findFitIn2Dimensions(distanceMatrix):
    centerMassDistance = __GetCMDists(distanceMatrix)
    matrixA = __buildMatrixA(distanceMatrix,centerMassDistance)

    eigenvalues, eigenvectors = np.linalg.eig(matrixA)

    #The default sort is decreasing by absolute value of eigenvalue
    #Resort so that we are sorting by decreasing real value of the eigenvalue
    ### This is based on the referenced paper, which says we should only use positive valued eigenvalues
    ### Negative eigenvalues would cause issues with the square root later anyway
    sorter = np.flip(np.argsort(eigenvalues))
    eigenvalues = eigenvalues[sorter]
    eigenvectors = eigenvectors[:,sorter]

    #keep only the two largest eigenvalues
    eigenvalues = eigenvalues[:2]
    #keep only the two largest eigenvectors
    eigenvectors = eigenvectors[:,:2]
    if min(eigenvalues) < 0:
        raise Exception('Not enough positive eigenvalues for matrixA!!!')

    bestFitMatrix = np.zeros(eigenvectors.shape)
    bestFitMatrix[:,0] = (eigenvalues[0]**.5)*eigenvectors[:,0]
    bestFitMatrix[:,1] = (eigenvalues[1]**.5)*eigenvectors[:,1]

    return bestFitMatrix

def jsCompatibleConnectivityMatrixFit(connectivityMatrix):
    #connectivity matrix is extressed as a list of lists

    connectivityMatrix = np.array(connectivityMatrix)
    if connectivityMatrix.shape[0] != connectivityMatrix.shape[1]:
        raise exception('The connectivity matrix bust be square!')

    distanceMatrix = buildDistanceMatrix(connectivityMatrix)
    bestFitMatrix = findFitIn2Dimensions(distanceMatrix)

    return bestFitMatrix.tolist()


if __name__ == "__main__":

    # connectivityMatrix = pd.read_csv('test.csv').set_index('index')

    # distanceMatrix = buildDistanceMatrix(connectivityMatrix.values)
    # bestFitMatrix = findFitIn2Dimensions(distanceMatrix)
    # print(bestFitMatrix)

    test = [[0,1,0],
            [1,0,1],
            [0,1,0]]
    print(jsCompatibleConnectivityMatrixFit(test))
            

