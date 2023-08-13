function addParagraph(text) {
    let para = document.createElement("p");
    para.innerHTML = text;

    let element = document.getElementById("main");
    element.appendChild(para);
}

function addTable(tableContents) {
    console.log('here')
    console.log(tableContents)
    console.log(tableContents.length)
    let table = document.createElement("table");
    for (let row_idx = 0; row_idx < tableContents.length; row_idx++) {
        row = table.insertRow(row_idx);
        for (let col_idx = 0; col_idx < tableContents[0].length; col_idx++) {
            col = row.insertCell(col_idx);
            col.innerHTML = tableContents[row_idx][col_idx];
        }
    }

    let element = document.getElementById("main");
    element.appendChild(table);
}

async function main(){
    addParagraph('Starting graph solver. It will take a moment to download and setup the browser version of python...')
    let pyodide = await loadPyodide();
    addParagraph('...Pyodide (broswer-compatible python interpretter) downloaded')

    await pyodide.loadPackage("micropip");
    const micropip = pyodide.pyimport("micropip");
    await micropip.install('numpy');
    addParagraph('...Numpy downloaded and installed in browser')

    pyodide.runPython(`
    import numpy as np
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
        #connectivity matrix is expressed as a list of lists
        
        #remember the "to_py()" to convert the output to a python variable
        connectivityMatrix = np.array(connectivityMatrix.to_py())
        if connectivityMatrix.shape[0] != connectivityMatrix.shape[1]:
            raise exception('The connectivity matrix bust be square!')

        distanceMatrix = buildDistanceMatrix(connectivityMatrix)
        bestFitMatrix = findFitIn2Dimensions(distanceMatrix)

        return bestFitMatrix.tolist()
    `);
    addParagraph('...Custom python functions defined')

    addParagraph('Finding an initial guess for the best 2-dimensional representation of this graph...')
    let connectivityMatrixFit = pyodide.globals.get('jsCompatibleConnectivityMatrixFit');
    let test = [[0,1,0,0,1,0,0,0,1,0],
                [1,0,0,0,0,0,1,0,0,1],
                [0,0,0,0,1,0,1,1,0,0],
                [0,0,0,0,0,1,0,0,1,0],
                [1,0,1,0,0,0,1,0,0,1],
                [0,0,0,1,0,0,0,0,1,1],
                [0,1,1,0,1,0,0,1,0,0],
                [0,0,1,0,0,0,1,0,1,0],
                [1,0,0,1,0,1,0,1,0,1],
                [0,1,0,0,1,1,0,0,1,0]];

    // Remember the "toJs() to convert the output to a javascript variable"
    let result = connectivityMatrixFit(test).toJs()

    addParagraph('...Found initial guess')
    result.splice(0,0,['x','y']) //add "x" and "y" column labels
    addTable(result)
}
main()