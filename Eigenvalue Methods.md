# Math 4610 Fundamentals of Computational Mathematics Software Manual Eigenvalue Methods

**Table of Contents**
1. [largestEigenvallue](#largestEigenvallue)
2. [largestEigenvalueLeslie](#largestEigenvalueLeslie)
3. [smallestEigenvalue](#smallestEigenvalue)
4. [shiftedEigenvalue](#shiftedEigenvalue)
5. [largestTwoEigenvalues](#largestTwoEigenvalues)



<hr>

<a id="largestEigenvallue"></a>

**Routine Name:**         largestEigenvalue  

**Author:** Carter Green

**Language:** Python. The code can be ran using a python interpreter.

**Description/Purpose:** This method finds the largest value of a matrix using the power method. Note that it does assume that the matrix is positive definite.

**Input:** There are 4 inputs. Note that the last 2 have default values if none is given. The first input is the matrix. The second input is a vector that the power method will iterate on. If it is an eigenvector then the power method will fail. The third input is the tolerance for how precise the method is. The default value is 0.00001. The fourth input is the maximum iteration counter. The default value is 30.

**Output:** There is one output: the approximation of the largest eigenvalue

**Usage/Example:**

    x1 = np.array([1,1,3,5])
    A = np.array([[77, 1, 2, 3],
                  [0, 76, 6, 4],
                  [2, 4, 111, 7],
                  [8, 5, 6, 32]])
    print("Largest Eigenvalue", largestEigenvalue(A, x1))



Output from the lines above:

    Largest Eigenvalue 112.48232255224971

  


     

**Implementation/Code:** The following is the code for largestEigenvalue. Assume that the files that are imported from are in the same folder as largestEigenvalue

    #Find Largest Eigenvalue of a matrix
    from twoNormLength import *
    from matrixTimesVector import *
    from dotProduct import *
    import numpy as np
    
    
    def largestEigenvalue(A, x0, tol = .00001, maxiter = 30):
    
        #initialize
        error = 10*tol
        x0 = np.array(x0)
        length = twoNormLength(x0)
        x0 = 1/length * x0
        y = np.array(matrixTimesVector(A, x0))
        lambda0 = dotProduct(x0, y)
        iter = 0
    
        #loop
        while(error > tol and iter < maxiter):
            x1 = 1/twoNormLength(y) * y
            y = np.array(matrixTimesVector(A, x1))
            lambda1 = np.dot(x1, y)
            error = np.absolute(lambda1 - lambda0)
            iter = iter +1
            lambda0 = lambda1
                   
        
        return lambda1

     

**Last Modified:** November/2023






<hr>

<a id="largestEigenvalueLeslie"></a>

**Routine Name:**          largestEigenvalueLeslie 

**Author:** Carter Green

**Language:** Python. The code can be ran using a python interpreter.



**Description/Purpose:** This method finds an approximation of the largest eigenvalue of a Leslie matrix using the power method.

**Input:** There are 4 inputs. Note that the last 2 have default values if none is given. The first input is the matrix. The second input is a vector that the power method will iterate on. If it is an eigenvector then the power method will fail. The third input is the tolerance for how precise the method is. The default value is 0.00001. The fourth input is the maximum iteration counter. The default value is 30.

**Output:** There is one output: the approximation of the largest eigenvalue


**Usage/Example:**

    # Leslie Matrix
    A = np.array(
        [[1,1,2,1],
         [.1,0,0,0],
         [0,.05,0,0],
         [0,0,.1,0]])
    x1 = np.array(makeVectorOne(4))
    print(largestEigenvalueLeslie(A, x1))



Output from the lines above:

    1.0995897202241027


     

**Implementation/Code:** The following is the code for largestEigenvalueLeslie. Assume that the files that are imported from are in the same folder as largestEigenvalueLeslie.

    #largest eigenvalue of a Leslie Matrix
    
    from twoNormLength import *
    from matrixTimesVector import *
    from dotProduct import *
    import numpy as np
    
    def largestEigenvalueLeslie(A, x0, tol = .00001, maxiter = 30):
    
        #initialize
        error = 10*tol
        x0 = np.array(x0)
        length = twoNormLength(x0)
        x0 = 1/length * x0
        
        y = []
        sum = 0
        for i in range(len(x0)):
            sum = sum + A[0][i] * x0[i]
        y.append(sum)
        for i in range(len(x0)-1):
            y.append(A[i+1][i] * x0[i])
        y = np.array(y)
        
        lambda0 = dotProduct(x0, y)
        iter = 0
    
    
        #loop
        while(error > tol and iter < maxiter):
            x1 = 1/twoNormLength(y) * y
            
            y = []
            sum = 0
            for i in range(len(x0)):
                sum = sum + A[0][i] * x1[i]
            y.append(sum)
            for i in range(len(x0)-1):
                y.append(A[i+1][i] * x1[i])
            y = np.array(y)
    
                     
            lambda1 = dotProduct(x1, y)
            error = abs(lambda1 - lambda0)
            iter = iter +1
            lambda0 = lambda1
    
        return lambda1

     

**Last Modified:** November/2023





<hr>

<a id="smallestEigenvalue"></a>

**Routine Name:**      smallestEigenvalue     

**Author:** Carter Green

**Language:** Python. The code can be ran using a python interpreter.

**Description/Purpose:** This method finds the smallest value of a matrix using the inverse power method. The method assumes the matrix is positive definite.

**Input:** There are 4 inputs. Note that the last 2 have default values if none is given. The first input is the matrix. The second input is a vector that the inverse power method will iterate on. If it is an eigenvector then the inverse power method will fail. The third input is the tolerance for how precise the method is. The default value is 0.0000000001. The fourth input is the maximum iteration counter. The default value is 30.


**Output:** There is one output - the smallest eigenvalue of the matrix.

**Usage/Example:**

    x1 = np.array([1,1,3,5])
    A = np.array([[77, 1, 2, 3],
                  [0, 76, 6, 4],
                  [2, 4, 111, 7],
                  [8, 5, 6, 32]])
    print("Smallest Eigenvalue", smallestEigenvalue(A, x1, .0000001))




Output from the lines above:

    Smallest Eigenvalue 30.64630336852988


     

**Implementation/Code:** The following is the code for smallestEigenvalue. Assume that the files that are imported from are in the same folder as smallestEigenvalue.


    #Inverse Power Method
    from twoNormLength import *
    from matrixTimesVector import *
    from dotProduct import *
    from LUfactorization import *
    from forwardSubDiagEquals1 import *
    from backSubstitution import *
    from copyMatrix import *
    import numpy as np
    
    
    
    def smallestEigenvalue(A0, x0, tol = .0000000001, maxiter = 50):
    
        #initialize
        A = copyMatrix(A0)
        error = 10*tol
        x0 = np.array(x0)
        length = twoNormLength(x0)
        x0 = 1/length * x0
        LUfactorization(A)
        y2 = forwardSubDiagEquals1(A, x0)
        y = backSubstitution(A, y2)
        y = np.array(y)
        lambda0 = dotProduct(x0, y)
        iter = 0
    
        #loop
        while(error > tol and iter < maxiter):
            x1 = 1/twoNormLength(y) * y
            y2 = forwardSubDiagEquals1(A, x1)
            y = backSubstitution(A, y2)
            y = np.array(y)
            lambda1 = dotProduct(x1, y)
            error = abs(lambda1 - lambda0)
            iter = iter +1
            lambda0 = lambda1
              
        return 1/lambda1

     

**Last Modified:** November/2023




<hr>

<a id="shiftedEigenvalue"></a>

**Routine Name:**           shiftedEigenvalue

**Author:** Carter Green

**Language:** Python. The code can be ran using a python interpreter.

**Description/Purpose:** This method finds the smallest eigenvalue of A-rI, where A is a matrix and r is a shift to that matrix. It can be used to find eigenvalues of A other than the largest or smallest ones.

**Input:** There are 5 inputs. Note that the last 2 have default values if none is given. The first input is the matrix. The second input is the number that the matrix is being shifted by. The third input is a vector that the inverse power method will iterate on. If it is an eigenvector then the inverse power method will fail. The fourth input is the tolerance for how precise the method is. The default value is 0.0000000001. The fifth input is the maximum iteration counter. The default value is 30.

**Output:** There is one output: the smallest eigenvalue of A-rI.

**Usage/Example:**

    x1 = np.array([1,1,3,5])
    A = np.array([[77, 1, 2, 3],
                  [0, 76, 6, 4],
                  [2, 4, 111, 7],
                  [8, 5, 6, 32]])
    print("A shifted Eigenvalue: ", shiftedEigenvalue(A, 70.0, x1, .0000001))



Output from the lines above:

  
    A shifted Eigenvalue:  75.51431157760199


     

**Implementation/Code:** The following is the code for shiftedEigenvalue. Assume that the files that are imported from are in the same folder as shiftedEigenvalue.


    #Find Smallest eigenvalue of a shifted matrix, A-If
    
    from twoNormLength import *
    from matrixTimesVector import *
    from dotProduct import *
    from LUfactorization import *
    from forwardSubDiagEquals1 import *
    from backSubstitution import *
    from copyMatrix import *
    import numpy as np
    
    
    def shiftedEigenvalue(A0, shift, x0, tol = .0000000001, maxiter = 30):
    
        #initialize
        A = copyMatrix(A0)
        for i in range(len(A)):
            A[i][i] = A[i][i] - shift
        A = np.array(A)
        error = 10*tol
        x0 = np.array(x0)
        length = twoNormLength(x0)
        x0 = 1/length * x0
        LUfactorization(A)
        y2 = forwardSubDiagEquals1(A, x0)
        y = backSubstitution(A, y2)
        y = np.array(y)
        lambda0 = dotProduct(x0, y)
        iter = 0
    
        #loop
        while(error > tol and iter < maxiter):
            x1 = 1/twoNormLength(y) * y
            y2 = forwardSubDiagEquals1(A, x1)
            y = backSubstitution(A, y2)
            y = np.array(y)
            lambda1 = dotProduct(x1, y)
            error = abs(lambda1 - lambda0)
            iter = iter +1
            lambda0 = lambda1
              
        return 1/lambda1 + shift

     

**Last Modified:** November/2023




<hr>

<a id="largestTwoEigenvalues"></a>

**Routine Name:**           largestTwoEigenvalues

**Author:** Carter Green

**Language:** Python. The code can be ran using a python interpreter.



**Description/Purpose:** This method finds the largest 2 eigenvalues of a matrix using the power method and the shifted inverse power method.

**Input:** There are 4 inputs. Note that the last 2 have default values if none is given. The first input is the matrix. The second input is a vector that the inverse power method will iterate on. If it is an eigenvector then the inverse power method will fail. The third input is the tolerance for how close the computed eigenvalue must be to a previous computation of that same eigenvalue. The default value is 0.1. The fourth input is the maximum iteration counter. The default value is 30.

**Output:** There are two outputs, placed into one array. The first output is the second largest eigenvalue. The second output is the largest eigenvalue.

**Usage/Example:**

    A = np.array([[77, 1, 2, 3],
                  [0, 76, 6, 4],
                  [2, 4, 111, 7],
                  [8, 5, 6, 32]])
    x1 = np.array([1,1,3,5])
    print("The largest 2 eigenvalues are ", largestTwoEigenvalues(A, x1))




Output from the lines above:

    The largest 2 eigenvalues are  [77.3673525801705, 112.48232255224971]
  
     

**Implementation/Code:** The following is the code for largestTwoEigenvalues. Assume that the files that are imported from are in the same folder as largestTwoEigenvalues

    
    #find Largest 2 eigenvalues
    
    from largestEigenvalue import *
    from twoNormDistance import *
    from smallestEigenvalue import *
    from shiftedEigenvalue import *
    import numpy as np
    
    
    def largestTwoEigenvalues(A0, x0, tol = .1, maxiter = 30):
        #find largest and smallest eigenvalue
        A = copyMatrix(A0)
        largestEig = largestEigenvalue(A, x0)
        A = copyMatrix(A0)
        smallestEig = smallestEigenvalue(A, x0, .0000001)
        ans = []
        ans.append(1)
        ans.append(largestEig)
        iter = 0
    
        #loop until it is known the largest 2 are found
        found = False
        a = smallestEig
        b = largestEig
        lastEig = smallestEig
        while(not found and iter < maxiter):
            
            #check the midpoint of the current largest and second largest eigenvalue
            #if a new eigenvalue is not gotten, then we have the largest & second largest
            c = (a+b)/2
            A = copyMatrix(A0)
            shiftedEig = shiftedEigenvalue(A, c, x0, .0000000001)
            if(np.absolute(shiftedEig -  a) < tol or np.absolute(shiftedEig - b) < tol):
                ans[0] = lastEig
                found = True
            else:
                a = c
                lastEig = shiftedEig
                iter = iter + 1
                
        #return the array with the two eigenvalues
        return ans
                
                
        

     

**Last Modified:** November/2023


