# Math 4610 Software Manual for Every Method

**Table of Contents**
1. [largestEigenvallue](#largestEigenvallue)
2. [largestEigenvalueLeslie](#largestEigenvalueLeslie)
3. [smallestEigenvalue](#smallestEigenvalue)
4. [shiftedEigenvalue](#shiftedEigenvalue)
5. [largestTwoEigenvalues](#largestTwoEigenvalues)
7. [jacobiIteration](#jacobiIteration)
8. [jacobiLeslie](#jacobiLeslie)
9. [makeLeslie](#makeLeslie)
10. [gaussSeidelIteration](#gaussSeidelIteration)
11. [testCode](#testCode)
12. [BisectionMethod](#BisectionMethod)
13. [Newton'sMethod](#Newton'sMethod)
14. [SecantMethod](#SecantMethod)
15. [HybridMethod](#HybridMethod)
16. [HybridMethodKnownRoot](#HybridMethodKnownRoot)
17. [LogisticODERoots](#LogisticODERoots)
18. [TestRoutine](#TestRoutine)
19. [upperTriangular](#upperTriangular)
20. [backSubstitution](#backSubstitution)
21. [LUfactorization](#LUfactorization)
22. [forwardSubstitution](#forwardSubstitution)
23. [forwardSubDiagEquals1](#forwardSubDiagEquals1)
24. [matrixTimesVector](#matrixTimesVector)
25. [makeMatrix](#makeMatrix)
26. [makeVectorOne](#makeVectorOne)
28. [TestingOfMatrixOperations](#TestingOfMatrixOperations)
29. [smaceps](#smaceps)
30. [dmaceps](#dmaceps)
31. [l2-norm](#l2-norm)
32. [l1-norm](#l1-norm)
33. [linf-norm](#linf-norm)
34. [l2-distance](#l2-distance)
35. [linf-distance](#linf-distance)
36. [linreg](#linreg)
37. [fordiff](#fordiff)
38. [backdiff](#backdiff)
39. [upperTriangular](#upperTriangular1)
40. [backSubstitution](#backSubstitution1)
41. [jacobiOpenMP](#jacobi)
42. [power_methodOpenMP](#power_method)



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

<hr>

<a id="jacobiIteration"></a>

**Routine Name:**  jacobiIteration        

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This code gets an approximation for the solution of a linear system of equations. This method will always work if a system is diagonally dominant. Otherwise, it might work or it might diverge.

**Input:** This code has 5 inputs. The first one is the matrix that we will iteratite on. The second input is the solution vector we are solving for. The third input is the guess of the vector we are solving for. The fourth input is the minimum error and it has a default value of .00001. The fifth input is the maximum iteration count and it has a default value of 70.

**Output:** This code has 2 outputs in duple format. The first output is the vector x that solves the system Ax = y. The second output is the number of FLOPS that it took to compute that vector.

**Usage/Example:** 

    #test Jacobi Iteration
    A = makeMatrix(100)
    x = makeVector(100)
    y = matrixTimesVector(A, x)
    x0 = makeVectorOne(100)
    x1, counterJ = jacobiIteration(A, y, x0)
    print("Error for Jacobi: ", twoNormDistance(x, x1))



Output from the lines above:

    Error for Jacobi:  2.849335286387407e-08
  


     

**Implementation/Code:** The following is the code for jacobiIteration. Assume that all the files imported are in the same folder as the jacobiIteration method.

    from matrixTimesVector import *
    from twoNormLength import *
    #solve with Jacobi Iteration
    
    def jacobiIteration(A, b, x0, tol = .00001, maxiter = 70):
        #Initialize
        counter = 0
        error = 10*tol
        iter = 0
    
        #Loop
        while(error > tol and maxiter > iter):
            #calculate residual
            r = []
            Ax = matrixTimesVector(A, x0)
            counter += (2*len(x0) - 1) * len(x0)
            for i in range(len(b)):
                r.append(b[i] - Ax[i])
                counter += 1
    
            #calculate D^{-1}(r_k) + x_k
            x1 = []
            for i in range(len(b)):
                x1.append(1/A[i][i] * r[i] + x0[i])
                counter += 2
    
            #calculate error
            error = twoNormLength(r)
            counter += 2*len(x0) + 1
            
            iter += 1
            x0 = x1
    
    
        #print out number of FLOPS and return
        #print("Number of FLOPS for Jacobi: ",  counter)
        return x1, counter
        



   

**Last Modified:** December/2023



<hr>

<a id="jacobiLeslie"></a>

**Routine Name:**     jacobiLeslie     

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This method computes an approximation of x in the linear system Ax=b, when A is a Leslie Matrix. 

**Input:** This code has 5 inputs. The first one is the matrix that we will iteratite on. The second input is the solution vector we are solving for. The third input is the guess of the vector we are solving for. The fourth input is the minimum error and it has a default value of .00001. The fifth input is the maximum iteration count and it has a default value of 70.

**Output:** This code has 2 outputs in duple format. The first output is the vector x that solves the system Ax = y. The second output is the number of FLOPS that it took to compute that vector.

**Usage/Example:**

    #test Jacobi on a Leslie Matrix
    A1 = makeLeslie(5)
    x = makeVector(5)
    y = matrixTimesVector(A1, x)
    x0 = makeVectorOne(5)
    x1 = jacobiLeslie(A1, y, x0)
    print("Error for Jacobi Iteration on a Leslie Matrix: ", twoNormDistance(x, x1))



Output from the lines above:

    Error for Jacobi Iteration on a Leslie Matrix:  1.5700924586837752e-16
     

**Implementation/Code:** The following is the code for jacobiLeslie. Assume that all the imported files are in the same folder as the method.

    from matrixTimesVector import *
    from twoNormLength import *
    
    #solve a Leslie matrix with Jacobi Iteration
    def jacobiLeslie(A, b, x0, tol = .00001, maxiter = 70):
        #Initialize
        error = 10*tol
        iter = 0
        
        #permute rows
        a = []
        B = []
        for i in range(len(A) -1):
            a.append(A[i+1])
            B.append(b[i+1])
        a.append(A[0])
        B.append(b[0])
        A = a
        b = B
    
        #Loop
        while(error > tol and maxiter > iter):
            #calculate residual
            r = []
            Ax = matrixTimesVector(A, x0)
            for i in range(len(b)):
                r.append(b[i] - Ax[i])
    
            #calculate D^{-1}(r_k) + x_k
            x1 = []
            for i in range(len(b)):
                x1.append(1/A[i][i] * r[i] + x0[i])
    
            #calculate error
            error = twoNormLength(r)
            iter += 1
    
            x0 = x1
    
    
    
        return x1


**Last Modified:** December/2023



<hr>

<a id="makeLeslie"></a>

**Routine Name:**         makeLeslie 

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This method makes a Leslie Matrix

**Input:** This method has 1 input: the size of the Leslie Matrix

**Output:** This method will output a pseudo-randomly generated Leslie Matrix

**Usage/Example:**

    A1 = makeLeslie(5)



Output from the lines above:

There won't be a visible output from running that code, but A1 will have a Leslie Matrix stored in it.


     

**Implementation/Code:** The following is the code for makeLeslie.

    #make a Leslie Matrix
    def makeLeslie(n):
        A = []
        #make top row
        row = []
        for i in range(n):
            row.append(2*random.random()+n)
        A.append(row)
        #make the rest
        for i in range(n-1):
            rows = []
            for k in range(i):
                rows.append(0)
            rows.append(random.random())
            for j in range(n-1-i):
                rows.append(0)
            A.append(rows)
        #return matrix
        return A
   

**Last Modified:** December/2023




<hr>

<a id="gaussSeidelIteration"></a>

**Routine Name:**          gaussSeidelIteration

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This method performs Gauss-Seidel iteration to find an approximate solution for the x vector that solves Ax = b.

**Input:** This code has 5 inputs. The first one is the matrix that we will iteratite on. The second input is the solution vector we are solving for. The third input is the guess of the vector we are solving for. The fourth input is the minimum error and it has a default value of .00001. The fifth input is the maximum iteration count and it has a default value of 70.


**Output:** This code has 2 outputs in duple format. The first output is the vector x that solves the system Ax = y. The second output is the number of FLOPS that it took to compute that vector.


**Usage/Example:**

    A = makeMatrix(100)
    x = makeVector(100)
    y = matrixTimesVector(A, x)
    x0 = makeVectorOne(100)
    x1, counterG = gaussSeidelIteration(A, y, x0)
    print("Error for Gauss-Seidel: ", twoNormDistance(x, x1))



Output from the lines above:

    Error for Gauss-Seidel:  2.827835298297449e-09


     

**Implementation/Code:** The following is the code for gaussSeidelIteration. Assume all imported files are in the same folder as the method.

    from matrixTimesVector import *
    from twoNormLength import *
    from backSubstitution import *
    
    #Gauss Seidel Iteration
    def gaussSeidelIteration(A, b, x0, tol = .00001, maxiter = 70):
        #Initialize
        error = 10*tol
        iter = 0
        counter = 0
    
        #Loop
        while(error > tol and maxiter > iter):
            #calculate residual
            r = []
            Ax = matrixTimesVector(A, x0)
            counter += (2*len(x0) - 1) * len(x0)
            for i in range(len(b)):
                r.append(b[i] - Ax[i])
                counter += 1
    
            #calculate x1
            DPlusUInverse = backSubstitution(A, r)
            counter += (2*len(x0) - 1) * len(x0)
            x1 = []
            for i in range(len(b)):
                x1.append(DPlusUInverse[i] + x0[i])
                counter +=1
    
            #calculate error
            error = twoNormLength(r)
            counter += 2*len(x0) + 1
            iter += 1
    
            x0 = x1
    
        #print out number of FLOPS and return
        #print("Number of FLOPS for Gauss-Seidel: ",  counter)
        return x1, counter
    
            
            


   

**Last Modified:** December/2023


<hr>

<a id="testCode"></a>

**Routine Name:**    testCode    

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This method tests the methods explained in this file and then compares them to an LU-factorization method of solving. Make sure that if running this code you have downloaded the tabulate methods from the internet.

**Input:** This method has no inputs.

**Output:** This method has no outputs.

**Usage/Example:**

    main()



Output from the lines above:

    Error for Jacobi:  2.641006762881539e-08
    Error for Jacobi Iteration on a Leslie Matrix:  1.5700924586837752e-16
    Error for Gauss-Seidel:  1.8219216975537628e-09
    Error for LU-Factorization:  5.674097786814442e-15
      Jacobi error    Jacobi FLOPS    G.S. error    G.S. FLOPS     LU error    LU FLOPS
    --------------  --------------  ------------  ------------  -----------  ----------
       2.54299e-08          550827   1.91272e-09        402010  3.18877e-15      711350
       3.00448e-08          550827   2.41212e-09        402010  2.95938e-15      711350
       2.91761e-08          550827   2.59832e-09        402010  2.99736e-15      711350
       3.14929e-08          550827   2.10264e-09        402010  2.8011e-15       711350
       2.29135e-08          550827   1.89959e-09        402010  2.6401e-15       711350
       3.28202e-08          550827   2.48891e-09        402010  2.9154e-15       711350
       2.92467e-08          550827   1.98246e-09        402010  3.13384e-15      711350
       2.83322e-08          550827   2.26965e-09        402010  2.7862e-15       711350
       1.70223e-08          571228   2.61696e-09        402010  3.27929e-15      711350
       1.99432e-08          550827   1.72574e-09        402010  3.27535e-15      711350
       2.64453e-08          550827   1.95739e-09        402010  2.97894e-15      711350
       3.06034e-08          550827   2.5241e-09         402010  2.48905e-15      711350
       2.44806e-08          550827   2.01409e-09        402010  2.91769e-15      711350
       2.7189e-08           550827   2.56603e-09        402010  2.47772e-15      711350
       2.38508e-08          550827   2.21612e-09        402010  2.60713e-15      711350
       1.8408e-08           571228   2.48226e-09        402010  2.6364e-15       711350
       3.14352e-08          550827   2.47517e-09        402010  3.08188e-15      711350
       2.7744e-08           550827   2.4814e-09         402010  3.13492e-15      711350
       2.99935e-08          550827   1.7817e-09         402010  3.16867e-15      711350
       2.97153e-08          550827   1.90869e-09        402010  3.09439e-15      711350
    

  


     

**Implementation/Code:** The following is the code for testCode. Assume all files imported are in the same folder as the method.

        from jacobiIteration import *
        from matrixTimesVector import *
        from twoNormDistance import *
        from makeMatrix import *
        from jacobiLeslie import *
        from gaussSeidelIteration import *
        from LUfactorization import *
        from forwardSubDiagEquals1 import *
        from backSubstitution import *
        from tabulate import tabulate
        
        def LUFactor(A, y, x0):
            counter = 0
            counter += LUfactorization(A)
            
            y2 = forwardSubDiagEquals1(A, y)
            counter += (2*len(x0) - 1) * len(x0)
            x2 = backSubstitution(A, y2)
            counter += (2*len(x0) - 1) * len(x0)
            return x2, counter
        
        #test code for LinearSystems
        def main():
            #test Jacobi Iteration
            A = makeMatrix(100)
            x = makeVector(100)
            y = matrixTimesVector(A, x)
            x0 = makeVectorOne(100)
            x1, counterJ = jacobiIteration(A, y, x0)
            print("Error for Jacobi: ", twoNormDistance(x, x1))
        
            #test Jacobi on a Leslie Matrix
            A1 = makeLeslie(5)
            x = makeVector(5)
            y = matrixTimesVector(A1, x)
            x0 = makeVectorOne(5)
            x1 = jacobiLeslie(A1, y, x0)
            print("Error for Jacobi Iteration on a Leslie Matrix: ", twoNormDistance(x, x1))
        
            #test Gauss-Seidel
            A = makeMatrix(100)
            x = makeVector(100)
            y = matrixTimesVector(A, x)
            x0 = makeVectorOne(100)
            x1, counterG = gaussSeidelIteration(A, y, x0)
            print("Error for Gauss-Seidel: ", twoNormDistance(x, x1))
        
            #test LU Factorization
            n = 100
            A = makeMatrix(n)
            x = makeVectorOne(n)
            y = matrixTimesVector(A,x)
            x1L, counterL = LUFactor(A, y, x0)
            print("Error for LU-Factorization: ", twoNormDistance(x, x1L))
        
        
            #compare the 3 methods in a table
            data = []
            for i in range(20):
                data.append([])
                A = makeMatrix(100)
                x = makeVector(100)
                y = matrixTimesVector(A, x)
                x0 = makeVectorOne(100)
        
                #Jacobi
                x1J, counterJ = jacobiIteration(A, y, x0)
                data[i].append(twoNormDistance(x,x1J))
                data[i].append(counterJ)
        
                #Gauss-Seidel
                x1G, counterG = gaussSeidelIteration(A, y, x0)
                data[i].append(twoNormDistance(x,x1G))
                data[i].append(counterG)
        
                #LU Factorazation
                x1L, counterL = LUFactor(A, y, x0)
                data[i].append(twoNormDistance(x,x1L))
                data[i].append(counterL)
                
            print(tabulate( data, headers=["Jacobi error", "Jacobi FLOPS", "G.S. error", "G.S. FLOPS", "LU error", "LU FLOPS"] ))
        
        main()



**Last Modified:** December/2023

<hr>
<a id="BisectionMethod"></a>

**Routine Name:**    Bisection Method

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** Given two endpoints with opposite signs of function values (one is positive and one is negative), this method will solve for an x-value where the function equals zero. This does have the requirements that the function is continuous over the interval given.

**Input:** This method has 4 inputs. The first one is the function (defined in its own method) that the root is being solved for. The next one is a, the start point of the interval. The third input is b, the end point of the interval. The last input is the tolerance, or how precise, we want the method to get.

**Output:** This method has 2 possible outputs. If the endpoints given have the same sign on the function values, the method will return NULL and print an error message to the screen. If the endpoints given are valid (and the function is continuous over the interval), then the approximate x-value where the function has a zero will be returned.

**Usage/Example:**
The following is a snippet from a test code: Assume that the test file and the Bisection method file are in the same folder so that the import works.

    from BisectionMethod import *
    def function1(x):
      return x*x - 5*x+6
    val = bisection(function1, 1.5, 2.6, .01)
    print("Testing Bijection Method", val)



Output from the lines above:

    Testing Bijection Method 2.00703125
  


     

**Implementation/Code:** The following is the code for the Bisection Method:

    import math

    #Bisection method
    def bisection(function, a, b, tol):
        #initialize variables
        fa = function(a)
        fb = function(b)
        if(fa == 0): return a
        if(fb == 0): return b
        if(fa * fb > 0):
            print("Bisectino method failed")
            return;
        n = (int) (math.log(tol/(b-a))/math.log(1/2))+1
    
        #loop
        for i in range(n):
            c = 1/2*(a+b)
            fc = function(c)
            if(fa*fc < 0):
               b=c
               fb = fc
            elif (fc == 0): return c
            else:
                a=c
                fa = fc
            i = i + 1
    
    
        #return final value
        return c

   

**Last Modified:** November/2023



<hr>

<a id="Newton'sMethod"></a>

**Routine Name:**          Newton's Method

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This routine returns the x-value approximately at which a function has a zero. AKA it solves for a root. In order for the routine to work, the function must be twice continuous and differentiable, f'(x) cannot equal zero, and the given guess of the root must be close to the root in question.

**Input:** This routine has 5 inputs. The first one is a pointer to the function that must be defined in another routine. The second input is the pointer to the derivative of the function. This must also be defined in another routine. The third input is the guess for the root. The fourth input is the tolerance, or how much error we are allowing the routine to have. The last input is the maximum number of iterations that we will let the routine run before stopping it.

**Output:** If the conditions are met, this routine will return an approximation for the root. If the conditions aren't met, the routine may crash with an error message or return a garbage value of an approximation.

**Usage/Example:**
The following is a snippet from a test code: Assume that the test file and the Newton's method file are in the same folder so that the inport works.

    from NewtonsMethod import *
    def function1(x):
        return x*x - 5*x+6

    def derivativeFunction1(x):
        return 2*x-5
        
    val = newtonsMethod(function1, derivativeFunction1, 4, .01, 30)
    print("Testing Newton's Method", val)
    

Output from the lines above:

  
    Testing Newton's Method 3.0000152590218967

     

**Implementation/Code:** The following is the code for Newton's Method

    #Newton's Method for root finding
    def newtonsMethod(function, derivative, x0, tol, maxiter):
        error = 10*tol
    
        #loop
        i=0
        while (i < maxiter and error > tol):
            x1 = x0 - function(x0) / derivative(x0)
            error = abs(x1-x0)
            x0 = x1
            i = i + 1
            
        return x0


   

**Last Modified:** November/2023




<hr>

<a id="SecantMethod"></a>

**Routine Name:**          Secant Method

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This routine finds and returns an approximation for the root of a function. In order for this routine to work, f'(x) cannot equal zero, f (the function) must be twice continuous and differentiable, and the two guesses for the root must be different and sufficiently close to the actual root

**Input:** This routine has five inputs. The first one is a pointer to the routine that defines the function. The second input is the first guess for the root. The third input is the second guess for the root. The fourth input is the tolerance for how exact the routine will get. The fifth input is the maximum number of iterations through the method the function will go through before stopping.

**Output:** If the conditions are met, this routine will return an approximation for a root of a function. If the conditions are not met, then either the program will crash or a garbage approximation will be returned

**Usage/Example:**
The following is a snippet from a test code: Assume that the test file and the secant method file are in the same folder so that the inport works.

    from SecantMethod import *
    def function1(x):
        return x*x - 5*x+6
    val = secantMethod(function1, .01, 0, .01, 30)
    print("Testing Secant Method", val)


Output from the lines above:

    Testing Secant Method 1.9949042720195709

     

**Implementation/Code:** The following is the code for Secant Method

    #Secant Method for finding roots
    def secantMethod(function, x1, x0, tol, maxiter):
        error = 10*tol
    
        #loop
        i=0
        while (i < maxiter and error > tol):
            x2 = x0 - function(x0) * (x1-x0) / (function(x1) - function(x0)) 
            error = abs(x2-x1)
            x0 = x1
            x1 = x2
            i = i + 1
            
        return x0

  

**Last Modified:** November/2023




<hr>

<a id="HybridMethod"></a>

**Routine Name:**         Hybrid Method 

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This routine will return an approximation for the root of a function. Because it utilizes both the bisection method and the secant method, the conditions for both of those must be met. i.e. the function must be twice continuous and differentiable, the derivative cannot equal zero, and the function values of initial guesses for the root must have opposite signs.

**Input:** This routine has 5 inputs. The first one is a pointer to the routine that contains the function. The second one is the first guess for the root. The third input is the second guess for the root. The fourth input is the tolerance for the error of the routine. The fifth input is the maximum number of iterations the routine will go through before stopping.

**Output:** If the conditions are met, then this routine will return an approximation for the root of a function. If the conditions are not met, then either the program will crash or an error message will be printed or a garbage approximation will be returned.

**Usage/Example:**
The following is a snippet from a test code: Assume that the test file and the hybrid method file are in the same folder so that the import works.

    from HybridMethod import *
    def function1(x):
        return x*x - 5*x+6
    val = hybridMethod(function1, 2.4, 1111, .01, 30)
    print("Testing a hybrid Method", val)


Output from the lines above:

    Testing a hybrid Method 3.0001240687574207

  


     

**Implementation/Code:** The following is the code for the hybrid method

    #Hybrid method for root-finding compining bisection and secant
    def hybridMethod(function, a, b, tol, maxiter):
    
        #initialize
        fa = function(a)
        fb = function(b)
        if(fa*fb >= 0):
            print("Hybrid Method Failed")
            return
        error = 10*tol
        iter = 0
        #to give some values to x0 and x1 that are within a and b
        x0 = 1/3*(a+b)
        x1 = 1/2*(a+b)
    
        #loop
        while(error > tol and iter < maxiter):
            x2 = x0 - function(x0) * (x1-x0) / (function(x1) - function(x0))
    
            #if x2 < a or >b
            if(abs(x2-x1) > .5*(b-a)):
                iter = iter + 1
                c=x1
                fc = function(c)
    
                #reduce the range by an order of magnitude
                for i in range(4):
                    if(fa*fc < 0):
                        b=c
                        fb=fc
                    else:
                        a=c
                        fa=fc
                    c = .5*(a+b)
                    fc = function(c)
                error = .5*(a+b)
                x0 = x1
                x1 = .5*(a+b)
                    
            #else secant method
            else:
                while(iter < maxiter and error > tol):
                    error = abs(x2 - x1)
                    x0 = x1
                    x1 = x2
                    x2 = x0 - function(x0) * (x1-x0) / (function(x1) - function(x0))
                    iter = iter + 1
    
        return x1
            
            

**Last Modified:** November/2023





<hr>

<a id="HybridMethodKnownRoot"></a>

**Routine Name:**          Hybrid Method Known Root

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This routine returns a root. In this routine, it is assumed that there is a root on the given interval. So, the routine searches for an interval on which the bisection method can be applied and then uses a hybrid method involving the bisection method and the secant method. Because it utilizes the secant method, the conditions for that must be met. i.e. the function must be twice continuous and differentiable and the derivative cannot equal zero. Additionally, there must be a root on the interval given.

**Input:** This routine has 5 inputs. The first one is a pointer to the routine that contains the function. The second one is the first guess for the root. The third input is the second guess for the root. The fourth input is the tolerance for the error of the routine. The fifth input is the maximum number of iterations the routine will go through before stopping. 

**Output:** If the conditions are met, then this routine will return an approximation for the root of a function. If the conditions are not met, then either the program will crash or an error message will be printed or a garbage approximation will be returned.

**Usage/Example:**
The following is a snippet from a test code: Assume that the test file and the Hybrid Method Known Root file (in this case named HybridMethodAlt) are in the same folder so that the import works.

    from HybridMethodAlt import *
    def function1(x):
        return x*x - 5*x+6
    val = hybridMethodAlt(function1, -200, 3001, .01, 30)
    print("Testing alternate hybrid Method", val)
      

Output from the lines above:


    Testing alternate hybrid Method 1.999990234293511
        

**Implementation/Code:** The following is the code for Hybrid Method Known Root (names hybridMethodAlt in the code)

    #Hybrid method for root-finding compining bisection and secant
    #Assume there is a root in the interval given
    def hybridMethodAlt(function, a, b, tol, maxiter):
    
        #initialize
        fa = function(a)
        fb = function(b)
        foundInterval = False
        numIntervals = 1
    
        #find an interval that the bisection method will work
        #start with 1 interval. If not work then double again and again
            #until interval is found
        while (not foundInterval):
            for i in range(1, numIntervals+1):
    
                #only check odd numbered interval borders, because even ones
                    #already checked
                if(i%2 == 1 and not foundInterval):
                    c = i/numIntervals*(b-a)+a
                    fc = function(c)
                    if(fa * fc < 0):
                        
                        b = c
                        fb = fc
                        foundInterval = True               
    
            #update numIntervals
            numIntervals = numIntervals * 2
        
        error = 10*tol
        iter = 0
        #to give some values to x0 and x1 that are within a and b
        x0 = 1/3*(a+b)
        x1 = 1/2*(a+b)
    
        #loop
        while(error > tol and iter < maxiter):
            x2 = x0 - function(x0) * (x1-x0) / (function(x1) - function(x0))
    
            #if x2 < a or >b
            if(abs(x2-x1) > .5*(b-a)):
                iter = iter + 1
                c=x1
                fc = function(c)
    
                #reduce the range by an order of magnitude
                for i in range(4):
                    if(fa*fc < 0):
                        b=c
                        fb=fc
                    else:
                        a=c
                        fa=fc
                    c = .5*(a+b)
                    fc = function(c)
                error = .5*(a+b)
                x0 = x1
                x1 = .5*(a+b)
                    
            #else secant method
            else:
                while(iter < maxiter and error > tol):
                    error = abs(x2 - x1)
                    x0 = x1
                    x1 = x2
                    x2 = x0 - function(x0) * (x1-x0) / (function(x1) - function(x0))
                    iter = iter + 1
                    
        return x1
    
    #def function(x):
        #replace this with your function
   

**Last Modified:** November/2023





<hr>

<a id="LogisticODERoots"></a>

**Routine Name:**       LogisticODERoots   

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This routine solves for the roots of a logistic equation minus some value, if it exists. If it doesn't exist, then it will print out that the logistic equation never reaches that value. Because this routine is a hybrid method combining the bisection method and the secant method, the conditions for both of those methods must be met, which they will because a logistic equation is twice continuous and differentiable. Note that the logistic equation is the solution to y' = \alpha * y - \beta * y^2

**Input:** This routine has 4 inputs. The first one is the value for \alpha. The second one is the value for \beta. The third input is P0, or the initial condition for the logistic equation. The fourth input is the value called Pinfty, which is the one we are solving whether or not the logistic equation reaches.

**Output:** The routine will output one of two cases. If the logistic equation doesn't ever reach the given Pinfty value, then the routine will output a message saying stuff. If the logistic equation does reach the given Pinfty value, then the routine will output the x-value where it reaches it.

**Usage/Example:**
The following is a snippet from a test code: Assume that the test file and the logistic ODE file are in the same folder so that the import works.

    from LogisticODERoots import *
    val1 = logisticODERoots(.1, .001,2,29.75)
    print(val1)
    val2 = logisticODERoots(.1, .001, 2, 115.35)
    print(val2)


Output from the lines above:

    30.32588873949352
    The population density of the pest species will never reach Pinfty
    None


     

**Implementation/Code:** The following is the code for LogisticODEroots

    import math
    #solve Logistic equation
    #\(\frac{dP}{dt} = \alpha\ P - \beta\ P\)
    def logisticODERoots(alpha, beta, P0, Pinfty):
    
        #C = P0/(alpha-beta*P0)
        C = P0/(alpha-beta*P0)
        
    
        #Check if Pinfty is a solution
        maxVal = alpha/beta
        if(Pinfty > maxVal):
            print("The population density of the pest species will never reach Pinfty")
            return
    
        #else
        #initialize
        a = 0
        b = 1000
        tol = .001
        fa = logisticFunction(a, alpha, beta, C, Pinfty)
        fb = logisticFunction(b, alpha, beta, C, Pinfty)
        if(fa*fb >= 0):
            print("Hybrid Method Failed")
            return
        error = 10*tol
        #to give some values to x0 and x1 that are within a and b
        x0 = 1/3*(a+b)
        x1 = 1/2*(a+b)
    
        #loop
        while(error > tol):
            x2 = x0 - logisticFunction(x0, alpha, beta, C, Pinfty) * (x1-x0) /(logisticFunction(x1, alpha, beta, C, Pinfty) -
                 logisticFunction(x0, alpha, beta, C, Pinfty))
    
            #if x2 < a or >b
            if(abs(x2-x1) > .5*(b-a)):
                c=x1
                fc = logisticFunction(c, alpha, beta, C, Pinfty)
    
                #reduce the range by an order of magnitude
                for i in range(4):
                    if(fa*fc < 0):
                        b=c
                        fb=fc
                    else:
                        a=c
                        fa=fc
                    c = .5*(a+b)
                    fc = logisticFunction(c, alpha, beta, C, Pinfty)
                error = .5*(a+b)
                x0 = x1
                x1 = .5*(a+b)
                    
            #else secant method
            else:
                while(error > tol):
                    error = abs(x2 - x1)
                    x0 = x1
                    x1 = x2
                    x2 = x0 - logisticFunction(x0, alpha, beta, C, Pinfty) * (x1-x0)/(logisticFunction(x1, alpha, beta, C, Pinfty) -
                           logisticFunction(x0, alpha, beta, C, Pinfty))
    
        return x1
                    
                    
    
    
    def logisticFunction(t, alpha, beta, C, Pinfty):
        num = alpha*C*math.exp(alpha*t)
        den = 1+beta*C*math.exp(alpha*t)
        return num/den - Pinfty
    

   

**Last Modified:** November/2023




<hr>

<a id="TestRoutine"></a>

**Routine Name:**          Test Routine

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This routine tests all the other routines.

**Input:** This routine has no inputs.

**Output:** This routine will output the values of all the other routines given (as well as 2 that aren't. They are function iteration routines)

**Usage/Example:**

    main()





Output from the lines above:

    g(x) = x + f(x). Initial guess is  3.85
    3.85
    5.422500000000001
    13.71350625000001
    139.20622866878932
    18827.549185512107
    354401304.1361356
    1.2560028295578845e+17
    
    g(x) = x - f(x). Initial guess is  3.85
    3.85
    2.2774999999999985
    2.4779937499999973
    2.7275094749609354
    2.9257489137639356
    2.994486776192765
    2.9999696043632524
    2.9999999990761057
    2.9999999999999996
    
    g(x) = x + f(x). Initial guess is  1.1
    1.1
    2.81
    2.6561
    2.4304672099999993
    2.185302018885183
    2.034336838202924
    2.001179018457774
    2.0000013900845235
    2.0000000000019322
    1.9999999999999996
    
    g(x) = x - f(x). Initial guess is  1.1
    1.1
    -0.6099999999999999
    -10.0321
    -166.83563041000002
    -28841.141356762128
    -831984487.6088753
    -6.921981926137097e+17
    
    g(x) = x + f(x). Initial guess is  1.9
    1.9
    2.0099999999999993
    2.0000999999999993
    2.0000000099999995
    2.0000000000000004
    2.0000000000000004
    2.0000000000000004
    
    g(x) = x - f(x). Initial guess is  1.9
    1.9
    1.7900000000000005
    1.5359000000000012
    0.8564111900000042
    -1.5949729863571975
    -18.113776745352382
    -442.791568452583
    
    g(x) = x + f(x). Initial guess is  10.1
    10.1
    67.60999999999999
    4306.672099999997
    18530203.888518386
    343368382029252.3
    1.1790184577738519e+29
    1.3900845237714323e+58
    
    g(x) = x - f(x). Initial guess is  10.1
    10.1
    -47.40999999999999
    -2538.168099999999
    -6457532.312457605
    -41699762311633.945
    -1.7388701768470173e+27
    -3.0236694919279773e+54
    
    g(x) = x + epsilon * f(x). Initial guess is  3.1 . Epsilon is 0.1
    3.1
    3.111
    3.1233321000000003
    3.1371863906890414
    3.152787040336974
    3.170400132340165
    3.190343766084336
    
    g(x) = x - epsilon * f(x). Initial guess is  3.1 . Epsilon is 0.1
    3.1
    3.089
    3.0793079
    3.070748135699759
    3.0631727922592837
    3.056456432865172
    3.0504920566974687
    
    g(x) = x + epsilon * f(x). Initial guess is  2.9 . Epsilon is 5
    2.9
    2.4500000000000006
    1.2124999999999981
    8.250781250000019
    172.35820617675893
    144430.15923410846
    104296888188.12802
    
    g(x) = x - epsilon * f(x). Initial guess is  2.9 . Epsilon is 5
    2.9
    3.349999999999999
    0.9875000000000074
    -9.200781249999878
    -692.4921905517435
    -2415761.966830107
    -29179592211755.47
    
    g(x) = x + epsilon * f(x). Initial guess is  1.9 . Epsilon is 0.01
    1.9
    1.9011
    1.9021868121
    1.9032606181762717
    1.9043215970744654
    1.9053699246715845
    1.9064057739364353
    
    g(x) = x - epsilon * f(x). Initial guess is  1.9 . Epsilon is 0.01
    1.9
    1.8988999999999998
    1.8977867878999999
    1.8966601803717218
    1.895519990992231
    1.8943660301793306
    1.893198105125323
    
    g(x) = x + epsilon * f(x). Initial guess is  10.1 . Epsilon is 0.01
    10.1
    10.6751
    11.3409226001
    12.12004172430959
    13.042993752084165
    14.152040924649022
    15.507241501745979
    
    g(x) = x - epsilon * f(x). Initial guess is  10.1 . Epsilon is 0.01
    10.1
    9.524899999999999
    9.0339077999
    8.609488288524059
    8.238729817047933
    7.911899617917182
    7.621513043173061
    7.361714084659143
    7.127851446249412
    6.916181356163884
    6.723654778458589
    6.547762181582629
    6.386418394796123
    6.237875915402026
    6.100658751812601
    5.973511317342554
    5.8553585086254865
    5.745274201411432
    5.642456154987966
    5.546205848127748
    5.455912147436071
    5.37103798120247
    5.2911093903073985
    5.215706474020777
    5.144455857490394
    5.077024389668241
    5.013113842618791
    4.952456430759169
    4.89481100531145
    4.839959807799842
    4.787705688778655
    4.737869715593953
    4.690289106954228
    4.644815443233804
    4.601313110378461
    4.559657942499977
    4.519736034098945
    4.481442697624567
    4.44468154598487
    4.409363682831929
    4.375406986098754
    4.342735472463674
    4.311278732248915
    4.280971425789942
    4.251752833595139
    4.223566453695053
    4.196359640492024
    4.170083280193122
    4.144691498565316
    4.120141397310785
    4.096392815837984
    4.073408115613392
    4.051151984630612
    4.029591259836377
    4.008694765614699
    3.988433166656767
    3.9687788337407284
    3.9497057211162807
    3.931189254337908
    3.913206227520585
    3.895734709105555
    3.8787539553235355
    3.8622443306303325
    3.8461872344669876
    3.830565033764569
    3.8153610006738
    3.800559255052864
    3.7861447112938276
    3.7721030291099367
    3.75842056794323
    3.7450843446850044
    3.7320819944312076
    3.719401734021192
    3.7070323281318527
    3.6949630577202988
    3.6831836906271365
    3.6716844541694758
    3.6604560095680516
    3.6494894280666257
    3.6387761686142563
    3.6283080569922186
    3.618077266277483
    3.6080762985438177
    3.598297967709872
    
    Testing Bijection Method 2.00703125
    Testing Bijection Method 2.996966552734375
    Testing Newton's Method 3.0000152590218967
    Testing Newton's Method 1.9999847409781035
    Testing Secant Method 1.9949042720195709
    Testing Secant Method 3.0079437482011504
    Testing a hybrid Method 3.0001240687574207
    Testing a hybrid Method 1.9999929119280258
    Testing alternate hybrid Method 1.999990234293511
    30.32588873949352
    The population density of the pest species will never reach Pinfty
    None
    41.753458007592016
    The population density of the pest species will never reach Pinfty
    None
    The population density of the pest species will never reach Pinfty
    None

  


     

**Implementation/Code:** The following is the code for Test Routine

        #import functions
        from BijectiveMethod import *
        from FixedPointIterationProb1 import *
        from FixedPointIterationProb2 import *
        from NewtonsMethod import *
        from SecantMethod import *
        from HybridMethod import *
        from HybridMethodAlt import *
        from LogisticODERoots import *
        
        #test functions
        def main():
            #test FixedPointIterationProb1
            fixedPointIteration(function1, 3.85, .0000001)
            fixedPointIteration(function1, 1.1 , .0000001)
            fixedPointIteration(function1, 1.9, .0000001)
            fixedPointIteration(function1, 10.1, .0000001) 
        
            #test FixedPointIterationProb2
            fixedPointIterationEpsilon(function1, 3.1, .01, .1)
            fixedPointIterationEpsilon(function1, 2.9, .01, 5)
            fixedPointIterationEpsilon(function1, 1.9, .01, .01)
            fixedPointIterationEpsilon(function1, 10.1, .01, .01)
        
            #test BijectionMethod
            val = bisection(function1, 1.5, 2.6, .01)
            print("Testing Bijection Method", val)
            val = bisection(function1, 2.5, 100.6, .01)
            print("Testing Bijection Method", val)
        
            #test NewtonsMethod
            val = newtonsMethod(function1, derivativeFunction1, 4, .01, 30)
            print("Testing Newton's Method", val)
            val = newtonsMethod(function1, derivativeFunction1, 1, .01, 30)
            print("Testing Newton's Method", val)
        
            #test SecantMethod
            val = secantMethod(function1, .01, 0, .01, 30)
            print("Testing Secant Method", val)
            val = secantMethod(function1, 7.01, 80, .01, 30)
            print("Testing Secant Method", val)
        
            #test HybridMethod
            val = hybridMethod(function1, 2.4, 1111, .01, 30)
            print("Testing a hybrid Method", val)
            val = hybridMethod(function1, -72.4, 2.84, .01, 30)
            print("Testing a hybrid Method", val)
        
            #test HybridMethodAlt
            val = hybridMethodAlt(function1, -200, 3001, .01, 30)
            print("Testing alternate hybrid Method", val)
        
            #test logisticODERoots
            val1 = logisticODERoots(.1, .001,2,29.75)
            print(val1)
            val2 = logisticODERoots(.1, .001, 2, 115.35)
            print(val2)
            val3 = logisticODERoots(.1, .0001, 2, 115.346)
            print(val3)
            val4 = logisticODERoots(.1, .001, 2, 155.346)
            print(val4)
            val5 = logisticODERoots(.1, .01, 100, 155.346)
            print(val5)
            
        main()
        

   

**Last Modified:** November/2023
<hr>

<a id="upperTriangular"></a>

**Routine Name:**       upperTriangular    

**Author:** Carter Green

**Language:** python. This code can be imported using import statements   

**Description/Purpose:** This code takes a matrix "a" and a vector "b" of the same size and performs Gaussian elimination to reduce "a" to an upper triangular matrix.

**Input:** This code has two inputs. a - the nxn matrix. b - the n vector

**Output:** This code has no outputs. It directly manipulates a and b.

**Usage/Example:**
This function is commonly used in addition to a backward substitution routine in order to solve a system of equations. Here is some code that creates a matrix A and a vector y and solves for the vector x2 that is the solution.

    n = 10
    A = makeMatrix(n)
    x = makeVectorOne(n)
    y = matrixTimesVector(A, x)
    upperTriangular(A, y)
    x2 = back_sub(A, y)



Output from the lines above:

There will be no output from the lines above. The final vector x2 will be the solution to Ax = b.
  


     

**Implementation/Code:** The following is the code for upperTriangular

    #reduce a square matrix to an upper triangular matrix
    def upperTriangular(a, b):
        #find length of matrix
        n = len(b)

        #loop through all pivots
        for k in range(n - 1):

            #loop through all rows beneath pivot
            for i in range(k + 1, n):

                #calculate factor and then loop through all elements in the row after
                    #to help reduce matrix
                factor = -a[i][k] / a[k][k]
                for j in range(k + 1, n):
                    a[i][j] += factor * a[k][j]

                #fix the solution vector
                b[i] += factor * b[k]
     

**Last Modified:** October/2023




<hr>

<a id="backSubstitution"></a>

**Routine Name:**      backSubstitution     

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This code is used to solve an upper triangular matrix for the vector x such that Ax = b.

**Input:** This code has two inputs. a - the nxn matrix. b - the n vector

**Output:** This code has out output - the vector x that solves ax = b

**Usage/Example:**
This code is commonly used in combination with some sort of Gaussian Elimination that will make a matrix upper triangular. Here is an example:

    n = 10
    A = makeMatrix(n)
    x = makeVectorOne(n)
    y = matrixTimesVector(A, x)
    upperTriangular(A, y)
    x2 = back_sub(A, y)




Output from the lines above:

There will be no output from this code - the final vector x2 will be the solution to Ax = b

  


     

**Implementation/Code:** The following is the code for backSubstitution

    def backSubstitution(a, b):
    n = len(b)
    x = [0] * n  # Initialize solution vector x

    # Solve for the last term manually
    x[n - 1] = b[n - 1] / a[n - 1][n - 1]

    # Loop back up through the array
    for i in range(n - 2, -1, -1):
        # Calculate each x[i]
        total = 0.0
        for j in range(i + 1, n):
            total += a[i][j] * x[j]
        x[i] = (b[i] - total) / a[i][i]

    return x
     

**Last Modified:** October/2023



<hr>

<a id="LUfactorization"></a>

**Routine Name:**     LUfactoriaztion     

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This program will take a matrix A and reduce it into it's LU factorization. The L will be in the lower triangular section of A with an implied diagonal of 1's and the U will be the upper triangular section of A.

**Input:** There is one input for this function. A is the nxn matrix that is being reduced to it's LU factorization

**Output:** There is no output for this function. Instead the matrix A is being directly manipulated.

**Usage/Example:**
This function is used to solve systems of equations. It is especially good when the same matrix will be solved for multiple vectors. Here is an example where a matrix is created and a b-vector is created and then solved for using an LU factorization. It is commonly used with a backsubstitution and forward substitution routine to solve for Ax=b.

    n = 10
    A = makeMatrix(n)
    x = makeVectorOne(n)
    y = matrixTimesVector(A,x)
    LUfactorization(A)
    y2 = forwardSubDiagEquals1(A, y)
    x2 = backSubstitution(A, y2)
    print(twoNormDistance(x, x2))




Output from the lines above:
This will print out the error between the the actual vector of ones and the x2 vector that is gotten from the LU factorization and substitution. Here is an example of an output:
  
    7.195067539997724e-16

     

**Implementation/Code:** The following is the code for LUfactorization

    #reduce matrix to LU factorization
    def LUfactorization(A):
        #find length of matrix
        n = len(A)

        #loop through all pivots
        for k in range(n - 1):

            #loop through all rows beneath pivot
            for i in range(k + 1, n):

                #calculate U section of LU matrix
                #calculate factor and then loop through all elements in the row after
                    #to help reduce matrix
                factor = -A[i][k] / A[k][k]
                for j in range(k + 1, n):
                    A[i][j] += factor * A[k][j]

                #calculate L section of LU matrix
                #Note: for L section of LU matrix, a[i][j] = -factor from the
                    #jth pivot and ith row
                A[i][k] = -factor

   

**Last Modified:** October/2023




<hr>

<a id="forwardSubstitution"></a>

**Routine Name:**          forwardSubstitution

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This code is used to solve a lower triangular matrix. So given a lower triangular matrix A and a vector b, this function will solve for x in Ax=b.

**Input:** This function has 2 inputs. It requires a lower triangular matrix A and a vector b.

**Output:** This function will return the vector x that is the solution to Ax=b.

**Usage/Example:**
This code is often used in combination with an LU factorization routine and a back substitution routine to solve a matrix Ax = b. This is done by factoring A into an LU factorization, getting LUx=b, then letting y=Ux then solving Ly=b and Ux=y. Note that this forward substitution routine would give the wrong output values because the LUfactorization routine assumes that the diagonal is ones, and puts the U part of the LU factorization in the diagonal. So, the forwardSubDiagEquals1 substituion (which is very similar) is used below:

    n = 10
    A = makeMatrix(n)
    x = makeVectorOne(n)
    y = matrixTimesVector(A,x)
    LUfactorization(A)
    y2 = forwardSubDiagEquals1(A, y)
    x2 = backSubstitution(A, y2)
    print(twoNormDistance(x, x2))





Output from the lines above:
This will print out the error between the the actual vector of ones and the x2 vector that is gotten from the LU factorization and substitution. Here is an example of an output:

      7.195067539997724e-16


     

**Implementation/Code:** The following is the code for forwardSubstitution

    #forward substitution
    def forwardSubstitution(A, b):
    
        n = len(A)
        #initialize solution vector x
        x = [0] * n

        #solve for the first term manually
        x[0] = b[0] / A[0][0]

        #Loop down through the array
        for i in range(1, n):
            #calculate each x[i]
            total = 0
            for j in range(0, i):
                total += A[i][j] * x[j]
            x[i] = (b[i] - total) / A[i][i]

        return x
    
   

**Last Modified:** October/2023







<hr>

<a id="forwardSubDiagEquals1"></a>

**Routine Name:**    forwardSubDiagEquals1      

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This code is used to solve a lower triangular matrix when the diagonal is assumed to be 1's, but may in fact be something else. So given a lower triangular matrix A and a vector b, this function will solve for x in Ax=b.

**Input:** This function has 2 inputs. It requires a lower triangular matrix A and a vector b.

**Output:** This function will return the vector x that is the solution to Ax=b.

**Usage/Example:**
This code is often used in combination with an LU factorization routine and a back substitution routine to solve a matrix Ax = b. This is done by factoring A into an LU factorization, getting LUx=b, then letting y=Ux then solving Ly=b and Ux=y.

    n = 10
    A = makeMatrix(n)
    x = makeVectorOne(n)
    y = matrixTimesVector(A,x)
    LUfactorization(A)
    y2 = forwardSubDiagEquals1(A, y)
    x2 = backSubstitution(A, y2)
    print(twoNormDistance(x, x2))



Output from the lines above:
This will print out the error between the the actual vector of ones and the x2 vector that is gotten from the LU factorization and substitution. Here is an example of an output:

      7.195067539997724e-16

  


     

**Implementation/Code:** The following is the code for 

    #forward substitution
    def forwardSubDiagEquals1(A, b):
        #assume the diagonal is 1
        n = len(A)
        #initialize solution vector x
        x = [0] * n

        #solve for the first term manually
        x[0] = b[0] 

        #Loop down through the array
        for i in range(1, n):
            #calculate each x[i]
            total = 0
            for j in range(0, i):
                total += A[i][j] * x[j]
            x[i] = (b[i] - total)

        return x
    

   

**Last Modified:** October/2023






<hr>

<a id="matrixTimesVector"></a>

**Routine Name:**          matrixTimesVector

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This function will take a matrix A and a vector x and perform the matrix multiplication Ax = b.

**Input:** This function has two inputs. The first input is the nxn matrix A. The second input is the vector x of length n.

**Output:** This function will return the vector b (of length n) that is the result of matrix multiplication Ax

**Usage/Example:**
This snippet of code uses other function to make a diagonally dominant matrix of size n and a vector of ones. Then, this function is called to matrix multiply the matrix by the vector.

    n = 10
    A = makeMatrix(n)
    x = makeVectorOne(n)
    y = matrixTimesVector(A, x)




Output from the lines above:
There is no output from the code above, but the variable y will be the vector that is the result of Ax=y.
  


     

**Implementation/Code:** The following is the code for matrixTimesVector

    #multiply a matrix by a vector
    def matrixTimesVector(A,x):
        #check they're the same length
        if(len(A) != len(x)):
            print("Vector and matrix are different sizes")
            return
        y = []
        for i in range(len(A)):
            sum = 0

            #sum up each dot product
            for j in range(len(A[0])):
                sum += A[i][j] * x[j]
            y.append(sum)

        #return vector
        return y
        

   

**Last Modified:** October/2023






<hr>

<a id="makeMatrix"></a>

**Routine Name:**      MakeMatrix    

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This function makes a diagonally dominant matrix of size n. It does this by randomly assigning numbers between 0 and 1 to all the entries in the matrix except for the diagonal, which is assigned the value n.

**Input:** There is one input to this function, n, which is the size of the diagonally dominant matrix that will be created.

**Output:** There is one output to this function, y, which is the diagonally dominant nxn matrix.

**Usage/Example:**
This snippet of code creates a diagonally dominant matrix.

    n = 10
    A = makeMatrix(n)
    




Output from the lines above:
There will be no output from the lines above. However, the variable A will be the diagonally dominant matrix of size n.

  


     

**Implementation/Code:** The following is the code for makeMatrix

    import random
    #make a diagonally dominant matrix of size n
    def makeMatrix(n):
        random.seed()
        y = []
        #make the rows
        for i in range(n):
            x = []
            #make the columns
            for j in range(n):

                #if on diagonal
                if(i == j):
                x.append(n)

                #else
                else:
                    x.append(random.random())
            y.append(x)

        #return the array
        return y

   

**Last Modified:** October/2023





<hr>

<a id="makeVectorOne"></a>

**Routine Name:**          makeVectorOne

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This function makes a vector of length n that is made entirely of ones.

**Input:** There is one input to this function, n, which is the length of the vector.

**Output:** There is one output to this funcion, x, which is the vector of length n.

**Usage/Example:**
This snippet of code will create a vector of length n that has all entries equal to 1.

    n = 10
    x = makeVectorOne(n)



Output from the lines above:
There is no output from the lines above. However, x will be the vector of length n that has all entries equal to 1.

  


     

**Implementation/Code:** The following is the code for makeVectorOne

   
    #make a vector of size n of ones
    def makeVectorOne(n):
        x = []
        for i in range(n):
            x.append(1)
        return x


**Last Modified:** October/2023





<hr>

<a id="TestingOfMatrixOperations"></a>

**Routine Name:**      TestingOfMatrixOperations    

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** This function tests all the other functions involved with solving a system of equations. This includes: upperTriangular, backSubstitution, LUfactorization, and forwardSubDiagEquals1. Note that for the import statements to work, this .py file must be in the same folder as all other .py files that include the functions called.

**Input:** This routine has no inputs.

**Output:** This routine will output 30 errors for the computing gaussian elimination and solving for a vector of ones and 30 errors for computing LU factorization and solving for a vector of ones.

**Usage/Example:**
The only example is the code itself. It is designed for testing purposes only. To see an example, look at what is in the code.




Output from the lines above:
The following are example outputs for what this function will output:

    7.364386412590295e-16
    7.447602459741819e-16
    7.364386412590295e-16
    6.753223014464259e-16
    4.440892098500626e-16
    6.753223014464259e-16
    5.087681048627601e-16
    8.455206652451151e-16
    7.447602459741819e-16
    6.280369834735101e-16
    6.473657049138938e-16
    5.874748045952207e-16
    7.447602459741819e-16
    9.019494489765868e-16
    7.771561172376096e-16
    7.447602459741819e-16
    6.377745716588144e-16
    3.1401849173675503e-16
    7.771561172376096e-16
    5.874748045952207e-16
    1.0295784775289034e-15
    8.95090418262362e-16
    8.15843973306311e-16
    7.021666937153402e-16
    9.485749680535094e-16
    8.3820000221454525e-16
    9.019494489765868e-16
    5.874748045952207e-16
    7.850462293418876e-16
    6.753223014464259e-16
    
    
    5.438959822042073e-16
    8.08254562088053e-16
    7.108895957933346e-16
    6.280369834735101e-16
    6.753223014464259e-16
    6.661338147750939e-16
    5.087681048627601e-16
    7.108895957933346e-16
    6.280369834735101e-16
    5.438959822042073e-16
    9.805224261780596e-16
    6.377745716588144e-16
    7.447602459741819e-16
    4.965068306494546e-16
    6.377745716588144e-16
    7.364386412590295e-16
    6.753223014464259e-16
    5.438959822042073e-16
    1.0053497077208614e-15
    4.440892098500626e-16
    5.438959822042073e-16
    8.08254562088053e-16
    5.438959822042073e-16
    7.850462293418876e-16
    4.965068306494546e-16
    7.447602459741819e-16
    7.691850745534255e-16
    7.771561172376096e-16
    7.364386412590295e-16
    7.021666937153402e-16
  


     

**Implementation/Code:** The following is the code for TestingOfMatrixOperations

    
    from makeMatrix import *
    
    from matrixTimesVector import *
    
    from upperTriangular import *
    
    from twoNormDistance import *
    
    from backSubstitution import *
    
    from forwardSubstitution import *
    
    from LUfactorization import *
    
    from forwardSubDiagEquals1 import *
    
    from backSubstitution import *
    
    
    
    def testUpperTriangular():
        for i in range(30):
            n = 10
            A = makeMatrix(n)
            x = makeVectorOne(n)
            y = matrixTimesVector(A, x)
            upperTriangular(A, y)
            x2 = backSubstitution(A, y)
            print(twoNormDistance(x, x2))
    
    
    def testLUfactorization():
        for i in range(30):
            n = 10
            A = makeMatrix(n)
            x = makeVectorOne(n)
            y = matrixTimesVector(A,x)
            LUfactorization(A)
            y2 = forwardSubDiagEquals1(A, y)
            x2 = backSubstitution(A, y2)
            print(twoNormDistance(x, x2))
    
    
    def main():
        testUpperTriangular()
        print()
        print()
        testLUfactorization()
        
    main()

   

**Last Modified:** October/2023



<hr>
<a id="smaceps"></a>

**Routine Name:**    smaceps  

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc smaceps.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc smaceps.c -o smaceps.exe

**Description/Purpose:** This routine will compute the single precision value for the machine epsilon or the number of digits
in the representation of real numbers in single precision. This is a routine for analyzing the behavior of any computer. This
usually will need to be run one time for each computer.

**Input:** There are no inputs needed in this case. 

**Output:** This routine returns a single precision value for the number of decimal digits that can be represented on the
computer being queried.

**Usage/Example:**

This routine has no parameters than need to be called. It will return single precision machine epsilon as a float. Here is an example of it being used in code:

      float a = smaceps();
      printf("smaceps = %10.10e\n", a);

Output from the lines above:

      smaceps = 5.9604644775e-08

The value is the number of digits that can be represented in a float number. The number of decimal digits that can be represented is roughly eight (E-08 on the
end of the second value).

**Implementation/Code:** The following is the code for smaceps()

     float smaceps(){
	//initial conditions
	float one = 1;
	float a = 1;
	float appone = one + a;
	
	//get a small enough that the computer can't tell 1+a apart from 1
	while(one != appone){
		a *= (float) 1/2;
		appone = (float) one + (float) a;
	}
	return a;
    }


**Last Modified:** October/2023




<hr>
<a id="dmaceps"></a>

**Routine Name:**   dmaceps      

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc dmaceps.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc dmaceps.c -o dmaceps.exe

**Description/Purpose:** This routine will calculate the double precision value for machine epsilon. This will usually need to be run one time on each computer.

**Input:** This program takes no inputs. 

**Output:** This program returns double precision machine epsilon for the computer the code is ran on. It returns it as a double.

**Usage/Example:**
    
    double b = dmaceps();
	printf("dmaceps = %10.10e\n", b);


Output from the lines above:

    dmacpes = 1.11002230246e-16     

**Implementation/Code:** The following is the code for dmaceps()

     double dmaceps(){
	//initial conditions
	double one = 1;
	double a = 1;
	double appone = one + a;
	
	//get a small enough that the computer can't tell 1+a apart from 1
	while(one != appone){
		a *= (double) 1/2;
		appone = (double) 1 + (double) a;
	}
	return a;
}

**Last Modified:** October/2023



<hr>
<a id="l2-norm"></a>

**Routine Name:**  l2-norm         

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc l2-norm.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc l2-norm.c -o l2norm.exe

**Description/Purpose:** This program calculates the length of a vector in the l2-norm (Euclidean norm).

**Input:** This program takes 2 inputs. The first is the vector as an array. The second is the size of the array

**Output:** This program outputs the l2-norm of a vector as a double.

**Usage/Example:**

    double vector[] = {1.2,5,7};
	double distance2 = l2_norm(vector, sizeof(vector)/sizeof(vector[0]));
	printf("l-2 norm is %f. l-2 norm should be 8.685\n", distance2);

Output from the lines above:

     l-2 norm is 8.685620. L-2 norm should be 8.685

**Implementation/Code:** The following is the code for l2norm(double vector[], int n)

    double l2_norm(double vector[], int n){
	double listSquared = 0;
	
	//add each term squared to a sum
	for(int i = 0; i < n ; i++){
		listSquared += vector[i] * vector[i];
	}
	
	//return square root of the sum
	return sqrt(listSquared);
    }

     

**Last Modified:** October/2023

<hr>
<a id="l1-norm"></a>
**Routine Name:**   l1-norm        

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc l1-norm.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc l1-norm -o l1-norm.exe

**Description/Purpose:** This vector calculates and returns the length of a vector in l1.

**Input:** This program has 2 inputs. The first one is the vector as an array. The second is the length of that array

**Output:** This program outputs the length of the vector in the l1-norm.

**Usage/Example:**

    double vector[] = {1.2,5,7};
    double distance1 = l1_norm(vector, sizeof(vector)/sizeof(vector[0]));
	printf("l1-norm is %f. l1-norm should be 13.2\n", distance1);

Output from the lines above:

    L1-norm is 13.200000. L1-norm should be 13.2.
     

**Implementation/Code:** The following is the code for l1_norm(double vector[], int n)

    double l1_norm(double vector[], int n){
	double listSummed = 0;
	
	//add each term's absolut to a sum
	for(int i = 0; i < n ; i++){
		listSummed += fabs(vector[i]);
	}
	
	//return square root of the sum
	return listSummed;
}
     

**Last Modified:** October/2023


<hr>

<a id="linf-norm"></a>

**Routine Name:**   linf-norm        

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc linf-norm.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc linf-norm.c -o linf-norm.exe

**Description/Purpose:** This program calculates and returns the length of vector in the l-inf norm. 

**Input:** This program has 2 inputs. The first one is the vector as an array. The second is the length of the array/vector.

**Output:** This program outputs the length of the vector in the l-inf norm.

**Usage/Example:**

    double vector[] = {1.2,5,7};
    double distanceinf = linf_norm(vector, sizeof(vector)/sizeof(vector[0]));
	printf("linf-norm is %f. linf-norm should be 7\n", distanceinf);

Output from the lines above:

     linf-norm is 7.000000. linf-norm should be 7

**Implementation/Code:** The following is the code for linf_norm(double vector[], int n)

    double linf_norm(double vector[], int n){
	double largestVal = 0;
	
	//find largest number
	for(int i = 0; i < n ; i++){
		if(largestVal < fabs(vector[i])){
			largestVal = fabs(vector[i]);
		}
	}
	
	//return square root of the sum
	return largestVal;
    }

     

**Last Modified:** October/2023


<hr>

<a id="l2-distance"></a>

**Routine Name:**    l2-distance       

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc l2-distance.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc l2-distance.c -o l2-distance.exe

**Description/Purpose:** This program calculates and returns the distance between 2 vectors in the l-2 norm.

**Input:** This program has 3 inputs. The first one is the first vector, as an array. The second input is the second vector, as an array. The third input is the length of the arrays/vectors

**Output:** This program outputs the distance between the two vectors as defined by the l2-norm

**Usage/Example:**

    double vector[] = {1.2,5,7};
    double vector2[] = {2.5,4,7.0004};
	double distance2_2 = l2_distance(vector, vector2, sizeof(vector)/sizeof(vector[0]));
	printf("l2-distance is %f. l2-distance should be 1.640\n", distance2_2);

Output from the lines above:

     l2-distance is 1.640122. l2-distance should be 1.640

**Implementation/Code:** The following is the code for l2_distance(double vector1[], double vector2[], int n)

    double l2_distance(double vector1[], double vector2[], int n){
	double listSquared = 0;
	
	//add each term squared to a sum
	for(int i = 0; i < n ; i++){
		listSquared += (vector2[i] - vector1[i]) * (vector2[i] - vector1[i]);
	}
	
	//return square root of the sum
	return sqrt(listSquared);
    }
     

**Last Modified:** October/2023



<hr>

<a id="l1-distance"></a>

**Routine Name:**     l1-distance      

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc l1-distance.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc l1-distance.c -o l1-distance.exe

**Description/Purpose:** This program calculates and returns the distance between two vectors in the l1-norm.

**Input:** This program has 3 inputs. The first two are the vectors that are going to be compared. They are inputted as arrays. The third one is the length of the arrays/vectors.

**Output:** This program outputs the distance between the two vectors, as defined in the l1-norm.

**Usage/Example:**

    double vector[] = {1.2,5,7};
    double vector2[] = {2.5,4,7.0004};
    double distance1_2 = l1_distance(vector, vector2, sizeof(vector)/sizeof(vector[0]));
	printf("l1-distance is %f. l1-distance should be 2.3004\n", distance1_2);


Output from the lines above:

    l1-distance is 2.300400. l1-distance should be 2.3004

     

**Implementation/Code:** The following is the code for l1_distance(double vector1[], double vector2[], int n)

    double l1_distance(double vector1[], double vector2[], int n){
	double listSummed = 0;
	
	//add each term's absolut to a sum
	for(int i = 0; i < n ; i++){
		listSummed += fabs(vector2[i] - vector1[i]);
	}
	
	//return square root of the sum
	return listSummed;
    }
     

**Last Modified:** October/2023



<hr>

<a id="linf-distance"></a>

**Routine Name:**      linf-distance     

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc linf-distance.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc linf-disstance.c -o linf-distance.exe

**Description/Purpose:** This program calculates and returns the distance between two vectors as defined by the linf-norm.

**Input:** This program has 3 inputs. The first two are the vectors that are going to be compared. The third input is the length of the vectors/arrays.

**Output:** This program will output the distance between the two vectors, as defined by the linf-norm.

**Usage/Example:**

    double vector[] = {1.2,5,7};
    double vector2[] = {2.5,4,7.0004};
    double distanceinf_2 = linf_distance(vector, vector2, sizeof(vector)/sizeof(vector[0]));
	printf("linf-disstance is %f. l-inf distance should be 1.3\n", distanceinf_2);

Output from the lines above:

     linf-distance is 1.300000. l-inf distance should be 1.3

**Implementation/Code:** The following is the code for linf_distance(double vector1[], double vector2[], int n)

    double linf_distance(double vector1[], double vector2[], int n){
	double largestVal = 0;
	
	//find largest number
	for(int i = 0; i < n ; i++){
		if(largestVal < fabs(vector1[i] - vector2[i])){
			largestVal = fabs(vector1[i] - vector2[i]);
		}
	}
	
	//return square root of the sum
	return largestVal;
    }    

     

**Last Modified:** October/2023



<hr>

<a id="linreg"></a>

**Routine Name:**     linreg      

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc linreg.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc linreg.c -o linreg.exe

**Description/Purpose:** This program calculates and returns the linear regression equation for a pair of vectors of the same list (can be thought of as a list of ordered pairs).

**Input:** This program takes in 4 inputs. The first one, (*vectorA) is the independent variable vector. The second one, (*vectorY) is the dependent variable vector. The third one (*vectorX) is the output vector. The fourth input is the size (length) of the variable vectors.

**Output:** The program has no outputs. Rather, it stores the values that would be output into an output vector

**Usage/Example:**

    double vectorA[] = {1, 2, 3, 4};
    double vectorY[] = {1, 2, 98, 6};
    int n = sizeof(vectorA) / sizeof(vectorA[0]);
    double vectorX[2];

    linreg(vectorA, vectorY, vectorX, n);

    printf("Resulting vectorX: [%lf, %lf]\n", vectorX[0], vectorX[1]);



Output from the lines above:

    Resulting vectorX: [1.600000, 0.600000]
     

**Implementation/Code:** The following is the code for linreg(const double *vectorA, const double *vectorY, double *vectorX, int n)

    void linreg(const double *vectorA, const double *vectorY, double *vectorX, int n) {
    double y1 = 0, y2 = 0;
    double a = n, b = 0, c = 0, d = 0;
    double detA, x1, x2;
    
    // Calculate y1 and y2
    for (int i = 0; i < n; i++) {
        y1 += vectorY[i];
        y2 += vectorA[i] * vectorY[i];
    }

    // Calculate a, b, c, and d
    for (int i = 0; i < n; i++) {
        b += vectorA[i];
        d += vectorA[i] * vectorA[i];
    }
    c = b;

    // Calculate the determinant of (A^T * A)
    detA = 1 / (a * d - b * c);

    // Calculate x vector
    x1 = detA * (d * y1 - b * y2);
    x2 = detA * (a * y2 - c * y1);

    // Store the results in the output vector
    vectorX[0] = x1;
    vectorX[1] = x2;
    }
     

**Last Modified:** October/2023



<hr>

<a id="fordiff"></a>

**Routine Name:**    fordiff       

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc fordiff.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc fordiff.c -o fordiff.exe

**Description/Purpose:** This program approximates the derivative of a function at a point a using the forward difference.

**Input:** This program has 3 inputs. The first one is the pointer to the function that's derivative will be approximated. The second input is the point at which the derivative will be approximated. The third input is the step value of h that will be used in the approximation.

**Output:** This program will output the approximation of the function at point a.

**Usage/Example:**

    //Example function to be used with fordiff
    double exampleFunc(double x) {
    return x * x; // You should replace this with your actual function
    }

    int main() {
    double a = 2.0; // Value of 'a'
    double h = 0.001; // Value of 'h'

    // Using the forward difference function with the example function
    double result = fordiff(exampleFunc, a, h);

    printf("Result: %lf\n", result);
    }


Output from the lines above:

    Result: 4.001000
     

**Implementation/Code:** The following is the code for backdiff(double (*func)(double), double a, double h)

    double backdiff(double (*func)(double), double a, double h) {
    return (func(a) - func(a - h)) / h;
    }
     

**Last Modified:** October/2023



<hr>

<a id="backdiff"></a>

**Routine Name:**           backdiff

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc backdiff.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc backdiff.c -o backdiff.exe

**Description/Purpose:** This program approximates the derivative of a function at a point a using the backward difference.

**Input:** This program has 3 inputs. The first one is the pointer to the function that's derivative will be approximated. The second input is the point at which the derivative will be approximated. The third input is the step value of h that will be used in the approximation.

**Output:** This program will output the approximation of the function at point a.

**Usage/Example:**

    // Example function to be used with backdiff
    double exampleFunc(double x) {
    return x * x; // Replace this with your actual function
    }

    int main() {
    double a = 2.0; // Value of 'a'
    double h = 0.001; // Value of 'h'

    // Using the backward difference function with the example function
    double result = backdiff(exampleFunc, a, h);

    printf("Result: %lf\n", result);

    return 0;
    }



Output from the lines above:

    Result: 3.999000


     

**Implementation/Code:** The following is the code for backdiff(double (*func)(double), double a, double h)

    double backdiff(double (*func)(double), double a, double h) {
    return (func(a) - func(a - h)) / h;
    }
     

**Last Modified:** October/2023




<hr>

<a id="upperTriangular1"></a>

**Routine Name:**           upperTriangular

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc upperTriangular.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc upperTriangular.c -o upperTriangular.exe

**Description/Purpose:** This code will reduce a matrix to upper triangular form.

**Input:** This program has 3 inputs. The first one is a NxN matrix. N must be given/defined at compilation or dynamic memory allocation should be used. The second input is the solution vector. The third input is the length of the array and solution vector.

**Output:** This program has no outputs. Instead, it directly modifies the vectors.

**Usage/Example:**
In this case, N was defined to be 3

    #define N 3
    
    double a[][N] = {{2, 1, 1}, {1, 3, 2}, {1, 0, 0}};
    double b[] = {4, 5, 6};
    int n = 3;

    upper_triangular(a, b, n);

    // Print the upper triangular matrix a
    printf("Upper Triangular Matrix a:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }

    // Print the updated vector b
    printf("Updated Vector b:\n");
    for (int i = 0; i < n; i++) {
        printf("%lf\n", b[i]);
    }

    return 0;
    }



Output from the lines above:

    Upper Triangular Matrix a:
    2.000000 1.000000 1.000000
    0.000000 2.500000 1.500000
    0.000000 0.000000 -0.400000

    Updated Vector b:
    4.000000
    3.000000
    -0.500000


     

**Implementation/Code:** The following is the code for upper_triangular(double[][N] a, double[] b, int n)

    void upper_triangular(double[][N] a, double[] b, int n){
	
	//loop through all pivots
	for(int k = 0; k < (n-1); k++){
		
		//loop through all rows beneath the pivots
		for(int i = k+1; i < n; i++){
			
			//calculate factor
			double factor = (-a[i][k])/a[k][k];
			
			//loop through all elements in row (that are right of pivots
			for(int j = k+1; j < n; j++){
				
				//calculate new term
				a[i][j] = a[i][j] + factor * a[k][j];
			}
			
			//calculate new term in solution vector
			b[i] = b[i] + factor * b[k];
		}
	}
    }    

     

**Last Modified:** October/2023




<hr>

<a id="backSubstitution1"></a>

**Routine Name:**           backSubstitution

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc backSubstitution.c -o

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc backSubstitution.c -o backSubstitution.exe

**Description/Purpose:** This code takes in an upper triangular matrix and solves it for the x values

**Input:** This code has 4 inputs. The first one is the upper triangular matrix a. Then there's the output vector b. The third input is the x-matrix that is being solved for. The fourth input is the length of the upper triangular matrix and output vector

**Output:** This code has no outputs because it directly changes the vectors themselves

**Usage/Example:**

    #define N 3
    int main() {
    double a[N][N] = {{2.0, 3.0, 5.0},
                     {0.0, 4.0, 7.0},
                     {0.0, 0.0, 6.0}};
    double b[N] = {10.0, 14.0, 18.0};
    double x[N];  // The solution vector

    back_sub(a, b, x, N);

    printf("Solution Vector x:\n");
    for (int i = 0; i < N; i++) {
        printf("x[%d] = %lf\n", i, x[i]);
    }

    return 0;
    }


Output from the lines above:

    Solution Vector x:
    x[0] = 1.000000
    x[1] = 2.000000
    x[2] = 3.000000



**Implementation/Code:** The following is the code for back_sub(double[][N] a, double[] b, double[] x, int n)

    void back_sub(double[][N] a, double[] b, double[] x, int n){
	//solve first term manually
	x[n-1] = b[n-1]/a[n-1][n-1];
	
	//loop back up through array
	for(int i = n-2; i >=0; i--){
		
		//calculate each x[n]
		double sum = 0.0;
		for(int j = i+1; j<n; j++){
			
				sum = sum + a[i][j] *x[j];
		}
		x[i] = (b[i] - sum)/a[i][i];
	}
    }
     

**Last Modified:** October/2023


<hr>

<a id="jacobi"></a>

**Routine Name:**     jacobi     

**Author:** Carter Green

**Language:** C. This program can be compiled with the gcc compiler. Note that it can be compiled with the -fopenmp addition. 

**Description/Purpose:** This program solves the linear system Ax = b for x.

**Input:** This method has no input.

**Output:** This method will output a lot of x-vectors to the screen as well as residuals. Then it will output the number of times it solved the system as well as the total time taken to do so.

**Usage/Example:**

    ./jacobi.exe

Output from the lines above:

There will be 1000 copies of an output vector. Then there will be the total time it took to print them.
An example of the output vector is 

    x[0] = 1.000270          res[0] = -0.032290
    x[1] = 1.000268          res[1] = -0.033692
    x[2] = 1.000262          res[2] = -0.031283
    x[3] = 1.000229          res[3] = -0.024364
    x[4] = 1.000265          res[4] = -0.028551
    x[5] = 1.000255          res[5] = -0.031552
    x[6] = 1.000268          res[6] = -0.033071
    x[7] = 1.000253          res[7] = -0.033491
    x[8] = 1.000256          res[8] = -0.030943
    x[9] = 1.000264          res[9] = -0.032340
    x[10] = 1.000260         res[10] = -0.039420
    x[11] = 1.000246         res[11] = -0.028835
    x[12] = 1.000273         res[12] = -0.032660
    x[13] = 1.000275         res[13] = -0.033664
    x[14] = 1.000250         res[14] = -0.037572
    x[15] = 1.000250         res[15] = -0.023670
    x[16] = 1.000237         res[16] = -0.024610
    x[17] = 1.000247         res[17] = -0.029123
    x[18] = 1.000252         res[18] = -0.031054
    x[19] = 1.000268         res[19] = -0.036663
    x[20] = 1.000246         res[20] = -0.036405
    x[21] = 1.000244         res[21] = -0.033246
    x[22] = 1.000263         res[22] = -0.029514
    x[23] = 1.000260         res[23] = -0.034821
    x[24] = 1.000245         res[24] = -0.032074

  An example of the total time taken is 

    Time taken to run program 1000 times: 0.124000


**Implementation/Code:** The following is the code for jacobi

    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <time.h>
    #ifdef _OPENMP
    	#include <omp.h>
    #else
    	#define omp_get_num_threads() 1
    #endif
    
    #define matrix_dimension 25
    
    int n = matrix_dimension;
    float sum;
    
    int main()
    {
      clock_t timestart = clock();
      int numTimeToRun;
     for(numTimeToRun = 0; numTimeToRun < 1000; numTimeToRun++)
     {
      
      float A[n][n];
      float ones[n];
      float x0[n];
      float b[n];
      
      float tol = 0.0001;
      float error = 10.0 * tol;
      float x1[n];
      float res[n];
      int maxiter = 100;
      int iter = 0;
      //
      // create a matrix
      //
      srand((unsigned int)(time(NULL)));
      float a = 5.0;
      
      for(int i=0; i<n; i++)
      {
        for(int j=0; j<n; j++)
        {
          A[i][j] = ((float)rand()/(float)(RAND_MAX) * a);
        }
        x0[i] = ((float)rand()/(float)(RAND_MAX) * a);
      }
      //
      // modify the diagonal entries for diagonal dominance
      // --------------------------------------------------
      //
      for(int i=0; i<n; i++)
      {
        sum = 0.0;
        for(int j=0; j<n; j++)
        {
          sum = sum + fabs(A[i][j]);
        }
        A[i][i] = A[i][i] + sum;
      }
      //
      // generate a vector of ones
      // -------------------------
      //
      #pragma omp parallel num_threads(4)
      { 
      #pragma omp for
      for(int j=0; j<n; j++)
      {
        ones[j] = 1.0;
      }
      //
      // use the vector of ones to generate a right hand side for the testing
      // operation in the code
      // ---------------------
      // 
      #pragma omp critical
      #pragma omp master
      {
      
      for(int i=0; i<n; i++)
      {
        sum = 0.0;
        for(int j=0; j<n; j++)
        {
          sum = sum + A[i][j];
        }
        b[i] = sum;
      }
      }
      
    	
    
    //
    // Jacobi iteration test
    // ---------------------
    //
      
      #pragma omp critical
      #pragma omp master
      {
      for(int i=0; i<n; i++)
      {
        sum = b[i];
        for(int j=0; j<n; j++)
        {
          sum = sum - A[i][j] * x0[i];
        }
        res[i] = sum;
      } 
      }
      
      }
    
      //
      // loop starts here for Jacobi
      // ---------------------------
      //
      while ( error > tol && iter < maxiter) 
      {
        for(int i=0; i<n; i++)
        {
          x1[i] = x0[i] + res[i] / A[i][i];
        }
    	
      
        //
        // compute the error
        // -----------------
        //
        sum = 0.0;
        for(int i=0; i<n; i++)
        {
          float val = x1[i] - x0[i];
          sum = sum + val * val;
        }
    	
      
        error = sqrt(sum);
        //
        // reset the input for the next loop through
        // -----------------------------------------
        //
    
        for(int i=0; i<n; i++)
        {
          x0[i] = x1[i];
        } 
    	
    	
        //
        // compute the next residual
        // -------------------------
        //
    	{
        for(int i=0; i<n; i++)
        {
          sum = b[i];
          for(int j=0; j<n; j++)
          {
            sum = sum - A[i][j] * x0[j];
          }
          res[i] = sum;
        }
    	}
    	
        //
        // update the iteration counter
        // ----------------------------
        //
        iter++;
      //
      // end of loop
      // -----------
      }
      
    {
    for(int i=0; i<n; i++)
       printf("x[%d] = %6f \t res[%d] = %6f\n", i, x1[i], i, res[i]);
    }
    
    
    
     }
    clock_t timestop = clock();
    float deltaTime = ((float)timestop - (float)timestart)/CLOCKS_PER_SEC
    printf("Time taken to run program %d times: %f\n", numTimeToRun, deltaTime);
    
    
    }
   

**Last Modified:** December/2023

<hr>

<a id="power_method"></a>

**Routine Name:**         power_method 

**Author:** Carter Green

**Language:** C. This code can be compiled with the gcc compiler. Note that it can be compiled with the -fopenmp addition.  

**Description/Purpose:** This method iteratively solves for the largest eigenvalue of a matrix using the power method.

**Input:** This method has no inputs.

**Output:** This method will output 7000 eigenvalues of a matrix and then the total time taken to solve all of them.

**Usage/Example:**

        ./power_method.exe

Output from the lines above:

There will be 7000 copies of the following, with different numbers:

        in 20 iterations, eigenvalue = 51.694553
        
Then there will be the following, with a new number each time it is run

        Time taken to run program 7000 times: 0.280999988317489624023437500000


**Implementation/Code:** The following is the code for power_method

        #include <stdio.h>
        #include <stdlib.h>
        #include <math.h>
        #include <time.h>
        
        #ifdef _OPENMP
          #include <omp.h>
        #else
          #define omp_get_num_threads() 0
          #define omp_set_num_threads(int) 0
          #define omp_get_thread_num() 0
        #endif
        
        #define matrix_dimension 10
        
        
        
        int n = matrix_dimension;
        float ynrm;
        
        int main()
        {
          clock_t timestart = clock();
          int numTimeToRun;
        for(numTimeToRun = 0; numTimeToRun < 7000; numTimeToRun++)
         {  
        	
          float A[n][n];
          float v0[n];
          float v1[n];
          float y[n];
          //
          // create a matrix
          //
          srand((unsigned int)(time(NULL)+numTimeToRun));
          float a = 5.0;
          #pragma omp parallel num_threads(4)
          {
          #pragma omp master
          for(int i=0; i<n; i++)
          {
            for(int j=0; j<n; j++)
            {
              A[i][j] = ((float)rand()/(float)(RAND_MAX) * a);
            }
            v0[i] = ((float)rand()/(float)(RAND_MAX) * a);
          }
          //
          // modify the diagonal entries for diagonal dominance
          // --------------------------------------------------
          //
          #pragma omp master
          for(int i=0; i<n; i++)
          {
            float sum = 0.0;
            for(int j=0; j<n; j++)
            {
              sum = sum + fabs(A[i][j]);
            }
            A[i][i] = A[i][i] + sum;
          }
          //
          // generate a vector of ones
          // -------------------------
          //
          #pragma omp for
          for(int j=0; j<n; j++)
          {
            v0[j] = 1.0;
          }
          } //end parallel section
        //
        // power iteration test
        // --------------------
        //
        
        
          float tol = 0.0000001;
          float error = 10.0 * tol;
          float lam1, lam0;
          int maxiter = 100;
          int iter = 0;
          
          
          //
          // loop starts here for the power method
          // -------------------------------------
          //
          while ( error > tol && iter < maxiter ) 
          {
            for(int i=0; i<n; i++)
            {
              y[i] = 0;
              for(int j=0; j<n; j++)
              {
                y[i] = y[i] + A[i][j] * v0[j];
              }
            } 
            //
            // compute the norm of the output
            // ------------------------------
            //
            ynrm = 0.0;
            for(int i=0; i<n; i++)
            {
              ynrm = ynrm + y[i] * y[i];
            }
            ynrm = sqrt(ynrm);
            //
            // normalize the vector from the computation above
            // -----------------------------------------------
            //
        	
            for(int i=0; i<n; i++)
            {
              v1[i] = y[i] / ynrm;
            } 
            //
            // compute the next product
            // ------------------------
            //
            for(int i=0; i<n; i++)
            {
              y[i] = 0.0;
              for(int j=0; j<n; j++)
              {
                y[i] = y[i] + A[i][j] * v1[j];
              }
            }
            //
            // the following is an approximation of the eigenvalue using the Rayleigh
            // rayleigh quotient
            // -----------------
            //
        	
            lam1 = 0.0;
            for(int i=0; i<n; i++)
            {
              lam1 = lam1 + v1[i] * y[i];
            }
            //
            // error computations and updates for the next iteration
            // ------------------------------------------------------
            //
            error = fabs(lam1-lam0);
            lam0 = lam1;
            for(int i=0; i<n; i++)
            {
              v0[i] = v1[i];
            }
            //
            // update the iteration counter
            // ----------------------------
            //
            iter++;
          //
          // end of loop
          // -----------
        	}
          
        
          printf("in %d iterations, eigenvalue = %f\n", iter, lam1);
         }
          
        clock_t timestop = clock();
        float deltaTime = ((float)timestop - (float)timestart)/CLOCKS_PER_SEC;
        printf("Time taken to run program %d times: %-.30f\n",numTimeToRun, deltaTime);
        }	
         

**Last Modified:** December/2023
