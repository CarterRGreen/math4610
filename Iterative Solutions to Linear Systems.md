# Math 4610 Fundamentals of Computational Mathematics Iterative Solutions to Linear Systems

**Table of Contents**
1. [jacobiIteration](#jacobiIteration)
2. [jacobiLeslie](#jacobiLeslie)
3. [makeLeslie](#makeLeslie)
4. [gaussSeidelIteration](#gaussSeidelIteration)
5. [testCode](#testCode)


<hr>

<a id="jacobiIteration"></a>

**Routine Name:**  acobiIteration        

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




