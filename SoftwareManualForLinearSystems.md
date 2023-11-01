# Math 4610 Fundamentals of Computational Mathematics Software Manual Template File

**Table of Contents**
1. [upperTriangular](#upperTriangular)
2. [backSubstitution](#backSubstitution)
3. [LUfactorization](#LUfactorization)
4. [forwardSubstitution](#forwardSubstitution)
5. [forwardSubDiagEquals1](#forwardSubDiagEquals1)
6. [matrixTimesVector](#matrixTimesVector)
7. [makeMatrix](#makeMatrix)
8. [makeVectorOne](#makeVectorOne)
9. [TestingOfMatrixOperations](#TestingOfMatrixOperations)


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

<a id=""></a>

**Routine Name:**          

**Author:** Carter Green

**Language:** python. This code can be imported using import statements  

**Description/Purpose:** 

**Input:** 

**Output:** 

**Usage/Example:**




Output from the lines above:

  


     

**Implementation/Code:** The following is the code for 

   

**Last Modified:** October/2023







