# Math 4610 Fundamentals of Computational Mathematics Open MP Programs

**Table Of Contents**
1. [jacobi](#jacobi)
2. [](#)



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


