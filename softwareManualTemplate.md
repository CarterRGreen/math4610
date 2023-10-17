# Math 4610 Fundamentals of Computational Mathematics Software Manual Template File

**Table of Contents**
1. [smaceps](#smaceps)
2. [dmaceps](#dmaceps)
3. [l2-norm](#l2-norm)
4. [l1-norm](#l1-norm)
5. [linf-norm](#linf-norm)
6. [l2-distance](#l2-distance)
7. [linf-distance](#linf-distance)
8. [linreg](#linreg)
9. [fordiff](#fordiff)
10. [backdiff](#backdiff)
11. [upperTriangular](#upperTriangular)
12. [backSubstitution](#backSubstitution)

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

<a id="upperTriangular"></a>

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

<a id="backSubstitution"></a>

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

<a id=""></a>

**Routine Name:**           

**Author:** Carter Green

**Language:** C. The code can be compiled using the GCC C compiler (gcc).

For example,

    gcc 

will produce an executable **./a.exe** than can be executed. If you want a different name, the following will work a bit
better

    gcc 

**Description/Purpose:** 

**Input:** 

**Output:** 

**Usage/Example:**



Output from the lines above:

  


     

**Implementation/Code:** The following is the code for 
     

**Last Modified:** October/2023




