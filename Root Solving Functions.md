# Math 4610 Fundamentals of Computational Mathematics Root Solving Functions

**Table of Contents**
1. [Bisection Method](#Bisection Method)
2. [Newton's Method](#Newton's Method)
3. [Secant Method](#Secant Method)
4. [Hybrid Method](#Hybrid Method)
5. [Hybrid Method Known Root](#Hybrid Method Known Root)
6. [LogisticODERoots](#LogisticODERoots)



<hr>
<a id="Bisection Method"></a>

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

   

**Last Modified:** October/2023



<hr>

<a id="Newton's Method"></a>

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


   

**Last Modified:** October/2023




<hr>

<a id="Secant Method"></a>

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

  

**Last Modified:** October/2023




<hr>

<a id="Hybrid Method"></a>

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
            
            

**Last Modified:** October/2023





<hr>

<a id="Hybrid Method Known Root"></a>

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
   

**Last Modified:** October/2023





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








