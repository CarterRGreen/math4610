# Math 4610 Fundamentals of Computational Mathematics Root Solving Functions

**Table of Contents**
1. [BisectionMethod](#BisectionMethod)
2. [Newton'sMethod](#Newton'sMethod)
3. [SecantMethod](#SecantMethod)
4. [HybridMethod](#HybridMethod)
5. [HybridMethodKnownRoot](#HybridMethodKnownRoot)
6. [LogisticODERoots](#LogisticODERoots)
7. [TestRoutine](#TestRoutine)



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



