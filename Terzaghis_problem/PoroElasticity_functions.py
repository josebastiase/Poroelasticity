import numpy as np
from scipy.optimize import fsolve
from pychebfun import *
from functools import partial

## Mechanical functions

# Compute Biot modulus
def Biot_modulus(n, K, a, K_f):
    Biot_modulus = n / K_f + ((1.0 - a) * (a - n)) / K
    return Biot_modulus

# Compute compressibility
def Confined_compressibility(K, G):
    Confined_compressibility = 1.0 / (K + (4.0 / 3.0) * G)
    return Confined_compressibility

# Compute consolidation coefficient
def Consolidation_coeff(k, a, y, Biot_modulus, Confined_compressibility):
    Consolidation_coeff = k / (y * (Biot_modulus + a**2.0 * Confined_compressibility)) 
    return Consolidation_coeff

# Compute initial condition (given by compression of the pore structure)
def Initial_condition(a, Biot_modulus, Confined_compressibility, q):
    Initial_condition = ((a * Confined_compressibility) / (Biot_modulus + a ** 2.0 * Confined_compressibility)) * q 
    return Initial_condition

# Terzaghi's problem solution
def Terzaghis_problem(c_v, t, h, n_serie = None):
    
    ##--------------------------
    # Terzaghis_problem function:
    # Function to compute classical Terzaghi's problem of consolidation of soils
        # input
            # c_v: Consolidation coefficient (SI units) Float
            # t: Time (SI units) List
            # h: Depth (SI units) Float
            # n_serie: Optional value with the number of eigenvalues computed (the more the better) Integer
        
        # return
            # pressure: Array ([101, length of t])
            # Normalized lenght: lenght discretization / Total lenght
    ##--------------------------

    # Make time an array so we can looped over it later
    t = np.array(t)
    
    # Points where the solution is computed (change this parameter if you want more points un your solution)
    Z = np.linspace(0, h, 101) 
    
    # Empty array to store pressure solution. Rows = length scale, Columns = time 
    pressure = np.zeros((len(Z),len(t))) 
    
    n_serie = n_serie or 100
    
    # Time loop. Time, element form t list. j is used to append results in pressure array
    for j, time in enumerate(t):
        
        # Lenght loop. z is depth discretization, i is used to append results in pressure array
        for i, z in enumerate(Z):
            # Empty array to store sum of the loop below
            D = np.zeros(n_serie - 1) 

            # Compute sum from the inverse Fourier transfor solution. A, B, C are used for simplification.
            # solution is storaged in D, then the sum of the array is obtained.
            for k in range(1, n_serie): 
                A = ((-1.0) ** (k - 1.0)) / (2.0 * k - 1.0)
                B = np.cos((2.0 * k - 1.0) * (np.pi / 2.0) * (z / h))
                C = np.exp( - ((2.0 * k - 1.0) ** 2.0) * ((np.pi ** 2.0) / 4.0) * (( c_v * time) / (h ** 2.0)))
                D[k - 1] = (A * B * C)
            
            # Get solution
            pressure[i, j] = ((4.0 / np.pi) * D.sum())

            
    return pressure, -Z[::-1]/h


# Terzaghi's problem solution two layers (Verruijt's solution)
def Tezaghis_problem_two_layers(c_v, t, h, S, a, m_v, k, n_serie = None):
        
    ##--------------------------
    # Tezaghis_problem_two_layers function:
    # Function to obtain the two-layered solution of Terzaghi's problem. The solution was obtained by 
    # Verruijt, A. (2018). Numerical and analytical solutions of poroelastic problems. 
    # Geotechnical Research, 5(1), 39-50.
        # input
            # c_v: Consolidation coefficient (SI units) Float
            # t: Time (SI units) List
            # h: Depth of each layer (SI units) List
            # S: Biot modulus (SI units) Float
            # a: Biot coefficient (-) Float
            # m_v: Compressivility (SI units) Float
            # k: Hydraulic conductivity (SI units) array
            # n_serie: Optional value with the number of eigenvalues computed (the more the better) Integer
        
        # return
            # pressure: Array ([200, length of t])
            # Normalized lenght: lenght discretization / Total lenght
    ##--------------------------
    
    # System properties
    h = np.array(h)
    
    # Time properties
    t = np.array(t)
    
    # Pore structure properties
    c_v = np.array(c_v)
    k = np.array(k)
    
    # t_1 and t_2
    theta = h ** 2 / c_v
    
    # Compute alpha and betta constant
    m = (S + a**2 * m_v) # Constant value to simplificate expressions
    alpha = (k[1] * m)**.5 / (k[0] * m) ** .5 # Constant value to simplificate expressions
    betta = (theta[0] / theta[1]) ** .5 # Constant value to simplificate expressions
    
    n_serie = n_serie or 100
    
    # Objective function. Here eigenvalues are found. A function is defined, then Fsolve is used to find roots 
    def obj_fun(x):
        params = - alpha * np.sin(betta *x)* np.sin(x) + np.cos(betta *x) * np.cos(x)
        return params

    # Find roots
    X = np.zeros(n_serie) # Vector to store roots
    for i, x in enumerate(np.linspace(1, n_serie, 100)):
        X[i] = fsolve(obj_fun, x)
        
    # Not so efficient, but just to clean up the found values.    
    X = X[X >= 0]
    X = np.around(X, decimals=3)
    X = np.unique(X)

    # Empty arrays to save solutions
    pressure = np.zeros((200,len(t)))
    Z = np.zeros((200,len(t)))
    
    # Points where the solution is computed (change this parameter if you want more points un your solution)
    Z_1 = np.linspace(h[0], 0, 100)
    Z_2 = np.linspace(0, -h[1], 100)
    
    # Time loop. The top layers is labeled as *_1 whereas the second one *_2
    for n, time in enumerate(t):
        
        # Lenght loop. z is depth discretization
        for j, z_1 in enumerate(Z_1):
            F_1 = np.zeros_like(X)
            
            # Compute solution A, B, C, D, E are used to simplify code and F is for storage
            for i, x in enumerate(X):
                A_1 = np.cos(x) * np.cos(betta * x * z_1 / h[0])
                B_1 = alpha * np.sin(x) * np.sin(betta * x * z_1 / h[0])
                C_1 = (1 + alpha * betta) * np.cos(betta * x ) * np.sin(x)
                D_1 = (alpha + betta) * np.sin(betta * x) * np.cos(x)
                E_1 = (np.exp((- (x ** 2) * time) / theta[1])) / x
                F_1[i] = (((A_1-B_1)/(C_1+D_1))*E_1)
         
            pressure[j, n] = 2 * F_1.sum()
            Z[j, n] = z_1

        for j, z_2 in enumerate(Z_2):
            F_2 = np.zeros_like(X)
            
            for i, x in enumerate(X):
                A_2 = np.cos(x) * np.cos(x * z_2 / h[1])
                B_2 = np.sin(x) * np.sin(x * z_2 / h[1])
                C_2 = (1 + alpha * betta) * np.cos(betta * x) * np.sin(x)
                D_2 = (alpha + betta) * np.sin(betta * x) * np.cos(x)
                E_2 = (np.exp((- (x ** 2) * time) / theta[1])) / x
                F_2[i] = (((A_2 - B_2) / (C_2 + D_2)) * E_2)
            
            pressure[j + 100, n] = 2 * F_2.sum()
            Z[j + 100, n] = z_2
            
    return pressure, Z

# Terzaghi's problem solution multilayers
def Terzaghis_problem_multi(c_v, t, h, m, l, a, b, c, init, n_serie = None, space_sol = None):
    
    ##--------------------------
    # Terzaghis_problem_multi function:
    # Function to obtain the multi layers solution of Terzaghi's problem. The solution was obtained by 
    # Hickson, R. I., Barry, S. I., & Mercer, G. N. (2009). Critical times in multilayer diffusion. 
    # Part 1: Exact solutions. International Journal of Heat and Mass Transfer, 52(25-26), 5776-5783.
        # input
            # c_v: Consolidation coefficient (SI units) Float
            # t: Time (SI units) List
            # h: Depth of each layer (SI units) List
            # m: Number of layers Integer
            # l: Lenght of the layers (SI units) List
            # a, b, c : Muxed boundary conditions. Further information can be found in the paper above
            # n_serie: Optional value with the number of eigenvalues computed (the more the better) Integer
            # space_sol: Space discretization
        
        # return
            # pressure: Array ([space_sol, length of t])
            # Normalized lenght: lenght discretization
    ##--------------------------
    
    # Time properties
    t = np.array(t)
    init = np.array(init)
    
    # System properties
    l = np.array(l)
    z = np.zeros(m)
    
    # z stores the SS solution at the right boundaries 
    for i in range(m):
        if i == 0:
            z[0] = l[0]
        else:
            z[i] = z[i - 1] + l[i]
    
    # Diffusion properties. Notation D is used from diffusion coefficient
    D = np.array(c_v)
    D_av = sum(l / D)
    d = np.sqrt(D)
    
    # Boundary conditions
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    
    eigen_val = n_serie or 200
    space_sol = space_sol or m * 100
    
    # Compute steady state solution
    q = np.zeros(m)
    h = np.zeros(m)
    l_D = np.zeros(m)

    q[0] = ((a[0] * c[1] - a[1] * c[0]) * D[-1]) / \
    (a[0] * b[1] * D[0] - a[1] * b[0] * D[-1] + (a[0] * a[1] * D[0] * D[-1]) * D_av)

    h[0] = (c[0] - b[0] * q[0]) / (a[1] + 1E-8) # Compute h_0

    # Compute vectorized h and q
    for i in range(1, m):
        q[i] = D[0] * q[0] / D[i] 
        l_D[i] = l[i - 1] / D[i - 1]
        h[i] = h[0] + D[0] * q[0] * sum(l_D)
        
    # SS solution g    
    def steady_solution(i, x_i, x):
        g = init[i] - (q[i] * (x - x_i) + h[i])

        return g
    
    # SS solution w
    def w_solution(i, x_i, x):
        w = (q[i] * (x - x_i) + h[i])

        return w
    
   
    # Eigenfunction solution length-depend. mu_ is a eigenvalue and mu is the eigenvector
    def eigen_length(i, d, mu_, K1i, K2i, x_i, x): 
        params = K1i * np.sin((mu_ / d[i]) * (x - x_i)) + \
        K2i * np.cos((mu_ / d[i]) * (x - x_i)) 
    
        return params

    
    # Eigenfunction solution time-depend. K1i1 and K2i1 include other functions that are called later on
    def eigen_time(t, mu_):
        params = np.exp(- mu_ ** 2.0 * t)
        
        return params
    
    
    # Define K functions K_11 is always 1
    def K11(mu_):
    
        return 1


    def K21(b, a, d, mu_):
        params = - b[0] * mu_ / (a[0] * d[0])

        return params


    def K1i1(i, d, l, K1i, K2i, mu_):
        params = (d[i - 1] / d[i]) * (K1i(mu_) * np.cos(mu_ * l[i - 1] /\
        d[i - 1]) - K2i(mu_) * np.sin(mu_ * l[i - 1] / d[i - 1]))

        return params


    def K2i1(i, d, l, K1i, K2i, mu_):
        params = K1i(mu_) * (np.sin(mu_ * l[i - 1] / d[i - 1])) + K2i(mu_) * \
        (np.cos(mu_ * l[i - 1] / d[i - 1]))

        return params

    
    # Save K functions in an array. The purpose of this is to find K(mu_)
    K = np.empty((2, m), dtype=type(K2i1)) # Empty array to save the functions
    
    K[0, 0] = K11 # Save first function K_11
    K[1, 0] = partial(K21, b, a, d) # Partial is a function wrapper

    # After the first two are saved the others can be looped using partial
    for i in range(1, m):
        K[0, i] = partial(K1i1, i, d, l)  
        K[1, i] = partial(K2i1, i, d, l)  

    for n in range(1, m):
        K[0, n] = partial(K[0, n], K[0, n - 1], K[1, n - 1])
        K[1, n] = partial(K[1, n], K[0, n - 1], K[1, n - 1])

        
    
    # Function to wrap K_n to get K_n(mu_). Now eigenvalues are obtained
    def eigen_wrapper(K, d, l):

        def eigen(mu_):
            K1n = K[0, K.shape[1] - 1](mu_) # Get k1n to put it in the params function
            K2n = K[1, K.shape[1] - 1](mu_) # Get k2n to put it in the params function
            params = K1n * (a[1] * np.sin(mu_ * l[-1] / d[-1]) + \
            (mu_ * b[1] / d[-1]) * np.cos(mu_ * l[-1] / d[-1])) + \
            K2n * ((-mu_ * b[1] / d[-1]) * np.sin(mu_ * l[-1] / d[-1]) + \
            a[1] * np.cos(mu_ * l[-1] / d[-1]))
            return params

        return eigen
    
    # Get eigenvalues. Chebfunctions are used. Faster and better
    eigen_fun = chebfun(eigen_wrapper(K, d, l), [0, n_serie])
    mu = eigen_fun.roots() # eigenvector
    
    # Compute Sturm–Liouville and store the top solution in SLt and the bottom in SLb
    SLt = np.zeros((len(mu), m))
    SLb = np.zeros((len(mu), m))
    
    # slow part of the code. Get the sum of the integrals 
    for i in range(m):
        if i == 0:
            x_i = h[0]
            
            g_fun = chebfun(partial(steady_solution, i, x_i), [x_i, z[i]])
    
            for j, mu_ in enumerate(mu):
                K1i = K[0, i](mu_)
                K2i = K[1, i](mu_)
                eigen_length_fun = chebfun(partial(eigen_length, i, d, mu_, K1i, K2i, x_i), [x_i, z[i]])

                SLb[j, i] = eigen_length_fun.dot(eigen_length_fun)
                SLt[j, i] = g_fun.dot(eigen_length_fun)
                
        else:
            x_i = z[i - 1]
            g_fun = chebfun(partial(steady_solution, i, x_i), [x_i, z[i]])
            
            for j, mu_ in enumerate(mu):
                K1i = K[0, i](mu_)
                K2i = K[1, i](mu_)
                eigen_length_fun = chebfun(partial(eigen_length, i, d, mu_, K1i, K2i, x_i), [x_i, z[i]])

                SLb[j, i] = eigen_length_fun.dot(eigen_length_fun)
                SLt[j, i] = g_fun.dot(eigen_length_fun)

                

    SLt=SLt.sum(axis=1) 
    SLb=SLb.sum(axis=1)

    x = np.linspace(h[0], z[-1], int(space_sol))
    x_layer = int(space_sol / m)
    pressure = np.zeros((len(x),len(t)))

    #Time loop 
    for v, time in enumerate(t):

        #Layer loop
        for i in range(m):
            if i==0:
                x_i = h[0]
                x_disp = np.linspace(x_i, z[i], x_layer)

                #Lenght loop
                for j, r in enumerate(x_disp):
                    P = np.zeros_like(mu)
                    w = w_solution(i, x_i, r)

                    #Eigenvalues loop
                    for k, mu_ in enumerate(mu):
                        K1i = K[0, i](mu_)
                        K2i = K[1, i](mu_)
                        time_fun = eigen_time(time, mu_)
                        xi = eigen_length(i, d,mu_, K1i, K2i, x_i, r)
                        P[k] = (time_fun * xi * SLt[k] / SLb[k])

                    pressure[j, v] = w + P.sum()

            else:
                x_i = z[i - 1]
                x_disp = np.linspace(x_i, z[i], x_layer)

                #Lenght loop
                for j, r in enumerate(x_disp):
                    P = np.zeros_like(mu)
                    w = w_solution(i, x_i, r)

                    #Eigenvalues loop
                    for k, mu_ in enumerate(mu):
                        K1i = K[0, i](mu_)
                        K2i = K[1, i](mu_)
                        time_fun = eigen_time(time, mu_)
                        xi = eigen_length(i, d, mu_, K1i, K2i, x_i, r)
                        P[k] = (time_fun * xi * SLt[k] / SLb[k])

                    pressure[(i * 100) + j, v] = w + P.sum()
                    
                    
    return pressure, x

# Terzaghi's problem harmonic load
def Terzaghis_harmonic_load(c_v, t, h, t_0, n_serie = None):
    
    ##--------------------------
    # Terzaghis_harmonic_load function:
    # Function to compute Terzaghi's consolidation problem when the media is stressed under harmonic load
        # input
            # c_v: Consolidation coefficient (SI units) Float
            # t: Time (SI units) List
            # h: Depth of each layer (SI units) Float
            # t_0 : Time constante depending of boundary conditions
        
        # return
            # pressure: Array ([101, length of t])
            # Normalized lenght: lenght discretization
    ##--------------------------

    eigen_val = n_serie or 200
    
    # Time properties
    t = np.array(t)
    
    # Consolidation time
    t_c = h**2 / c_v
    
    # System's constant
    alpha = ((np.pi * t_c / t_0) ** .5) / h
    
    # Space discretization
    Z = np.linspace(0, h, 101)
    # Array to save solution
    pressure = np.zeros((len(Z), len(t)))
    
    # Time loop
    for l, time in enumerate(t):
        
        # Space loop
        for i, z in enumerate(Z):
            comp = ((np.cosh((1 + 1j) * alpha * z)) /\
                    (np.cosh((1 + 1j) * alpha * h))) - 1 # variable to easy code
            
            # Get the solution
            SteadySolution = - 0.5 * comp.imag * np.sin(2 * np.pi * time / t_0) \
            + 0.5 * comp.real * np.cos(2 * np.pi * time / t_0)
            eigen = np.arange(eigen_val)
            
            La = [] # Store the solution in a list (is so simple that speed is no concern)
            for k in eigen:
                A = (-1) ** k / (1 + 2 * k)
                B = np.cos((1 + 2 * k) * (np.pi * z / (2 * h)))
                C = np.exp(-(1 + 2 * k) ** 2 * np.pi ** 2 * time / (4 * t_c))
                D = (1 + (1 + 2 * k) ** 4 * np.pi ** 2 * t_0 ** 2 / (64 * t_c **2))
                La.append(A * B * C / D)

            La = np.array(La)
            
            pressure[i, l] = SteadySolution + (2 / np.pi) * La.sum()

    return pressure, Z / h
