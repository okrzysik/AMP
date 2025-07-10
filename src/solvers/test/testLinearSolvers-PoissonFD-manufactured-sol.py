import sympy as sym


# Some simple code to get C++ code for exact solution and corresponding source term of Poisson problems in 1D, 2D, and 3D. 
def main():

    # Just uncomment whichever one you want to use.
    #poisson_1D()
    #poisson_2D()
    #poisson_3D()
    
    diffusionTensor_3D()

def cxx_print_sol_and_rhs(u, f):
    print("-----------------")
    print("exact solution u:")
    print("-----------------")
    print( "u = ", sym.cxxcode( sym.simplify(u) ), ";" , sep = "")
    print("")
    
    print("--------------")
    print("source term f:")
    print("--------------")
    print( "f = ", sym.cxxcode( sym.simplify(f) ), ";", sep = "" )



# eqn in 1D
def poisson_1D():
    x = sym.symbols('x')

    # Exact solution
    u = sym.sin(2 * sym.pi * x) 

    # Differential operator applied to exact solution
    LU = -sym.diff( sym.diff(u, x), x) 
    
    print("\n\n1D Poisson problem:")
    cxx_print_sol_and_rhs(u, LU)


# Rotated anisotropic eqn in 2D
def poisson_2D():
    x, y, alpha, beta, gamma = sym.symbols('x y alpha beta gamma')

    # Exact solution
    u = ( sym.sin(2 * sym.pi * x) * sym.sin(4 * sym.pi * y) )

    # Differential operator applied to exact solution
    LU = -alpha*sym.diff( sym.diff(u, x), x) \
         -beta *sym.diff( sym.diff(u, y), y) \
         -gamma*sym.diff( sym.diff(u, x), y)
    
    print("\n\n2D Poisson problem:")
    cxx_print_sol_and_rhs(u, LU)


# # 3D equation with anisotropy epsilony in the y direction and epsilonz in the z direction
# def poisson_3D():

#     x, y, z, epsilony, epsilonz = sym.symbols('x y z epsilony epsilonz')

#     # Exact solution
#     u = ( sym.sin(2 * sym.pi * x) * sym.sin(4 * sym.pi * y) * sym.sin(6 * sym.pi * z) )

#     # Differential operator applied to exact solution
#     LU = -sym.diff( sym.diff(u, x), x) \
#          -epsilony*sym.diff( sym.diff(u, y), y) \
#          -epsilonz*sym.diff( sym.diff(u, z), z)
    
#     print("\n\n3D Poisson problem:")
#     cxx_print_sol_and_rhs(u, LU)

# 3D equation with anisotropy epsilony in the y direction and epsilonz in the z direction
def poisson_3D():

    x, y, z, alpha, beta, gamma, delta, epsilon, zeta = sym.symbols('x y z alpha beta gamma delta epsilon zeta')

    # Exact solution
    u = ( sym.sin(2 * sym.pi * x) * sym.sin(4 * sym.pi * y) * sym.sin(6 * sym.pi * z) )

    # Differential operator applied to exact solution
    LU = -alpha  *sym.diff( sym.diff(u, x), x) \
         -beta   *sym.diff( sym.diff(u, y), y) \
         -gamma  *sym.diff( sym.diff(u, z), z) \
         -delta  *sym.diff( sym.diff(u, x), y) \
         -epsilon*sym.diff( sym.diff(u, x), z) \
         -zeta   *sym.diff( sym.diff(u, y), z) \
    
    # print("\n\n3D Poisson problem:")
    # cxx_print_sol_and_rhs(u, LU)


def diffusionTensor_3D():

    cphi, sphi, cth, sth, cpsi, spsi, epsy, epsz = sym.symbols('cphi, sphi, cth, sth, cpsi, spsi, epsy, epsz')

    # Just the upper trinagular components
    d11 = epsy*(cphi*cth*spsi + cpsi*sphi)**2 + epsz*spsi**2*sth**2 +\
        (cphi*cpsi - cth*sphi*spsi)**2
    d12 = cpsi*epsz*spsi*sth**2 +\
        epsy*(cphi*cpsi*cth - sphi*spsi)*(cphi*cth*spsi + cpsi*sphi) +\
        (cphi*cpsi - cth*sphi*spsi)*(-cphi*spsi - cpsi*cth*sphi)
    d13 = -cphi*epsy*sth*(cphi*cth*spsi + cpsi*sphi) +\
        cth*epsz*spsi*sth + sphi*sth*(cphi*cpsi - cth*sphi*spsi)
    d22 = cpsi**2*epsz*sth**2 + epsy*(cphi*cpsi*cth - sphi*spsi)**2 +\
        (-cphi*spsi - cpsi*cth*sphi)**2
    d23 = -cphi*epsy*sth*(cphi*cpsi*cth - sphi*spsi) +\
        cpsi*cth*epsz*sth + sphi*sth*(-cphi*spsi - cpsi*cth*sphi)
    d33 = cphi**2*epsy*sth**2 + cth**2*epsz + sphi**2*sth**2

    print( "double d11 = ", sym.cxxcode( sym.simplify(d11) ), ";", sep = "" )
    print( "double d22 = ", sym.cxxcode( sym.simplify(d22) ), ";", sep = "" )
    print( "double d33 = ", sym.cxxcode( sym.simplify(d33) ), ";", sep = "" )
    print( "double d12 = ", sym.cxxcode( sym.simplify(d12) ), ";", sep = "" )
    print( "double d13 = ", sym.cxxcode( sym.simplify(d13) ), ";", sep = "" )
    print( "double d23 = ", sym.cxxcode( sym.simplify(d23) ), ";", sep = "" )




# Call the main method!
if __name__ == "__main__":
    main()