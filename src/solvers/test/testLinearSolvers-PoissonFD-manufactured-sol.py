import sympy as sym


# Some simple code to get C++ code for exact solution and corresponding source term of Poisson problems in 1D, 2D, and 3D. 
def main():

    # Just uncomment whichever one you want to use.
    #poisson_1D()
    #poisson_2D()
    poisson_3D()
    

def cxx_print(u, f):
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
    cxx_print(u, LU)


# Rotated anisotropic eqn in 2D
def poisson_2D():
    x, y, alpha, beta, gamma, g = sym.symbols('x y alpha beta gamma g')

    # Exact solution
    u = ( sym.sin(2 * sym.pi * x) * sym.sin(4 * sym.pi * y) )

    # Differential operator applied to exact solution
    LU = -alpha*sym.diff( sym.diff(u, x), x) \
         -beta *sym.diff( sym.diff(u, y), y) \
         -gamma*sym.diff( sym.diff(u, x), y)
    
    print("\n\n2D Poisson problem:")
    cxx_print(u, LU)


# 3D equation with anisotropy epsilony in the y direction and epsilonz in the z direction
def poisson_3D():

    x, y, z, epsilony, epsilonz = sym.symbols('x y z epsilony epsilonz')

    # Exact solution
    u = ( sym.sin(2 * sym.pi * x) * sym.sin(4 * sym.pi * y) * sym.sin(6 * sym.pi * z) )

    # Differential operator applied to exact solution
    LU = -sym.diff( sym.diff(u, x), x) \
         -epsilony*sym.diff( sym.diff(u, y), y) \
         -epsilonz*sym.diff( sym.diff(u, z), z)
    
    print("\n\n3D Poisson problem:")
    cxx_print(u, LU)




# Call the main method!
if __name__ == "__main__":
    main()