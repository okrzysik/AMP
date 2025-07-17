import sympy as sym

"""

"""


def main():

    # Just uncomment whichever one you want to use.
    #heat_1D()

    heat_2D()
    
    

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
def heat_1D():
    x, t, cxx = sym.symbols('x, t, cxx')

    # Exact solution
    u = sym.sin(1.5 * sym.pi * x) * sym.cos( 2 * sym.pi * t )

    # Differential operator applied to exact solution
    LU = sym.diff(u, t) - sym.diff( cxx*sym.diff(u, x), x) 
    
    print("\n\n1D Heat problem:")
    cxx_print_sol_and_rhs(u, LU)


    print("-----------------------------")
    print("exact solution gradient dudx:")
    print("-----------------------------")
    print( "dudx = ", sym.cxxcode( sym.simplify( sym.diff(u, x) ) ), ";", sep = "" )


# eqn in 2D
def heat_2D():
    x, y, t, cxx, cyy = sym.symbols('x, y, t, cxx, cyy')

    # Exact solution
    u = sym.sin(1.5 * sym.pi * x) * sym.sin(1.5 * sym.pi * y) * sym.cos( 2 * sym.pi * t )

    # Differential operator applied to exact solution
    LU = sym.diff(u, t) - sym.diff( cxx*sym.diff(u, x), x) - sym.diff( cyy*sym.diff(u, y), y) 
    
    print("\n\n2D Heat problem:")
    cxx_print_sol_and_rhs(u, LU)


    print("-----------------------------")
    print("exact solution gradient dudx:")
    print("-----------------------------")
    print( "dudx = ", sym.cxxcode( sym.simplify( sym.diff(u, x) ) ), ";", sep = "" )

    print("-----------------------------")
    print("exact solution gradient dudy:")
    print("-----------------------------")
    print( "dudy = ", sym.cxxcode( sym.simplify( sym.diff(u, y) ) ), ";", sep = "" )



# Call the main method!
if __name__ == "__main__":
    main()