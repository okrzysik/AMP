import sympy as sym

def main():

    
    solveForGhosts1D()
    #solveForGhosts2D()


def solveForGhosts1D():
    """The continuous BCs are:
        a1 * u - b1 * cxx * du/dx = r1 @ x = 0
        a2 * u + b2 * cxx * du/dx = r2 @ x = 1

    The boundary lies exactly half way between an interior node and a ghost node. The u on the boundary is the average of this ghost node and the interior node. The derivative on the boundary is approximated by a finite difference of these quantities

    On boundary k, we solve for u^k_{ghost} = alpha^k * u^k_{interior} + beta^k
    """

    print("2D boundaries")
    print("-------------")
    
    # First interior point on west, and east boundaries
    E1,   E2   = sym.symbols("E1,   E2")
    # Ghost point on west, east, south, and north boundaries
    Eg_1, Eg_2 = sym.symbols("Eg_1, Eg_2")
    # PDE coefficients
    h, cxx = sym.symbols("h, cxx")
    # Robin values
    a1, b1, r1  = sym.symbols("d_a1, d_b1, r1")
    a2, b2, r2  = sym.symbols("d_a2, d_b2, r2") 

    # 1
    print("\n-------------")
    print("Boundary 1")
    print("-------------")
    eqn1 = sym.Eq(a1 * (Eg_1 + E1)/2 - b1 * cxx * (E1 - Eg_1)/h, r1 )
    sol1 = sym.solve(eqn1, Eg_1)[0]
    print(f"Solution: Eg = {sol1}")
    sol1_poly = sym.Poly( sol1, E1 )
    coeff_1   = sol1_poly.coeff_monomial(1)  
    coeff_E1  = sol1_poly.coeff_monomial(E1)  
    print( "Constant term: beta1  = {}".format(coeff_1) )
    print( "      E1 term: alpha1 = {}".format(coeff_E1) )

    # 2
    print("\n-------------")
    print("Boundary 2")
    print("-------------")
    eqn2 = sym.Eq(a2 * (E2 + Eg_2)/2 + b2 * cxx * (Eg_2 - E2)/h, r2 )
    sol2 = sym.solve(eqn2, Eg_2)[0]
    print(f"Solution: Eg = {sol2}")
    sol2_poly = sym.Poly( sol2, E2 )
    coeff_2   = sol2_poly.coeff_monomial(1)  
    coeff_E2  = sol2_poly.coeff_monomial(E2)  
    print( "Constant term: beta2  = {}".format(coeff_2) )
    print( "      E2 term: alpha2 = {}".format(coeff_E2) )


def solveForGhosts2D():
    """The continuous BCs are:
        a1 * u - b1 * cxx * du/dx = r1 @ x = 0
        a2 * u + b2 * cxx * du/dx = r2 @ x = 1
        a3 * u - b3 * cyy * du/dy = r3 @ y = 0
        a4 * u + b4 * cyy * du/dy = r4 @ y = 1

    The boundary lies exactly half way between an interior node and a ghost node. The u on the boundary is the average of this ghost node and the interior node. The derivative on the boundary is approximated by a finite difference of these quantities

    On boundary k, we solve for u^k_{ghost} = alpha^k * u^k_{interior} + beta^k
    """

    print("2D boundaries")
    print("-------------")
    
    # First interior point on west, east, south, and north boundaries
    E1j,   E2j,   Ei3,   Ei4   = sym.symbols("E1j,   E2j,   Ei3,   Ei4")
    # Ghost point on west, east, south, and north boundaries
    Eg_1j, Eg_2j, Eg_i3, Eg_i4 = sym.symbols("Eg_1j, Eg_2j, Eg_i3, Eg_i4")
    # PDE coefficients
    h, cxx, cyy = sym.symbols("h, cxx, cyy")
    # Robin values
    a1, b1, r1  = sym.symbols("d_a1, d_b1, r1")
    a2, b2, r2  = sym.symbols("d_a2, d_b2, r2")
    a3, b3, r3  = sym.symbols("d_a3, d_b3, r3")
    a4, b4, r4  = sym.symbols("d_a4, d_b4, r4")   

    # 1
    print("\n-------------")
    print("Boundary 1")
    print("-------------")
    eqn1 = sym.Eq(a1 * (Eg_1j + E1j)/2 - b1 * cxx * (E1j - Eg_1j)/h, r1 )
    sol1 = sym.solve(eqn1, Eg_1j)[0]
    print(f"Solution: Eg = {sol1}")
    sol1_poly = sym.Poly( sol1, E1j )
    coeff_1   = sol1_poly.coeff_monomial(1)  
    coeff_E1j = sol1_poly.coeff_monomial(E1j)  
    print( "Constant  term: beta1   = {}".format(coeff_1) )
    print( "      E1j term: alpha1j = {}".format(coeff_E1j) )

    # 2
    print("\n-------------")
    print("Boundary 2")
    print("-------------")
    eqn2 = sym.Eq(a2 * (E2j + Eg_2j)/2 + b2 * cxx * (Eg_2j - E2j)/h, r2 )
    sol2 = sym.solve(eqn2, Eg_2j)[0]
    print(f"Solution: Eg = {sol2}")
    sol2_poly = sym.Poly( sol2, E2j )
    coeff_2   = sol2_poly.coeff_monomial(1)  
    coeff_E2j = sol2_poly.coeff_monomial(E2j)  
    print( "Constant  term: beta2   = {}".format(coeff_2) )
    print( "      E2j term: alpha2j = {}".format(coeff_E2j) )

    # 3
    print("\n-------------")
    print("Boundary 3")
    print("-------------")
    eqn3 = sym.Eq(a3 * (Eg_i3 + Ei3)/2 - b3 * cyy * (Ei3 - Eg_i3)/h, r3 )
    sol3 = sym.solve(eqn3, Eg_i3)[0]
    print(f"Solution: Eg = {sol3}")
    sol3_poly = sym.Poly( sol3, Ei3 )
    coeff_3   = sol3_poly.coeff_monomial(1)  
    coeff_Ei3 = sol3_poly.coeff_monomial(Ei3)  
    print( "Constant  term: beta3   = {}".format(coeff_3) )
    print( "      Ei3 term: alphai3 = {}".format(coeff_Ei3) )

    # 4
    print("\n-------------")
    print("Boundary 4")
    print("-------------")
    eqn4 = sym.Eq(a4 * (Ei4 + Eg_i4)/2 + b4 * cyy * (Eg_i4 - Ei4)/h, r4 )
    sol4 = sym.solve(eqn4, Eg_i4)[0]
    print(f"Solution: Eg = {sol4}")
    sol4_poly = sym.Poly( sol4, Ei4 )
    coeff_4   = sol4_poly.coeff_monomial(1)  
    coeff_Ei4 = sol4_poly.coeff_monomial(Ei4)  
    print( "Constant  term: beta4   = {}".format(coeff_4) )
    print( "      Ei4 term: alphai4 = {}".format(coeff_Ei4) )

    
















def west():
    """Eliminate Eg, then solve for the first interior DOF E0 in terms of E1 and E2
    
    The continuous BC is:
        a0 * u + b0 * du/dx = r0 @ x = 0
    """

    print("West boundary")
    print("-------------")
    a0, b0, r0, h, d, Eg, E0, E1, E2 = sym.symbols("d_a0, d_b0, d_r0, h, d, Eg, E0, E1, E2")   

    Eg = 3*E0 - 3*E1 + E2
    eqn0 = sym.Eq(a0 * (Eg + E0)/2 + b0 * d * (E0 - Eg)/h, r0 )
    sol0 = sym.solve(eqn0, E0)[0]
    print(f"Solution: E0 = {sol0}")

    sol0_poly = sym.Poly( sol0, E1, E2 )

    coeff_0  = sol0_poly.coeff_monomial(1)  
    coeff_E1 = sol0_poly.coeff_monomial(E1)  
    coeff_E2 = sol0_poly.coeff_monomial(E2)  

    print( "Constant term: beta0  = {}".format(coeff_0) )
    print( "      E1 term: alpha1 = {}".format(coeff_E1) )
    print( "      E2 term: alpha2 = {}".format(coeff_E2) )


def east():
    """Eliminate Eg, then solve for the last interior DOF E3 in terms of E2 and E1.
    
    The continuous BC is:
        a1 * u + b1 * du/dx = r1 @ x = 1
    """

    print("East boundary")
    print("-------------")
    a1, b1, r1, h, d, E1, E2, E3, Eg = sym.symbols("d_a1, d_b1, d_r1, h, d, E1, E2, E3, Eg")   

    Eg = 3*E3 - 3*E2 + E1
    eqn1 = sym.Eq(a1 * (Eg + E3)/2 + b1 * d * (Eg - E3)/h, r1 )
    sol1 = sym.solve(eqn1, E3)[0]
    print(f"Solution: E0 = {sol1}")

    sol0_poly = sym.Poly( sol1, E1, E2 )

    coeff_0  = sol0_poly.coeff_monomial(1)  
    coeff_E1 = sol0_poly.coeff_monomial(E1)  
    coeff_E2 = sol0_poly.coeff_monomial(E2)  

    print( "Constant term: beta1  = {}".format(coeff_0) )
    print( "      E1 term: alpha1 = {}".format(coeff_E1) )
    print( "      E2 term: alpha2 = {}".format(coeff_E2) )

# Call the main method!
if __name__ == "__main__":
    main()