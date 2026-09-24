import sage.all as sg
from blowupsWithSageRewrite import *
from sys import exit

if __name__ == "__main__":
    print("""Welcome to PlaceHolderName CLI! \n 
    This is a command line application which desingularizes projective hypersurfaces. \n
    Copyright (c) 2026 Henry Morin. All Rights Reserved. \n
    ---------------------------------------------------------- \n """)

    print('TBD Add Command List (Enter "exit" at any time to abort the program)')

    m2 = sg.macaulay2

    #Gets the number of dimensions from the user
    dimension = False
    while not dimension:
        dimension = input("Enter number of dimensions: ")
        if dimension == "exit":
            exit("User terminated program")

        try:
            dimension = int(dimension)
            if dimension <= 0:
                print(dimension, "is not a valid input, please enter a positive integer.") 
                dimension = False

        except:
            print(dimension, "is not a valid input, please enter a positive integer.")
            dimension = False

    prime = False

    prime = input("Characteristic of k (0 for QQ, p for GF(p)): ")

    try:
        prime = int(prime)
    except:
        prime = 0

    if not sg.is_prime(prime):
        prime = 0

    mode = False
    while not mode:
        mode = input("Hypersurface or double cover (h/d): ")
        if mode == "exit":
            exit("User terminated program")
        elif (mode != 'd') and (mode != 'h'):
            mode = 'h'

    # Sets N = dimension for convenience
    N = dimension

    # Sets R = QQ[x_0,...,x_N] to be the working ring
    if prime == 0:
        R = sg.PolynomialRing(sg.QQ, 'x', N + 1)
    else:
        R = sg.PolynomialRing(sg.GF(prime), 'x', N + 1)

    F = False
    R_m2 = m2(R)
    #F = m2('x0^2*x2 - x1^2*(x1 + x2)')
    #F = m2('x0^2*x2 - x1^3')
    #F = m2('x0^2*x3 - x1^2*x2')
    #F = m2('x0^7*x1^4 - x2^11')
    #F = m2('x0')
    #F = m2('x0^4*x1^2*x2^2-2*x0^2*x1^4*x2^2+x1^6*x2^2-2*x0^3*x1^2*x2^3+2*x0^2*x1^3*x2^3+2*x0*x1^4*x2^3-2*x1^5*x2^3+x0^2*x1^2*x2^4-2*x0*x1^3*x2^4+x1^4*x2^4-2*x0^5*x1*x2*x3+2*x0^3*x1^3*x2*x3+6*x0^4*x1*x2^2*x3-4*x0^3*x1^2*x2^2*x3-2*x0^2*x1^3*x2^2*x3+2*x0*x1^4*x2^2*x3-2*x1^5*x2^2*x3-6*x0^3*x1*x2^3*x3+6*x0^2*x1^2*x2^3*x3-2*x0*x1^3*x2^3*x3+2*x1^4*x2^3*x3+2*x0^2*x1*x2^4*x3-2*x0*x1^2*x2^4*x3+x0^6*x3^2-4*x0^5*x2*x3^2+2*x0^4*x1*x2*x3^2-2*x0*x1^4*x2*x3^2+6*x0^4*x2^2*x3^2-4*x0^3*x1*x2^2*x3^2+3*x0^2*x1^2*x2^2*x3^2+x1^4*x2^2*x3^2-4*x0^3*x2^3*x3^2+2*x0^2*x1*x2^3*x3^2-2*x0*x1^2*x2^3*x3^2+x0^2*x2^4*x3^2+2*x0^4*x1*x3^3-6*x0^3*x1*x2*x3^3-2*x0^2*x1^2*x2*x3^3+4*x0*x1^3*x2*x3^3+4*x0^2*x1*x2^2*x3^3-2*x0*x1^2*x2^2*x3^3-2*x0^4*x3^4+x0^2*x1^2*x3^4+4*x0^3*x2*x3^4+2*x0^2*x1*x2*x3^4-2*x0*x1^2*x2*x3^4-2*x0^2*x2^2*x3^4-2*x0^2*x1*x3^5+x0^2*x3^6') #Q1
    #F = m2('4*x0^4*x1^2*x2^2+x0^2*x1^4*x2^2-8*x0^3*x1^2*x2^3+4*x0^2*x1^3*x2^3+4*x0^2*x1^2*x2^4-4*x0*x1^3*x2^4+4*x0^5*x1*x2*x3-2*x0^4*x1^2*x2*x3+2*x0^3*x1^3*x2*x3-16*x0^4*x1*x2^2*x3+2*x0^3*x1^2*x2^2*x3-6*x0^2*x1^3*x2^2*x3+20*x0^3*x1*x2^3*x3-8*x0^2*x1^2*x2^3*x3-8*x0^2*x1*x2^4*x3+8*x0*x1^2*x2^4*x3+x0^6*x3^2-2*x0^5*x1*x3^2+x0^4*x1^2*x3^2-6*x0^5*x2*x3^2-4*x0^3*x1^2*x2*x3^2-2*x0^2*x1^3*x2*x3^2+13*x0^4*x2^2*x3^2+10*x0^3*x1*x2^2*x3^2+5*x0^2*x1^2*x2^2*x3^2+2*x0*x1^3*x2^2*x3^2-12*x0^3*x2^3*x3^2-4*x0^2*x1*x2^3*x3^2+4*x0*x1^2*x2^3*x3^2+4*x0^2*x2^4*x3^2-4*x0*x1*x2^4*x3^2-2*x0^5*x3^3+4*x0^4*x1*x3^3-2*x0^3*x1^2*x3^3+10*x0^4*x2*x3^3-6*x0^3*x1*x2*x3^3+8*x0^2*x1^2*x2*x3^3-16*x0^3*x2^2*x3^3+2*x0^2*x1*x2^2*x3^3-6*x0*x1^2*x2^2*x3^3+8*x0^2*x2^3*x3^3-4*x0*x1*x2^3*x3^3+x0^4*x3^4-2*x0^3*x1*x3^4+x0^2*x1^2*x3^4-4*x0^3*x2*x3^4+2*x0^2*x1*x2*x3^4-2*x0*x1^2*x2*x3^4+4*x0^2*x2^2*x3^4+x1^2*x2^2*x3^4') #V2
    #F = m2('x0^4*x1^4-2*x0^4*x1^3*x2+2*x0^3*x1^4*x2+x0^4*x1^2*x2^2-4*x0^3*x1^3*x2^2+x0^2*x1^4*x2^2+2*x0^3*x1^2*x2^3-2*x0^2*x1^3*x2^3+x0^2*x1^2*x2^4+2*x0^4*x1^3*x3-2*x0^4*x1^2*x2*x3+6*x0^3*x1^3*x2*x3-6*x0^3*x1^2*x2^2*x3+4*x0^2*x1^3*x2^2*x3-4*x0^2*x1^2*x2^3*x3+x0^4*x1^2*x3^2-2*x0^3*x1^3*x3^2-6*x0^2*x1^3*x2*x3^2+2*x0^3*x1*x2^2*x3^2-4*x0*x1^3*x2^2*x3^2+6*x0^2*x1*x2^3*x3^2+4*x0*x1*x2^4*x3^2-2*x0^3*x1^2*x3^3-2*x0^3*x1*x2*x3^3+2*x0^2*x1^2*x2*x3^3-2*x0^2*x1*x2^2*x3^3+6*x0*x1^2*x2^2*x3^3-6*x0*x1*x2^3*x3^3+4*x1^2*x2^3*x3^3-4*x1*x2^4*x3^3+x0^2*x1^2*x3^4+4*x0^2*x1*x2*x3^4+x0^2*x2^2*x3^4-2*x0*x1*x2*x3^5-2*x0*x2^2*x3^5+x2^2*x3^6') #H3
    #F = m2('-3*x0^2*x1^4*x2^2+6*x0^2*x1^3*x2^3+2*x0*x1^4*x2^3-3*x0^2*x1^2*x2^4-2*x0*x1^3*x2^4+x1^4*x2^4+2*x0^2*x1^4*x2*x3-6*x0^2*x1^2*x2^3*x3-6*x0*x1^3*x2^3*x3-2*x1^4*x2^3*x3+4*x0^2*x1*x2^4*x3+6*x0*x1^2*x2^4*x3-2*x1^3*x2^4*x3+x0^2*x1^4*x3^2-6*x0^2*x1^3*x2*x3^2-2*x0*x1^4*x2*x3^2+9*x0^2*x1^2*x2^2*x3^2+6*x0*x1^3*x2^2*x3^2+x1^4*x2^2*x3^2-4*x0^2*x1*x2^3*x3^2+4*x1^3*x2^3*x3^2-4*x0*x1*x2^4*x3^2+x1^2*x2^4*x3^2+2*x0*x1^3*x2*x3^3-6*x0*x1^2*x2^2*x3^3-2*x1^3*x2^2*x3^3+4*x0*x1*x2^3*x3^3-2*x1^2*x2^3*x3^3+x1^2*x2^2*x3^4') #Q9
    F = m2('x0*x1*x2*x3*(x0 + x1)*(x1 + x2)*(x2 + x3)*(x3 + x0)') #pg. 68-1 Myers
    #F = m2('x0*x1*x2*x3*(x0 + x1)*(x2 + x3)')

    while not F:
        F = input("Enter the defining polynomial: ")
        if F == "exit":
            exit("User terminated program")

        try:
            F = m2('F')

            if not F.isHomogeneous():
                print(F, "is not a valid input, please enter a homogeneous polynomial.") 
                F = False
        except:
            print(F, "is not a valid input, please enter a homogeneous polynomial.")
            F = False

    #Records affine chart information
    print(R_m2)
    F = F.sage()
    tree = ChartTree(F, N, mode, prime)

    #Records information about the singular locus
    ledger = SingularLedger(F, tree)
    print(ledger)

    if not ledger.components:
        print("F is smooth!")
        exit()
        
    #Asks the user for a blowup until all singularities are resolved
    while ledger.components:
        l = len(ledger.components)
        l_intersection = len(ledger.initial_intersections)
        admissible_idx = ledger.admissible_indices()
        
        valid = False
        i = 0
        while not valid:
            i = input(f'Blowup component (crepant options: {admissible_idx}): ')

            if i == 'exit':
                exit('User terminated program')

            if i == 'detail':
                ledger.print_state()
                continue

            is_intersection_point = False

            try:
                if i[0] =="i":
                    is_intersection_point = True
                    i = i[1:]

                i = int(i)
                invalid_cond = (is_intersection_point and i >= l_intersection) or (not is_intersection_point and i >= l) or (i < 0)
                if invalid_cond:
                    print(f'{i} is not a valid input, please enter an integer from 0 to {l - 1}')
                else:
                    if i not in admissible_idx:
                        confirmation = input("WARNING: The blowup you are about to compute is not admissible, are you sure you want to continue? yes / no:")
                        if confirmation != 'yes':
                            print(f"Blowup at index {i} was not computed.")
                            continue

                    valid = True

            except:
                print(f"{i} is not a valid input, please enter an integer from 0 to {l - 1}")

        blowup_component(ledger, i, is_intersection_point)
        print(ledger)
        print(ledger.tree)

        if ledger.admissible_indices() == []:
            print("No more admissible blowups exist!")

    print("All singularities have been resolved!")
        










