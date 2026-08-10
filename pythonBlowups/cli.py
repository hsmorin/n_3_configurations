import sage.all as sg
from blowupsWithSage import *
from sys import exit

if __name__ == "__main__":
    print("""Welcome to PlaceHolderName CLI! \n 
    This is a comand line application which desingularizes projective hypersurfaces. \n
    Copyright (c) 2026 Henry Morin. All Rights Reserved. \n
    ---------------------------------------------------------- \n """)

    print('TBD Add Comand List (Enter "exit" at any time to abort the program)')

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

    # Sets N = dimension for convince
    N = dimension

    # Sets R = QQ[x_0,...,x_N] to be the working ring
    variables = VarsList("x", N)
    R = m2(f'QQ[{variables.callable()}]')

    F = False
    F = m2('x0^2*x2 - x1^2*(x1 + x2)')
    #F = m2('x0^2*x2 - x1^3')
    F = m2('x0^7*x1^4 - x2^11')
    #F = m2('x0')
    #F = m2('x0^4*x1^2*x2^2-2*x0^2*x1^4*x2^2+x1^6*x2^2-2*x0^3*x1^2*x2^3+2*x0^2*x1^3*x2^3+2*x0*x1^4*x2^3-2*x1^5*x2^3+x0^2*x1^2*x2^4-2*x0*x1^3*x2^4+x1^4*x2^4-2*x0^5*x1*x2*x3+2*x0^3*x1^3*x2*x3+6*x0^4*x1*x2^2*x3-4*x0^3*x1^2*x2^2*x3-2*x0^2*x1^3*x2^2*x3+2*x0*x1^4*x2^2*x3-2*x1^5*x2^2*x3-6*x0^3*x1*x2^3*x3+6*x0^2*x1^2*x2^3*x3-2*x0*x1^3*x2^3*x3+2*x1^4*x2^3*x3+2*x0^2*x1*x2^4*x3-2*x0*x1^2*x2^4*x3+x0^6*x3^2-4*x0^5*x2*x3^2+2*x0^4*x1*x2*x3^2-2*x0*x1^4*x2*x3^2+6*x0^4*x2^2*x3^2-4*x0^3*x1*x2^2*x3^2+3*x0^2*x1^2*x2^2*x3^2+x1^4*x2^2*x3^2-4*x0^3*x2^3*x3^2+2*x0^2*x1*x2^3*x3^2-2*x0*x1^2*x2^3*x3^2+x0^2*x2^4*x3^2+2*x0^4*x1*x3^3-6*x0^3*x1*x2*x3^3-2*x0^2*x1^2*x2*x3^3+4*x0*x1^3*x2*x3^3+4*x0^2*x1*x2^2*x3^3-2*x0*x1^2*x2^2*x3^3-2*x0^4*x3^4+x0^2*x1^2*x3^4+4*x0^3*x2*x3^4+2*x0^2*x1*x2*x3^4-2*x0*x1^2*x2*x3^4-2*x0^2*x2^2*x3^4-2*x0^2*x1*x3^5+x0^2*x3^6')

    while not F:
        F = input("Enter the defining polynomial: ")
        if F == "exit":
            exit("User terminated program")

        try:
            F = m2(f'{F}')

            if ((not m2(f'ring {F.name()} === {R.name()}')) or (not F.isHomogeneous())):
                print(F, "is not a valid input, please enter a homogeneous polynomial.") 
                F = False
        except:
            print(F, "is not a valid input, please enter a homogeneous polynomial.")
            F = False

    #Records affine chart information
    tree = ChartTree(F, N)

    #Records information about the singular locus
    ledger = SingularLedger(F, tree)
    print(ledger)

    if not ledger.components:
        print("F is smooth!")
        exit()
        
    #Asks the user for a blowup until all singularities are resolved
    while ledger.components:
        l = len(ledger.components)
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

            try:
                i = int(i)
                if i >= l or i < 0:
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

        blowup_component(ledger, i)
        print(ledger)
        print(ledger.tree)

        if ledger.admissible_indices() == []:
            print("No more admissible blowups exist!")

    print("All singularities have been resolved!")
        










