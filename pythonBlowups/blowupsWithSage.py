import sage.all as sg
from copy import deepcopy
from math import prod
#The following functions only work if there is already a macaulay2 instance running with a ring QQ[x0,...,xN] and a homogeneous polynomial F initialized.

#TO DO: Str Methods, Error Handeling

#Alias for calling Macaulay2
m2 = sg.macaulay2

#Takes in a list of strings [s0, s1, ..., sn] and returns "s0, s1, ..., sn"
def callable(lst : list[str]):
    lst_comma = [i + ',' for i in lst]
    string = "".join(lst_comma)
    return string[:-1]
        

#Generates list of strings to name ring variables
class vars_list():
    def __init__(self, l : str, N : int):
        self.vars = [f'{l}{i}' for i in range(N + 1)]
    
    #Returns x0, ..., xN to run QQ(x0, ..., xN)
    def callable(self):
        return callable(self.vars)

#Node Class
class node():
    def __init__(self, data):
        self.data = data 
        self.children = []

    #Returns the path for each leaf under the node
    def depth_first_search(self):
        def recursive_step(current_path):
            current_node = current_path[-1]
            children = getattr(current_node, "children", [])

            if not children:
                yield list(current_path)
            else:
                for child in children:
                    yield from recursive_step(current_path + [child])

        yield from recursive_step([self])


class affine_chart():
    def __init__(self, f, phi, id, blowup_ideal = False):
        self.f = f
        self.phi = phi
        self.blowup_ideal = blowup_ideal
        self.id = id

    def __str__(self) -> str:
        return str(self.f)

#Takes in an ideal in A^n and gives the strict transform in a given chart
def strict_transform(sigma_i, I, blowup_ideal):
    return sigma_i(I).saturate(sigma_i(blowup_ideal))

#Tree of F in affine charts
class chart_tree():
    def __init__(self, F, N : int):
        new_vars = vars_list("y", N - 1)
        A = m2(f'QQ[{new_vars.callable()}]')
        R = F.ring()

        self.root = node(F)
        self.nodes = {"0" : self.root}
        self.N = N
        for i in range(N + 1):
            chartVars = new_vars.vars[0:i] + [str(1)] + new_vars.vars[i:]
            id = '0' + str(i)
            f = F(callable(chartVars))
            phi = m2(f'map({A.name()}, {R.name()}, {{{callable(chartVars)}}})')
            new_node = node(affine_chart(f, phi, id))
            self.nodes[id] = new_node
            self.root.children.append(new_node)

    #Returns a generator of lists of nodes
    def depth_first_search(self):
        yield from self.root.depth_first_search()

def I_in_chart(leaf, I):
    map = [(chart.data.phi, chart.data.blowup_ideal) for chart in leaf[1:]]
    for phi, blowup_ideal in map:
        if blowup_ideal:
            I = strict_transform(phi, I, blowup_ideal)
        else: 
            I = phi(I)
    return I


#Class which keeps track of the components of the singular locus
class singular_ledger():
    def __init__(self, F, tree):
        self.tree = tree
        L = F.ideal().singularLocus().ideal().decompose()
        self.components = []
        for I in L:
            # Computes and lists the intersection of I with each active chart
            active_charts = []
            for chart in tree.depth_first_search():
                id = chart[-1].data.id
                phi_I = I_in_chart(chart, I)
                if phi_I == m2('1'):
                    continue
                else:
                    active_charts.append(id)
            self.components.append((I, "0", active_charts))

    def __str__(self) -> str:
       return ''

#Tests if Z(I) is a graph over A^n where I \in k[A^N] for some N
def graphTest(I, indepVars, slice_ind : int):
    S = I.ring()
    S_vars = S.vars().entries().flatten()
    Chart_vars = S_vars.drop(N + slice_ind) 
    depVars = S_vars.select(m2(f'x -> not member(x, {indepVars.name()}'))

    E = m2(f'QQ[{depVars} | {indepVars}, MonomialOrder => Eliminate {depVars}')
    phi = m2(f'map({E.name()}, {S.name()}, {depVars.name()} | {indepVars.name()}')
    phi_I = phi(I)


    


#Blowups up A at I and returns list of maps A -> A corrisponding to affine charts
def blowup_affine_space(affine_chart, I):

    #Recover information about I
    f = affine_chart.f
    A = f.ring()
    A_vars = A.vars().entries().flatten()
    gens = I.mingens().entries().flatten()
    N = len(A_vars)
    r = len(gens)

    #Creates ring B = QQ[z0,...,z{n-1},t0,...,t{r-1}]
    Z_vars = vars_list('z', N - 1) #Avoids repeating variables
    B_vars = vars_list('t', r - 1)
    S = m2(f'QQ[{Z_vars.callable()},{B_vars.callable()}]')
    S_vars = S.vars().entries().flatten() 
    pi_star = m2(f'map({S.name()}, {A.name()}, {{{Z_vars.callable()}}})')
    pi_star_closure = m2(f'i -> {pi_star.name()}(i)')
    gens = gens.apply(pi_star_closure)

    #Compute blowup ideal blowup_ideal = (f_i * tj)
    blowup_ideal = m2('ideal()')
    for i in range(r):
        for j in range(i + 1, r):
            blowup_ideal = blowup_ideal.ideal(S_vars[N + i] * gens[j] - S_vars[N + j] * gens[i])

    #Creates new ring with correct monomialOrder
    E_vars = vars_list('e', N - 1)
    Et_vars = vars_list('et', r - 2)
    dep_vars = E_vars.vars[N - r + 1:]
    ind_vars = E_vars.vars[:N - r + 1] + Et_vars.vars
    Elim_ring = m2(f'QQ[{callable(dep_vars)}, {callable(ind_vars)}, MonomialOrder => {r - 1}]') 
    Elim_vars = Elim_ring.vars().entries().flatten()
    E_vars = Elim_vars.take(N)
    Et_vars = Elim_vars.drop(N)
    Elim_vars_re_ord = m2.join(E_vars, Et_vars)

    #Blf = (pi_star(f.ideal()) + blowup_ideal).saturate(pi_star(I) + blowup_ideal)
    for i in range(r):
        #Reorders vars to obtain the substitution map E <- S
        Elim_vars_re_ord = m2.join(Elim_vars.take(i), Elim_vars.drop(N + r - 2), (Elim_vars.drop(i)).take(N + r - i - 2))
        substitution = m2.insert(i + N, m2('1'), Elim_vars_re_ord)
        print(substitution)
        sub_to_elim = m2(f'map({Elim_ring.name()}, {S.name()}, {substitution})')

        #Writes blowup_ideal as a graph over independent variables in the i-th chart
        gens = sub_to_elim(blowup_ideal).gb().gens().entries().flatten()

        J = []
        solved = {}
        for g in gens:
            sup = m2.select(g.support(), m2(f'x -> member(x, {Elim_vars.drop(N)})'))
            if sup.length() == 0:
                J.append(g)
                pass
            elif sup.length() == 1:
                var = sup[0]
                if var.degree(g) == 1:
                    coeff = var.coefficient(g)
                    if coeff != 1 and coeff != -1:
                        pass
                    else:
                        fj = var - g if coeff == 1 else var + g
                        if solved.get(var, False):
                            pass
                            raise ValueError("Blowup not a graph in variables ti, xi")
                        else:
                            solved[var] = fj;
                        
        print(solved)


#Computes the affine chart tree corrisponding to blowing up component n
def blowup_component(ledger : singular_ledger, i : int):
    I, origin, charts = ledger.components[i]
    tree = ledger.tree
    N = tree.N
    
    for leaf in tree.nodes[origin].depth_first_search():
        chart = leaf[-1].data
        id = chart.id
        if id not in charts:
            continue
    
        f = chart.f
        A = f.ring()
        phi_I = I_in_chart(leaf, I)

        #Map sigma_i* = sub_i sigma_i iota 
        print(phi_I)
        B = blowup_affine_space(chart, phi_I)

        
        #Add thing to ledger
        #Win???
        
        


    



