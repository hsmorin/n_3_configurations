import sage.all as sg
from copy import deepcopy
from math import prod
from itertools import combinations
#The following functions only work if there is already a macaulay2 instance running with a ring QQ[x0,...,xN] and a homogeneous polynomial F initialized.

#TO DO: Str Methods, Error Handeling

#Alias for calling Macaulay2
m2 = sg.macaulay2

#Takes in a list of strings [s0, s1, ..., sn] and returns "s0, s1, ..., sn"
def callable(lst : list[str] | tuple[str]):
    lst_comma = [i + ',' for i in lst]
    string = "".join(lst_comma)
    return string[:-1]
        

#Generates list of strings to name ring variables
class VarsList():
    def __init__(self, l : str, N : int):
        self.vars = [f'{l}{i}' for i in range(N + 1)]
    
    #Returns x0, ..., xN to run QQ(x0, ..., xN)
    def callable(self):
        return callable(self.vars)

#Node Class
class Node():
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


class AffineChart():
    def __init__(self, f, phi, id, blowup_ideal = False):
        self.f = f
        self.phi = phi
        self.blowup_ideal = blowup_ideal
        self.id = id

    def __str__(self) -> str:
        return self.id + ": " + str(self.f)

#Takes in an ideal in A^n and gives the strict transform in a given chart
def strict_transform(sigma_i, I, blowup_ideal):
    return sigma_i(I).saturate(sigma_i(blowup_ideal))

#Tree of F in affine charts
class ChartTree():
    def __init__(self, F, N : int):
        new_vars = VarsList("y", N - 1)
        A = m2(f'QQ[{new_vars.callable()}]')
        R = F.ring()

        self.root = Node(F)
        self.nodes = {"0" : self.root}
        self.N = N
        for i in range(N + 1):
            chartVars = new_vars.vars[0:i] + [str(1)] + new_vars.vars[i:]
            id = '0' + str(i)
            f = F(callable(chartVars))
            phi = m2(f'map({A.name()}, {R.name()}, {{{callable(chartVars)}}})')
            new_node = Node(AffineChart(f, phi, id))
            self.nodes[id] = new_node
            self.root.children.append(new_node)

    #Returns a generator of lists of nodes
    def depth_first_search(self):
        yield from self.root.depth_first_search()

    def find_path_from_id(self, id : str) -> list[Node]:
        path = [self.root]
        for i in id[1:]:
            path.append(path[-1].children[i])
        return path


    def __str__(self) -> str:
        output = ""
        for leaf in self.depth_first_search():
            for chart in leaf:
                chart = chart.data
                if type(chart) == AffineChart:
                    output += chart.id + ": " + str(chart.f.sage()) + "\n"
                else:
                    output += "0\n"
            output += "\n"
        return output

def I_in_chart(leaf, I):
    map = [(chart.data.phi, chart.data.blowup_ideal) for chart in leaf[1:]]
    for phi, blowup_ideal in map:
        if blowup_ideal:
            I = strict_transform(phi, I, blowup_ideal)
        else: 
            I = phi(I)
    return I


#Class which keeps track of the components of the singular locus
class SingularLedger():
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

    def print_state(self):
        print("Singular components:\n")
        for i, component in enumerate(self.components):
            print("Component", i)
            print("Initial chart", component[1])
            print("Defining ideal", component[0])
            print("Active charts:", *component[2])


#Takes in a python list and outputs a M2 list
def to_m2_list(lst : list):
    return m2.join(*[m2(f'{{{elm.name()}}}') for elm in lst])

#Blowups up A at I and returns list of maps A -> A corrisponding to affine charts
def blowup_affine_space(chart, I):

    #Recover information about I
    f = chart.f
    id = chart.id
    A = f.ring()
    A_vars = A.vars().entries().flatten()
    gens = I.mingens().entries().flatten()
    N = len(A_vars)
    r = len(gens)

    if r <= 1:
        return [AffineChart(f, m2(f'id_{A.name()}'), id + '0')]

    #Creates ring B = QQ[z0,...,z{n-1},t0,...,t{r-1}]
    Z_vars = VarsList('z', N - 1) #Avoids repeating variables
    B_vars = VarsList('t', r - 1)
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
    E_vars = VarsList('e', N - 1)
    Et_vars = VarsList('et', r - 2)

    #Affine_cart_list
    charts = []

    for i in range(r):
        for dep_vars in combinations(E_vars.vars, r - 1):
            ind_vars = [var for var in E_vars.vars if var not in dep_vars] + Et_vars.vars
            Elim_ring = m2(f'QQ[{callable(dep_vars)},{callable(ind_vars)}, MonomialOrder => Eliminate {len(dep_vars)}]')
            Elim_vars = Elim_ring.vars().entries().flatten()

            combination_indices = [int(E[1:]) for E in dep_vars]
            E_vars_ord = []

            num_dep = 0
            num_ind = 0
            for j in range(N):
                if j in combination_indices:
                    E_vars_ord.append(Elim_vars[num_dep])
                    num_dep += 1
                else:
                    E_vars_ord.append(Elim_vars[r - 1 + num_ind])
                    num_ind += 1

            Et_vars_ord = Elim_vars.drop(N)
            E_vars_ord = m2.join(*[m2(f'{{{var.name()}}}') for var in E_vars_ord]) #Converts python list to M2 list
            Elim_vars_ord = m2.join(E_vars_ord, Et_vars_ord)
            substitution = m2.insert(i + N, m2('1'), Elim_vars_ord)
            sub_to_elim = m2(f'map({Elim_ring.name()}, {S.name()}, {substitution.name()})')

            #Writes blowup_ideal as a graph over independent variables in the i-th chart
            gens = sub_to_elim(blowup_ideal).gb().gens().entries().flatten()
            dep_vars_elim = Elim_vars.take(r - 1)

            relations = {}
            
            for g in gens:
                sup = m2.select(g.support(), m2(f'x -> member(x, {dep_vars_elim.name()})'))
                if sup.length() == 0:
                    break
                elif sup.length() == 1:

                    var = sup[0]
                    if var.degree(g) != 1:
                        break

                    coeff = var.coefficient(g)
                    if coeff == 1: 
                        fj = var - g
                    elif coeff == -1:
                        fj = var + g
                    else:
                        break
                    relations[var] = fj;
                else:
                    break

            if len(relations) == r - 1:
                #Generate Affine Chart
                pi = m2(f'map({Elim_ring.name()}, {A.name()}, {E_vars_ord.name()})')
                phi_i_sub = []
                for var in Elim_vars_ord:
                    if var.member(dep_vars_elim):
                        phi_i_sub.append(relations[var])
                    else:
                        phi_i_sub.append(var)

                phi_i_sub = to_m2_list(phi_i_sub)
                phi_i = m2(f'map({Elim_ring.name()}, {Elim_ring.name()}, {phi_i_sub.name()})')
                sub_vars = []
                j = 0
                for var in Elim_vars_ord:
                    if var.member(dep_vars_elim):
                        sub_vars.append(m2('1'))
                    else:
                        sub_vars.append(A_vars[j])
                        j += 1
                sub_i = m2(f'map({A.name()}, {Elim_ring.name()}, {(to_m2_list(sub_vars)).name()})')
                sigma_i = sub_i * phi_i * pi
                Blf = ((strict_transform(sigma_i, f.ideal(), I)).mingens().entries().flatten())[0]
                charts.append(AffineChart(Blf, sigma_i, id + str(i), blowup_ideal = I))
                break

        if len(charts) != i + 1:
            raise ValueError("Blowup ideal of I is not a graph")

    return charts

        
#Computes the affine chart tree corrisponding to blowing up component n
def blowup_component(ledger : SingularLedger, i : int):
    I, origin, active_charts = ledger.components[i]
    del ledger.components[i]
    tree = ledger.tree
    N = tree.N
    
    for leaf in tree.nodes[origin].depth_first_search():
        active_chart_node = leaf[-1]
        active_chart = active_chart_node.data
        id = active_chart.id
        if id not in active_charts:
            continue
    
        f = active_chart.f
        A = f.ring()
        phi_I = I_in_chart(leaf, I)

        blowup_charts = blowup_affine_space(active_chart, phi_I)
        for blowup_chart in blowup_charts:
            blowup_chart_node = Node(blowup_chart)
            tree.nodes[blowup_chart.id] = blowup_chart_node
            active_chart_node.children.append(blowup_chart_node)

        for j, component in enumerate(ledger.components):
            J, comp_origin, comp_active_charts = component

            if id not in comp_active_charts:
                continue
            
            component[2].remove(id)
            for blowup_chart_node in active_chart_node.children:
                blowup_chart = blowup_chart_node.data
                for sub_leaf in tree.nodes[comp_origin].depth_first_search():
                    if sub_leaf[-2].data.id != id:
                        continue
                    phi_J = I_in_chart(sub_leaf, J)
                    if not phi_J == m2(f'1_{A.name()}'):
                        component[2].append(blowup_chart.id)


        print("----------")
        print(f"Blowup of {id}:")
        for k, blowup_chart_node in enumerate(active_chart_node.children):
            blowup_chart = blowup_chart_node.data
            print("Chart", k)
            print("f =", blowup_chart.f.sage())
            K = blowup_chart.f.ideal().singularLocus().ideal().decompose()
            E = I_in_chart(leaf + [blowup_chart_node], I)
            for comp in K:
                if E.isSubset(comp):
                    ledger.components.append((comp, blowup_chart.id, [blowup_chart.id]))









                            





                
        

        
        

