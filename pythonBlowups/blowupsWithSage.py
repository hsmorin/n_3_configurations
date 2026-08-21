import sage.all as sg
from copy import deepcopy
from math import prod
from itertools import combinations
#The following functions only work if there is already a macaulay2 instance running with a ring QQ[x0,...,xN] and a homogeneous polynomial F initialized.

#TO DO: Str Methods, Error Handeling

#Alias for calling Macaulay2
m2 = sg.macaulay2

#Takes in a list of strings [s0, s1, ..., sn] and returns "s0, s1, ..., sn"
def callable(lst : list[str] | tuple[str, ...]):
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

    def __str__(self):
            return str(self.data)

def id_to_string(id : tuple[int, ...]):
    return "".join([str(i) for i in id])
    
class AffineChart():
    def __init__(self, f, phi, id, blownup_along = False, inclusion = False, pi = False, blowup_ideal = False):
        self.f = f
        self.phi = phi
        self.id = id
        self.blownup_along = blownup_along
        self.inclusion = inclusion
        self.pi = pi
        self.blowup_ideal = blowup_ideal
        self.fraction_field = self.f.ring().frac()
        if inclusion:
            self.blowup_fraction_field = inclusion.source().frac()
        else:
            self.blowup_fraction_field = False

    def __str__(self) -> str:
        f = str(self.f.sage())
        if len(f) <= 20:
            base = f'{id_to_string(self.id)}: ' + f
        else:
            base = f'{id_to_string(self.id)}: f omitted' 

        #Clutters output
        if self.blownup_along:
            base += f" [blown up along {self.blownup_along}]"
        return base

#Takes in a python list of ring elements and outputs a M2 list

def to_m2_list(lst : list):
    if not lst:
        return m2('{}')

    R = lst[0].ring()
    m2(f'use {R.name()}')
    return m2.join(*[m2(f'{{{elm.name()}}}') for elm in lst])

def generate_transition_map(source : AffineChart, target : AffineChart, tree : ChartTree):
    A = source.f.ring()
    A_frac = source.fraction_field
    B = target.f.ring()
    B_frac = target.fraction_field
    N = A.vars().entries().flatten().length().sage()
    source_id = source.id
    target_id = target.id

    trans_map = tree.transition_maps.get((source_id, target_id), False)
    if trans_map:
        return trans_map

    old_trans_map = tree.transition_maps.get((source_id[:-1], target_id), False)
    if old_trans_map:
        frac_blowdown_map = m2(f'map({A_frac.name()}, frac {A_frac.name()}, matrix {source.phi.name()})')
        return frac_blowdown_map * old_trans_map

    old_trans_map = tree.transition_maps.get((source_id, target_id[:-1]), False)

    if old_trans_map: #TO DO: Fix ordering?
        if not (target.pi and target.blownup_along):
            raise ValueError("Cannot find blowup information (blownup_along and pi)")
        B_vars = B.vars().entries().flatten()
        B_to_S_B = m2(f'map({B_frac.name()}, {target.blowup_fraction_field.name()}, {B_vars.name()} | {target.blownup_along.name()})')
        return old_trans_map * B_to_S_B * target.pi

    old_trans_map = tree.transition_maps.get((source_id[:-1], target_id[:-1]), False)
    if old_trans_map:
        if not (source.inclusion and target.inclusion and source.pi and target.pi):
            raise ValueError("Cannot find blowup information (maps pi and iota)")

        if source_id[:-1] == target_id[:-1]:
            return source.inclusion * target.pi
        A_S = source.inclusion
        B_S = target.inclusion
        A_S = A_S.source()
        B_S = B_S.source()

        vars_A_S = [v for v in A_S.vars().entries().flatten()]
        vars_B_S = [v for v in B_S.vars().entries().flatten()]
        t_vars = vars_A_S[N:]

        old_A = old_trans_map.source()
        old_B = old_trans_map.target()
        ones = [m2(f'1_{old_B.name()}') for t in t_vars]
        old_B_to_B_S = m2(f'map({old_B.name()}, {target.blowup_fraction_field.name()}, (flatten entries vars {old_B.name()}) | {to_m2_list(ones).name()})')
        A_S_to_old_A = m2(f'map({source.fraction_field.name()}, {old_A.name()}, take((flatten entries vars {A_S.name()}), {N})')
        lifted_transition_map = A_S_to_old_A * old_trans_map * old_B_to_B_S
        phi = m2(f'map({source.blowup_fraction_field.name()}, {target.blowup_fraction_field.name()}, take((flatten entries matrix {lifted_transition_map.name()}), {N}) | {to_m2_list(t_vars).name()})')
        return source.inclusion * phi * target.pi



#Takes in an ideal in A^n and gives the strict transform in a given chart
def strict_transform(sigma_i, I, blownup_along):
    return sigma_i(I).saturate(sigma_i(blownup_along.ideal()))

def projective_transition_map(A, i : int, j : int):
    m2(f"use {A.name()}")
    frac_A = m2(f'frac {A.name()}')

    if i == j:
        return m2(f'map({frac_A.name()}, {frac_A.name()})')
    A_vars = [v for v in A.vars().entries().flatten()]

    N = len(A_vars)

    inclusion_vars = []
    for k in range(N + 1):
        if k == i:
            inclusion_vars.append(m2(f"1_{A.name()}"))
        elif k > i: 
            inclusion_vars.append(A_vars[k - 1])
        else:
            inclusion_vars.append(A_vars[k])
    
    map_vars = [inclusion_vars[k] / inclusion_vars[j] for k in range(N + 1) if k != j]
    rational_map = m2(f"map({frac_A.name()}, {frac_A.name()}, {to_m2_list(map_vars).name()})")
    return rational_map

#Tree of F in affine charts
class ChartTree():

    def __init__(self, F, N : int):
        new_vars = VarsList("y", N - 1)
        A = m2(f'QQ[{new_vars.callable()}]')
        new_vars = [v for v in A.vars().entries().flatten()]
        R = F.ring()

        self.root = Node(F)
        self.nodes = {(0,) : self.root}
        self.N = N
        self.transition_maps = {}
        self.active_patches = {}

        for i in range(N + 1):
            m2(f'use {A.name()}')
            new_vars = [v for v in A.vars().entries().flatten()]
            chartVars = new_vars[0:i] + [m2(f"1_{A.name()}")] + new_vars[i:]
            id = (0, i)
            phi = m2(f'map({A.name()}, {R.name()}, {to_m2_list(chartVars)})')
            f = phi(F)
            new_node = Node(AffineChart(f, phi, id))
            self.nodes[id] = new_node
            self.root.children.append(new_node)
            self.transition_maps[(id, id)] = projective_transition_map(A, i, i)
            for j in range(i):
                self.transition_maps[(id, (0, j))] = projective_transition_map(A, i, j) 
                self.transition_maps[((0, j), id)] = projective_transition_map(A, j, i)

    #Returns a generator of lists of nodes
    def depth_first_search(self):
        yield from self.root.depth_first_search()

    def find_path_from_id(self, id : str) -> list[Node]:
        path = [self.root]
        for i in id[1:]:
            path.append(path[-1].children[int(i)])
        return path

    def find_sub_path_from_ids(self, start_id, end_id):
        path = [self.nodes[start_id]]
        n = len(start_id)
        if start_id != end_id[:n]:
            raise ValueError("Invalid call of find_sub_path_from_ids()\nPath not valid")

        for i in end_id[n:]:
            path.append(path[-1].children[int(i)])
        return path

    def get_active_patch_from_id(self, id : str):
        O = []
        for i in id:
            if i != 0:
                break

            O.append(0)

        path = self.find_sub_path_from_ids(O, id)
        active_patch = m2("ideal()")

        for chart_node in path:
            chart = chart_node.data
            if chart.id in self.active_patches:
                active_patch = self.active_patches[chart.id]
                continue
            active_patch = chart.phi(active_patch)
            A_vars = chart.f.ring().vars().entries().flatten()
            for i in range(int(chart.id[-1])):
                active_patch.ideal(m2(f"{A_vars[i].name()}"))

            self.active_patches[chart.id] = active_patch
        return active_patch

        

    def __str__(self) -> str:
        output = ""
        for leaf in self.depth_first_search():
            for chart in leaf:
                if isinstance(chart.data, AffineChart):
                    output += str(chart) + "\n"
                else:
                    output += "0\n"
            output += "\n"
        return output

def I_in_chart(leaf, I):
    map = [(chart.data.phi, chart.data.blownup_along) for chart in leaf[1:]]
    for phi, blownup_along in map:
        if blownup_along:
            I = strict_transform(phi, I, blownup_along)
        else: 
            I = phi(I)
    return I

class SubComponent():
    def __init__(self, I, origin : tuple[int, ...], active_charts : list[tuple[int, ...]]):
        self.defining_ideal = I
        self.origin = origin
        self.active_charts = active_charts

    def __str__(self) -> str:
        return id_to_string(self.origin) + ": " + str(self.defining_ideal)

class SingularComponent():
    def __init__(self, sub_components : list[SubComponent], mult_info = None):
        self.sub_components = sub_components
        self.mult_info = mult_info

    def __str__(self) -> str:
        if self.sub_components:
            return "".join([str(sub_comp) + "\n" for sub_comp in self.sub_components])[:-1]
        else:
            return "Empty component"


def order_along_ideal(f, I, max_order = 50):
    #Order of f vanishing along I
    #Largest k such that f in I^k
    fI = f.ideal()
    if not fI.isSubset(I):
        return 0

    if I == m2(f'ideal(1_{I.ring().name()})'):
        raise ValueError("Order_along_ideal: I is the unit ideal")

    k = 1
    Ik = I
    while fI.isSubset(Ik):
        k += 1
        if k > max_order:
            raise RuntimeError(f"order_along_ideal: exceeded max_order = {max_order}")
        Ik = Ik * I
    return k - 1



#Jacobian criterion
def is_smooth_center(I):
    J = I.jacobian()
    R = I.ring()
    minors_ideal = J.ideal()
    sing_locus = minors_ideal.saturate(R.vars().ideal()) #Saturate with the irrelevent_ideal
    return (sing_locus == m2(f'ideal(1_{R.name()})'))

def is_non_principle_linear_center(I):
    gens = I.mingens().entries().flatten()

    if m2(f"1_{I.ring().name()}") in gens:
        return False

    if len(gens) == 1:
        return False

    for f in gens:
        if f.degree().first() != m2("1"):
            return False

    return True

    
    

def discrepancy(f, I):
    #Returns {codim, mult, discrepancy, crepant, linear} for blowing up {f = 0} along center I
    c = I.codim().sage()
    m = order_along_ideal(f, I)
    a = (c - 1) - (m // 2)
    linear = is_non_principle_linear_center(I)
    return {"codim" : c, "mult" : m, "discrepancy" : a, "admissible" : (a == 0), "linear" : linear}

#Class which keeps track of the components of the singular locus
class SingularLedger():
    def __init__(self, F, tree):
        self.tree = tree
        L = F.ideal().singularLocus().ideal()
        L = [comp for comp in L.decompose()]
        self.components = []
        self.initial_intersections = []
        intersections = []

        for i, I in enumerate(L):
            # Computes and lists the intersection of I with each active chart
            active_charts = []
            for chart in tree.depth_first_search():
                id = chart[-1].data.id
                phi_I = I_in_chart(chart, I)
                if phi_I == m2(f'ideal(1_{phi_I.ring().name()})'):
                    continue

                active_charts.append(id)
            mult_info = discrepancy(F, I)
            if active_charts:
                self.components.append(SingularComponent([SubComponent(I, (0,), active_charts)], mult_info))

            if I.dim() == 0:
                continue

            for J in L[0:i]:
                if J.dim() == 0:
                    continue
                intersection = (I + J).radical()

                if intersection == m2(f'ideal(vars {intersection.ring().name()})'):
                    continue

                if intersection == m2(f'ideal(1_{intersection.ring().name()})'):
                    continue

                already_present = False
                for previous_comp in intersections:
                    if m2(f'{intersection.name()} == {previous_comp.name()}'):
                        already_present = True
                        break
                if not already_present:
                    intersections.append(intersection)

        for intersection in intersections:
            intersection_active_charts = []
            for chart in tree.depth_first_search():
                id = chart[-1].data.id
                intersection_in_chart = I_in_chart(chart, intersection)
                if intersection_in_chart == m2(f'ideal(1_{intersection_in_chart.ring().name()})'):
                    continue
                intersection_active_charts.append(id)
            mult_info = discrepancy(F, intersection)
            if intersection_active_charts:
                self.initial_intersections.append(SingularComponent([SubComponent(intersection, (0,), intersection_active_charts)], mult_info))

    def __str__(self) -> str: #Less detailed
        if not self.components:
            return "SingularLedger: no components"

        lines = [f"Singular Ledger: {len(self.components)} component(s)"]
        for i, comp in enumerate(self.components):
            lines.append(f"[{i}] {str(comp)} Admissible: {comp.mult_info["admissible"]}")

        if not self.initial_intersections:
            return "\n".join(lines)

        lines.append("")

        lines.append(f"Intersection Points: {len(self.initial_intersections)} point(s)")
        for i, comp in enumerate(self.initial_intersections):
            lines.append(f"[{i}] {str(comp)} Admissible: {comp.mult_info["admissible"]}")

        return "\n".join(lines)


    #Needs to be fixed
    def print_state(self): #More detailed
        print("Singular components:\n")
        for i, comp in enumerate(self.components):
            info = comp.mult_info
            origins = "".join([id_to_string(sub_comp.origin) + ", " for sub_comp in comp.sub_comps])[:-2]
            ideals = "".join([str(sub_comp.defining_ideal) + ", " for sub_comp in comp.sub_comps])[:-2]
            active_charts = "".join([str(sub_comp.active_charts) + ", " for sub_comp in comp.sub_comps])[:-2]
            print("Component", i)
            print("Initial chart(s):", origins)
            print("Defining ideal(s)", ideals)
            print("Active charts:", active_charts)
            linear_str = f", Linear = {info['linear']}" if info['linear'] is not None else ""
            print(f"  Codim = {info['codim']}, Mult = {info['mult']}, "
                  f"Discrepancy = {info['discrepancy']}, Admissible = {info['admissible']}{linear_str}")
            print()

    def admissible_indices(self):
        return [i for i, c in enumerate(self.components) if c.mult_info['admissible'] and c.mult_info['linear'] == True]

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
        raise ValueError(f"Cannot blowup ideal {I} in chart {id} with less than two generators")

    #Creates ring B = QQ[z0,...,z{n-1},t0,...,t{r-1}]
    Z_vars = VarsList(f'z{id_to_string(id)}_', N - 1) #Avoids repeating variables
    B_vars = VarsList(f't{id_to_string(id)}_', r - 1)
    S = m2(f'QQ[{Z_vars.callable()},{B_vars.callable()}]')
    S_vars = S.vars().entries().flatten() 
    pi_star = m2(f'map({S.name()}, {A.name()}, {{{Z_vars.callable()}}})')
    pi_star_closure = m2(f'i -> {pi_star.name()}(i)')
    gens_in_blowup = gens.apply(pi_star_closure)

    #Compute blowup ideal blowup_ideal = (f_i * tj)
    blowup_ideal = m2('ideal()')
    for i in range(r):
        for j in range(i + 1, r):
            blowup_ideal = blowup_ideal.ideal(S_vars[N + i] * gens_in_blowup[j] - S_vars[N + j] * gens_in_blowup[i])

    #Creates new ring with correct monomialOrder
    E_vars = VarsList(f'e{id_to_string(id)}_', N - 1)
    Et_vars = VarsList(f'et{id_to_string(id)}_', r - 2)

    #Affine_cart_list
    charts = []

    for i in range(r):
        for dep_vars in combinations(E_vars.vars, r - 1):
            ind_vars = [var for var in E_vars.vars if var not in dep_vars] + Et_vars.vars
            Elim_ring = m2(f'QQ[{callable(dep_vars)},{callable(ind_vars)}, MonomialOrder => Eliminate {len(dep_vars)}]')
            Elim_vars = Elim_ring.vars().entries().flatten()

            combination_indices = [int(e.partition('_')[2]) for e in dep_vars]

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
            E_vars_ord = to_m2_list(E_vars_ord) #Converts python list to M2 list
            Elim_vars_ord = m2.join(E_vars_ord, Et_vars_ord)

            substitution = m2.insert(i + N, m2('1'), Elim_vars_ord)
            sub_to_elim = m2(f'map({Elim_ring.name()}, {S.name()}, {substitution.name()})')

            S_vars_list = [S_vars[k] for k in range(N + r)]
            independent_indices = [k for k in range(N) if k not in combination_indices]
            S_vars_reord = [S_vars_list[k] for k in independent_indices]

            S_vars_reord = S_vars_reord + [S_vars_list[N + k] / S_vars_list[N + i] for k in range(r) if k != i]
            frac_S = m2(f'frac {S.name()}')
            S_vars_reord = [m2(f'promote({v.name()}, {frac_S.name()})') for v in S_vars_reord]

            pi_frac = m2(f'map(frac {S.name()}, frac {A.name()}, {to_m2_list(S_vars_reord).name()})')
            m2(f'use {S.name()}')
            m2(f'use {A.name()}')

            #Writes blowup_ideal as a graph over independent variables in the i-th chart
            gens_of_blowup = sub_to_elim(blowup_ideal).gb().gens().entries().flatten()
            dep_vars_elim = Elim_vars.take(r - 1)

            relations = {}
            
            for g in gens_of_blowup:
                sup = m2.select(g.support(), m2(f'x -> member(x, {dep_vars_elim.name()})'))
                if sup.length() == 0:
                    continue
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
                for var in Elim_vars:
                    if var.member(dep_vars_elim):
                        phi_i_sub.append(relations[var])
                    else:
                        phi_i_sub.append(var)

                phi_i_sub = to_m2_list(phi_i_sub)
                phi_i = m2(f'map({Elim_ring.name()}, {Elim_ring.name()}, {phi_i_sub.name()})')
                sub_vars = [m2(f'1_{A.name()}') for i in range(r - 1)]
                sub_vars += [var for var in A_vars]
                sub_i = m2(f'map({A.name()}, {Elim_ring.name()}, {(to_m2_list(sub_vars)).name()})')

                #Chart map A to A
                sigma_i = sub_i * phi_i * pi

                #Inclusion A to S
                iota = sub_i * phi_i * sub_to_elim
                iota_frac = m2(f'map(frac {A.name()}, frac {S.name()}, matrix {iota.name()})')
                m2(f'use {A.name()}')
                m2(f'use {S.name()}')

                Blf = ((strict_transform(sigma_i, f.ideal(), I)).mingens().entries().flatten())[0]
                charts.append(AffineChart(Blf, sigma_i, id + (i,), blownup_along = gens, inclusion = iota_frac, pi = pi_frac, blowup_ideal = blowup_ideal))
                break

        if len(charts) != i + 1:
            raise ValueError("Blowup ideal of I is not a graph")
    return charts

def translate_and_compare(phi, source_ideal, target_ideal):

    frac_source = m2(f'frac {source_ideal.ring().name()}')
    
    gens = source_ideal.mingens().entries().flatten()
    mapped = [phi(m2(f'promote({g.name()}, {frac_source.name()})')) for g in gens]

    numerators = [g.numerator() for g in mapped]
    denominators  = [g.numerator() for g in mapped]

    translated = m2(f'ideal {to_m2_list(numerators).name()}')
    denom_ideal = m2(f'ideal {to_m2_list(denominators).name()}')

    translated = translated.saturate(denom_ideal)

    return translated.radical() == target_ideal.radical()


def same_component(sub_comp : SubComponent, other_comp : SubComponent, tree : ChartTree): #TO DO: Fix
    A_origin = sub_comp.origin
    B_origin = other_comp.origin

    I = sub_comp.defining_ideal
    J = other_comp.defining_ideal

    if A_origin == B_origin:
        return I.radical() == J.radical()

    phi = tree.transition_maps.get((A_origin, B_origin), False)
    if phi:
        return translate_and_compare(phi, I, J)
    else:
        raise KeyError(f"ChartTree has no transition map from {A_origin} to {B_origin}")
    
    return False

#Computes the affine chart tree corrisponding to blowing up component n
def blowup_component(ledger : SingularLedger, i : int, is_intersection_point = False):
    if is_intersection_point:
        singular_component = ledger.initial_intersections[i]
        del ledger.initial_intersections[i]
    else:
        singular_component = ledger.components[i]
        del ledger.components[i]

    tree = ledger.tree
    N = tree.N
    
    for sub_component in singular_component.sub_components:
        I = sub_component.defining_ideal
        origin = sub_component.origin
        print("Origin 1:", id_to_string(origin))
        active_charts = sub_component.active_charts

        for end_id in active_charts:
            leaf = tree.find_sub_path_from_ids(origin, end_id)
            active_chart_node = leaf[-1]
            active_chart = active_chart_node.data
    
            f = active_chart.f
            A = f.ring()
            if f == m2(f'1_{A.name()}'):
                continue

            phi_I = I_in_chart(leaf, I)

            if phi_I == m2(f'ideal(1_{A.name()})'):
                continue

            print(f"Ideal: {phi_I}")
            blowup_charts = blowup_affine_space(active_chart, phi_I)
            for blowup_chart in blowup_charts:
                blowup_chart_node = Node(blowup_chart)
                tree.nodes[blowup_chart.id] = blowup_chart_node
                active_chart_node.children.append(blowup_chart_node)

            for component in ledger.components:
                for sub_component in component.sub_components:
                    J, sub_comp_origin, sub_comp_active_charts = (sub_component.defining_ideal, sub_component.origin, sub_component.active_charts)
    
                    if end_id not in sub_comp_active_charts:
                        continue
                
                    sub_component.active_charts.remove(end_id)
                    sub_leaf = tree.find_sub_path_from_ids(sub_comp_origin, end_id)
                    for blowup_chart_node in active_chart_node.children:
                        blowup_chart = blowup_chart_node.data
                        phi_J = I_in_chart(sub_leaf + [blowup_chart_node], J)
                        if not phi_J == m2(f'ideal(1_{A.name()})'):
                            sub_component.active_charts.append(blowup_chart.id)

            for intersection_component in ledger.initial_intersections:
                for sub_component in intersection_component.sub_components:
                    J, sub_comp_origin, sub_comp_active_charts = (sub_component.defining_ideal, sub_component.origin, sub_component.active_charts)
    
                    if end_id not in sub_comp_active_charts:
                        continue
                
                    sub_component.active_charts.remove(end_id)
                    sub_leaf = tree.find_sub_path_from_ids(sub_comp_origin, end_id)
                    for blowup_chart_node in active_chart_node.children:
                        blowup_chart = blowup_chart_node.data
                        phi_J = I_in_chart(sub_leaf + [blowup_chart_node], J)
                        if not phi_J == m2(f'ideal(1_{A.name()})'):
                            sub_component.active_charts.append(blowup_chart.id)

            print("----------")
            print(f"Blowup of {end_id}:")

            new_singular_sub_components = []

            for k, blowup_chart_node in enumerate(active_chart_node.children):
                blowup_chart = blowup_chart_node.data
                Blf = blowup_chart.f
                print("Chart", k)
                print("f =", Blf.sage())
                singular_locus = Blf.ideal().singularLocus().ideal().decompose()
                E = blowup_chart.phi(phi_I)
                for sub_comp in singular_locus:
                    if E.isSubset(sub_comp):
                        sub_comp_wrapper = SubComponent(sub_comp, blowup_chart.id, [blowup_chart.id])
                        new_singular_sub_components.append(sub_comp_wrapper)

            new_singular_components = []

            #Updates transition maps 
            for source_id, target_id in list(tree.transition_maps):
                if source_id == end_id and target_id == end_id:
                    for k, child in enumerate(active_chart_node.children):
                        for l in range(k + 1):
                            other_child = active_chart_node.children[l]
                            tree.transition_maps[(child.data.id, other_child.data.id)] = generate_transition_map(child.data, other_child.data, tree)
                            tree.transition_maps[(other_child.data.id, child.data.id)] = generate_transition_map(other_child.data, child.data, tree)
                elif source_id == end_id:
                    for child in active_chart_node.children:
                        tree.transition_maps[(child.data.id, target_id)] = generate_transition_map(child.data, tree.nodes[target_id].data, tree)
                elif target_id == end_id:
                    for child in active_chart_node.children:
                        tree.transition_maps[(source_id, child.data.id)] = generate_transition_map(tree.nodes[source_id].data, child.data, tree)
                else:
                    continue


            #Finds new singular components
            for sub_comp in new_singular_sub_components:
                if not new_singular_components:
                    chart = tree.nodes[sub_comp.origin].data
                    new_singular_components.append(SingularComponent([sub_comp], discrepancy(chart.f, sub_comp.defining_ideal)))
                    continue

                for comp in new_singular_components:
                    old_sub_comp = comp.sub_components[0]
                    if same_component(sub_comp, old_sub_comp, tree):
                        comp.sub_components.append(sub_comp)
                        break

                chart = tree.nodes[sub_comp.origin].data
                new_singular_components.append(SingularComponent([sub_comp], discrepancy(chart.f, sub_comp.defining_ideal)))
                
            ledger.components.extend(new_singular_components)










                            





                
        

        
        

