import sage.all as sg
from itertools import combinations
from sage.rings.polynomial.term_order import TermOrder

#The following functions only work if there is already a Macaulay2 instance running with a ring QQ[x0,...,xN] and a homogeneous polynomial F initialized.

#TO DO: Str Methods, Error Handling

#Alias for calling Macaulay2
m2 = sg.macaulay2
m2('needsPackage "RationalMaps"')

def ideal_to_string(I):
    G = I.gens()

    string = "Ideal ("
    for g in G:
        string += str(g) + ", "

    string = string[:-2] + ")"
    return string


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
    return "".join([str(i) for i in id])[:-1]
    
class AffineChart():
    def __init__(self, f, sigma, id, blowup_center = None, inclusion = None, restriction = None, blowup_ideal = None, tau = None, tau_inv = None, sigma_inv = None):
        self.f = f
        self.sigma = sigma
        self.sigma_inv = sigma_inv
        self.id = id
        self.blowup_center = blowup_center
        self.inclusion = inclusion
        self.restriction = restriction
        self.blowup_ideal = blowup_ideal
        self.tau = tau
        self.tau_inv = tau_inv
        self.skipped = False

    def __str__(self) -> str:
        f_str = str(self.f)
        if len(f_str) <= 40:
            base = f'{id_to_string(self.id)}: ' + f_str
        else:
            base = f'{id_to_string(self.id)}: f omitted'

        #Clutters output
        if self.blowup_center is not None:
            base += f" [Blownup along: {ideal_to_string(self.blowup_center)}]"
        return base

#Takes in a python list of ring elements and outputs a M2 list

def to_m2_list(lst : list):
    if not lst:
        return m2('{}')
    R = lst[0].ring()
    return m2.join(*[m2(f'{{{elm.name()}}}') for elm in lst])

def exceptional_divisor(chart :AffineChart):
    return chart.sigma(chart.blowup_center) 

def generate_transition_map(source : AffineChart, target : AffineChart, tree : ChartTree):
    A = source.f.parent()
    frac_A = A.fraction_field()
    B = target.f.parent()
    frac_B = B.fraction_field()
    N = len(A.gens())
    source_id = source.id
    target_id = target.id

    #Map already computed
    trans_map = tree.transition_maps.get((source_id, target_id), False)
    if trans_map:
        return trans_map

    #Source[:-1] to target map exists
    old_trans_map = tree.transition_maps.get((source_id[:-1], target_id), False)
    if old_trans_map:
        tau = tree.nodes[source_id[:-1]].data.tau
        frac_tau = tau.extend_to_fraction_field()
        frac_blowdown_map = source.sigma.extend_to_fraction_field()
        return frac_blowdown_map * frac_tau * old_trans_map

    #Source to target[:-1] exists
    old_trans_map = tree.transition_maps.get((source_id, target_id[:-1]), False)

    #TO DO: Test the remaining maps in the blowup (i.e. inclusion, pi, and restriction, and fix the remainder of this method)
    if old_trans_map:
        if target.sigma_inv is None:
            raise ValueError("Cannot find blowup information (blownup_along and pi)")
        tau_inv = tree.nodes[target_id[:-1]].data.tau_inv
        frac_tau_inv = tau_inv.extend_to_fraction_field()
        frac_sigma_inv = target.sigma_inv.extend_to_fraction_field()
        return old_trans_map * frac_tau_inv * frac_sigma_inv

    #Source[:-1] to target[:-1] exists
    old_trans_map = tree.transition_maps.get((source_id[:-1], target_id[:-1]), False)
    if old_trans_map:
        if (source.inclusion is None or source.restriction is None or target.inclusion is None or target.restriction is None):
            raise ValueError("Cannot find blowup information (maps pi and iota)")

        if source_id[:-1] == target_id[:-1]:
            return source.inclusion * target.restriction

        tau = target_id[:-1].data.tau
        tau_inv = source_id[:-1].data.tau_inv
        frac_tau = tau.extend_to_fraction_field()
        frac_tau_inv = tau_inv.extend_to_fraction_field()

        phi = frac_tau_inv * old_trans_map * frac_tau

        #Extracting ring info
        frac_S_i = source.restriction.target()
        frac_S_j = target.restriction.target()
        vars_S_i = frac_S_i.gens()
        x_vars_S_i = vars_S_i[:N]
        frac_A_i = phi.target()
        sub_map = frac_A_i.hom(x_vars_S_i, vars_S_i)

        #Gets vars for phi_tilde
        phi_tilde_vars = (sub_map * phi).im_gens() + vars_S_i[N:]

        #Constructs phi_tilde from phi such that phi_tilde(x,t) = phi(x,t)
        phi_tilde = frac_S_j.hom(phi_tilde_vars, frac_S_i)

        return source.inclusion * phi_tilde * target.restriction

#Takes in an ideal in A^n and gives the strict transform in a given chart
def strict_transform(sigma_i, I, blowup_center):
    R = blowup_center.ring()
    return sigma_i(I).saturation(sigma_i(R.ideal(blowup_center)))[0]

#Transition map (i,j) is the map from spec A_i to spec A_j
def projective_transition_map(A, i : int, j : int):
    frac_A = A.fraction_field()

    if i == j:
        return frac_A.Hom(frac_A).identity()

    vars_A = frac_A.gens()

    N = len(vars_A)

    inclusion_vars = vars_A[:i] + (frac_A(1),) + vars_A[i:]
    map_vars = [inclusion_vars[k] / inclusion_vars[j] for k in range(N + 1) if k != j]
    rational_map = frac_A.hom(map_vars, frac_A)
    return rational_map

#Tree of F in affine charts
class ChartTree():
    def __init__(self, F, N : int, mode : str = 'h', prime : int = 0):
        
        if prime == 0:
            A = sg.PolynomialRing(sg.QQ, 'x', N)
        else:
            A = sg.PolynomialRing(sg.GF(prime), 'x', N)

        vars_A = A.gens()
        R = F.parent()

        self.mode = mode
        self.root = Node(F)
        self.nodes = {(0,) : self.root}
        self.N = N
        self.transition_maps = {}
        self.active_patches = {}

        for i in range(N + 1):
            chartVars = vars_A[0:i] + (A(1),) + vars_A[i:]
            id = (0, i)
            phi = R.hom(chartVars, A)
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

    #TO DO: Fix
    def get_active_patch_from_id(self, id : str):
        O = (0,)
        for i in id:
            if i != 0:
                break

            O = O + (0,)

        path = self.find_sub_path_from_ids(O, id)
        active_patch = m2("ideal()")

        for chart_node in path:
            chart = chart_node.data
            if chart.id in self.active_patches:
                active_patch = self.active_patches[chart.id]
                continue
            active_patch = chart.sigma(active_patch)
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

#TO DO: Test

def I_in_chart(leaf, I):

    if isinstance(leaf[0].data, AffineChart):
        tau = leaf[0].data.tau
    else:
        tau = None

    if tau is not None:
        I = tau(I)
    
    for chart in leaf[1:]:
        chart = chart.data
        sigma, blowup_center, tau = chart.sigma, chart.blowup_center, chart.tau
        

        if blowup_center is not None:
            I = strict_transform(sigma, I, blowup_center)
        else: 
            I = sigma(I)

        if tau is not None:
            I = tau(I)
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
        self.blownup = False

    def __str__(self) -> str:
        if self.sub_components:
            return "".join([str(sub_comp) + "\n" for sub_comp in self.sub_components])[:-1]
        else:
            return "Empty component"


def order_along_ideal(f, I, max_order = 50):
    #Order of f vanishing along I
    #Largest k such that f in I^k
    R = f.parent()
    fI = R.ideal(f)
    if not fI <= I:
        return 0

    if I == R.ideal(1):
        raise ValueError("Order_along_ideal: I is the unit ideal")

    k = 1
    Ik = I
    while fI <= Ik:
        k += 1
        if k > max_order:
            raise RuntimeError(f"order_along_ideal: exceeded max_order = {max_order}")
        Ik = Ik * I
    return k - 1

def codim(I):
    R = I.ring()
    return R.ngens() - I.dimension()

#Jacobian criterion
def is_smooth_center(I):
    I = m2(I)
    R = I.ring()
    singular_locus = I.singularLocus().ideal()
    return singular_locus == m2(f"ideal 1_{R.name()}")

def is_non_principal_linear_center(minimal_generators : list[sg.polynomial]) -> bool:

    if not minimal_generators:
        return False

    R = minimal_generators[0].parent()

    if R("1") in minimal_generators:
        return False

    if len(minimal_generators) == 1:
        return False

    for f in minimal_generators:
        if f.degree() > 1:
            return False

    return True

def is_non_principal_complete_intersection(generators : list[sg.polynomial]) -> bool:
    r = len(generators)

    if r <= 1:
        return False

    I = sg.ideal(generators)
    
    return (codim(I) == r)

def is_non_principal_quadratic_center(minimal_generators : list[sg.poylnomial]) -> bool:
    if not minimal_generators:
        return False

    R = minimal_generators[0].parent()

    if R("1") in minimal_generators:
        return False

    if len(minimal_generators) == 1:
        return False

    is_quadratic = False
    for f in minimal_generators:
        d = f.degree()
        if d == 2:
            is_quadratic = True
            break
        if d > 3:
            return False

    return is_quadratic
    
    

def discrepancy_hyperplane(f, I):
    #Returns {codim, mult, discrepancy, crepant, linear} for blowing up {f = 0} along center I
    R = f.parent()
    N = len(R.gens())
    c = N - I.dimension()
    m = order_along_ideal(f, I)
    a = c - 1 - m
    return {"codim" : c, "mult" : m, "discrepancy" : a, "admissible" : (a == 0)}

def discrepancy_double_cover(f, I):
    #Returns {codim, mult, discrepancy, crepant, linear} for blowing up {f = 0} along center I
    R = f.parent()
    N = len(R.gens())
    c = N - I.dimension()
    m = order_along_ideal(f, I)
    a = (c - 1) - (m // 2)
    return {"codim" : c, "mult" : m, "discrepancy" : a, "admissible" : (a == 0)}

#Class which keeps track of the components of the singular locus
class SingularLedger():
    def __init__(self, F, tree):
        self.tree = tree
        R = F.parent()
        J = R.ideal(F.gradient())
        L = J.minimal_associated_primes()
        self.components = []
        self.initial_intersections = []
        intersections = []
        dimensions = []
        filtered_L = []
        for I in L:
            dim = I.dimension() - 1
            if dim <= -1:
                continue
            
            filtered_L.append(I)
            dimensions.append(dim)

        L = filtered_L
         

        # Computes and lists the intersection of I with each active chart
        for i, I in enumerate(L):
            active_charts = []
            for leaf in tree.depth_first_search():
                chart = leaf[-1].data
                A = chart.f.parent()
                id = chart.id
                sigma_I = I_in_chart(leaf, I)
                if sigma_I == 1:
                    continue

                active_charts.append(id)

            if tree.mode == 'h':
                mult_info = discrepancy_hyperplane(F, I)
            else:
                mult_info = discrepancy_double_cover(F, I)

            if active_charts:
                self.components.append(SingularComponent([SubComponent(I, (0,), active_charts)], mult_info))

            if dimensions[i] <= 0:
                continue

            for j, J in enumerate(L[0:i]):
                if dimensions[j] <= 0:
                    continue

                intersection = (I + J).radical()

                if intersection == R.ideal(R.gens()):
                    continue

                if intersection == R.ideal(1):
                    continue

                already_present = False
                for previous_comp in intersections:
                    if previous_comp == intersection:
                        already_present = True
                        break

                if not already_present:
                    intersections.append(intersection)

        for intersection in intersections:
            intersection_active_charts = []
            for leaf in tree.depth_first_search():
                chart = leaf[-1].data
                id = chart.id
                intersection_in_chart = I_in_chart(leaf, intersection)
                A = intersection_in_chart.ring()
                if intersection_in_chart == A.ideal(1):
                    continue
                intersection_active_charts.append(id)

            if tree.mode == 'h':
                mult_info = discrepancy_hyperplane(F, intersection)
            else:
                mult_info = discrepancy_double_cover(F, intersection)

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
            origins = "".join([id_to_string(sub_comp.origin) + ", " for sub_comp in comp.sub_components])[:-2]
            ideals = "".join([str(sub_comp.defining_ideal) + ", " for sub_comp in comp.sub_components])[:-2]
            active_charts = "".join([str(sub_comp.active_charts) + ", " for sub_comp in comp.sub_components])[:-2]
            print("Component", i)
            print("Initial chart(s):", origins)
            print("Defining ideal(s)", ideals)
            print("Active charts:", active_charts)
            linear_str = f", Linear = {info['linear']}" if info['linear'] is not None else ""
            print(f"  Codim = {info['codim']}, Mult = {info['mult']}, "
                  f"Discrepancy = {info['discrepancy']}, Admissible = {info['admissible']}{linear_str}")
            print()

    def admissible_indices(self):
        return [i for i, c in enumerate(self.components) if c.mult_info['admissible']]

#Blowups up A at I and returns list of maps A -> A corresponding to affine charts
def blowup_affine_space(chart, I):

    #Recover information about I
    f = chart.f
    id = chart.id
    A = f.parent()
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

                #Restriction A to S
                restriction = sub_i * phi_i * sub_to_elim
                restriction_frac = m2(f'map(frac {A.name()}, frac {S.name()}, matrix {iota.name()})')
                m2(f'use {A.name()}')
                m2(f'use {S.name()}')

                Blf = ((strict_transform(sigma_i, f.ideal(), I)).mingens().entries().flatten())
                charts.append(AffineChart(Blf, sigma_i, id + (i,), blownup_along = gens, inclusion = iota_frac, pi = pi_frac, blowup_ideal = blowup_ideal))
                break

        if len(charts) != i + 1:
            raise ValueError("Blowup ideal of I is not a graph")
    return charts

def translate_and_compare(phi, source_ideal, target_ideal):

    frac_source = source_ideal.ring().fraction_field()
    frac_map = frac_source.coerce_map_from(source_ideal.ring())
    
    gens = source_ideal.gens()
    mapped = [phi(frac_map(g)) for g in gens]

    numerators = [g.numerator() for g in mapped]
    denominators  = [g.denominator() for g in mapped]

    translated = sg.ideal(numerators)
    denom_ideal = sg.ideal(denominators)

    translated = translated.saturation(denom_ideal)[0]

    return translated.radical() == target_ideal.radical()


def same_component(sub_comp : SubComponent, other_comp : SubComponent, tree : ChartTree): 

    origin_A = sub_comp.origin
    origin_B = other_comp.origin

    I = sub_comp.defining_ideal
    J = other_comp.defining_ideal

    if origin_A == origin_B:
        return I.radical() == J.radical()

    phi = tree.transition_maps.get((origin_A, origin_B), False)
    if phi:
        return translate_and_compare(phi, J, I)
    else:
        raise KeyError(f"ChartTree has no transition map from {origin_A} to {origin_B}")
    
    return False

#Takes in a linear center I and returns the map which moves I to ideal(x0,x1,...,xi)
def lcov_of_linear_center(minimal_generators, R):
    n = R.ngens()
    r = len(minimal_generators)
    k = R.base_ring()
    V = sg.VectorSpace(k, n)
    variables = R.gens()

    coeff_rows = []
    b = []
    for f in minimal_generators:
        if f.degree() > 1:
            raise ValueError("All generators must be linear (degree 1) polynomials.")
        row = sg.vector(k, [f.coefficient(xi) for xi in variables])
        b.append(f.constant_coefficient())
        coeff_rows.append(row)
    b = b + [R(0)] * (n - r)

    M = sg.matrix(k, coeff_rows)
    b_vec = sg.vector(R, b)

    rnk = M.rank()
    if rnk != r:
        raise ValueError(
            f"The generators are not linearly independent, expected rank {r}, got rank {rnk}")

    V = sg.VectorSpace(k, n)
    current_space = V.subspace(coeff_rows)
    extension = list(coeff_rows) #Constructs copy

    for i in range(n):
        if len(extension) == n:
            break
        e_i = V.basis()[i]
        candidate = current_space + V.subspace([e_i])
        if candidate.dimension() > current_space.dimension():
            extension.append(e_i)
            current_space = candidate

    #Constructs map as matrices
    T = sg.matrix(k, extension)
    T_inv = T.inverse()

    #Wraps matrices as homomorphisms
    x_vec = sg.vector(R, list(variables))

    foward_images = list(T_inv * (x_vec - b_vec))
    tau = R.hom(foward_images, R)

    inverse_images = list(T * x_vec + b_vec)
    tau_inv = R.hom(inverse_images, R)

    return tau, tau_inv

def lcov_of_quadratic_center(minimal_generators, R):
    
    n = R.ngens()
    vars_R = R.gens()
    if n != 3:
        raise ValueError(f"Ring {R} is not generated by 3 variables")

    r = len(minimal_generators)
    if r != 2:
        raise ValueError(f"Ideal {minimal_generators} not of type (L, Q)")

    k = R.base_ring()
    V = sg.VectorSpace(k, n)
    variables = R.gens()


    L = [f for f in minimal_generators if f.degree() == 1]
    Q = [f for f in minimal_generators if f.degree() == 2]
    
    try:
        L = L[0]
        Q = Q[0]
    except:
        raise ValueError(f"Ideal {minimal_generators} not of type (L, Q)")

    sup = L.monomials()

    for i, v in enumerate(vars_R):
        if v in sup:
            non_zero_var = v
            non_zero_index = i
            break

    param_ring = sg.PolynomialRing(k, 'y', 4)
    param_vars = param_ring.gens()
    param_map_vars = list(param_vars)[:-1]

    coeff_of_non_zero_var = L.coefficient(non_zero_var)
    solve_for_non_zero_index = (coeff_of_non_zero_var * non_zero_var - L) / coeff_of_non_zero_var

    sub_vars = list(param_vars[:-1])
    sub_vars.insert(non_zero_index, param_ring(1))
    sub_map = R.hom(sub_vars, param_ring)
    solve_for_non_zero_index_in_param = sub_map(solve_for_non_zero_index)

    param_map_vars.insert(non_zero_index, solve_for_non_zero_index_in_param)
    param_hom = R.hom(param_map_vars, param_ring)

    tilde_Q = param_hom(Q)
    tilde_Q_homogenized = tilde_Q.homogenize(param_vars[2])
    conic = sg.Conic(tilde_Q_homogenized)

    rational_point_in_param = conic.rational_point()

    if rational_point_in_param[2] == 0:
        raise ValueError("Generated point at infinity :(")

    rational_point_in_param_dehom = list(rational_point_in_param)[:2]

    s, t = rational_point_in_param_dehom

    rational_point = [funt(s, t, 1) for funt in param_map_vars]



    return tau, tau_inv


    
   


#Takes in a chart, and a given number of generators r
#Returns the standard blowup of ideal(x0,...,x(r-1)) in charts
def blowup_affine_space_linear(chart : AffineChart, r : int):
    R = chart.f.parent()
    f = chart.f
    id = chart.id
    vars_R = R.gens()
    blowup_center = R.ideal(vars_R[:r])
    n = len(vars_R)
    if r > n:
        raise ValueError("Number of generators is greater than dimension")
    elif r <= 1:
        raise ValueError("Attempted to blowup a principal ideal")
    k = R.base_ring()

    #Computes map sigma_i : A^n to A^n by the diagram in Hartshorne p.29
    S = sg.PolynomialRing(k, n + r, 'z')
    vars_S = S.gens()
    charts = []

    #Construct blowup ideal
    blowup_ideal_gens = []
    for i in range(r):
        for j in range(i):
            blowup_ideal_gens.append(vars_S[i] * vars_S[j + n] - vars_S[j] * vars_S[i + n])
    blowup_ideal = S.ideal(blowup_ideal_gens)

    for i in range(r):
        vars_S = S.gens()

        #pi : R \to S
        pi = R.hom([vars_S[j] for j in range(n)], S)

        #iota : S \to S (spec Bl_I R \to spec S)
        iota_vars = []
        for j in range(r):
            if j == i:
                iota_vars.append(vars_S[j])
            else:
                iota_vars.append(vars_S[n + j] * vars_S[i])
        
        iota_vars = iota_vars + list(vars_S[r:])
        iota = S.hom(iota_vars, S)
        #print("pi", pi)
        #print("iota",iota)

        #chart_map : S \to R (spec R \to spec Bl_I R)
        chart_vars_z = [1] * (r)
        #print(vars_R)
        chart_vars_z[i] = vars_R[0]
        chart_vars_z += list(vars_R[1: n - r + 1])
        chart_vars_t = list(vars_R[n - r + 1:])
        chart_vars_t.insert(i, 1)
        #print("chart_vars_z", chart_vars_z)
        #print("chart_vars_t", chart_vars_t)
        chart_vars = chart_vars_z + chart_vars_t
        #print("chart_vars", chart_vars)
        chart_map = S.hom(chart_vars, R)

        #sigma_i : R \to R
        sigma_i = chart_map * iota * pi
        #print("sigma_i", sigma_i)

        #inclusion : S \to R (stored as fractional_map)
        #print("Inclusion", chart_map)
        inclusion_vars = (chart_map * iota).im_gens()
        inclusion = S.fraction_field().hom(inclusion_vars, R.fraction_field())

        #restriction : R \to S (fractional_map)
        frac_R = R.fraction_field()
        frac_S = S.fraction_field()

        restriction_vars_z = [vars_S[i]] + list(vars_S[r:n])
        restriction_vars_t = [v for j, v in enumerate(vars_S[n:]) if j != i]
        restriction_vars_t = [v / vars_S[n + i] for v in restriction_vars_t]
        restriction_vars = restriction_vars_z + restriction_vars_t
        restriction = frac_R.hom(restriction_vars, frac_S)

        #sigma_i_inv_vars = [vars_R[i]]
        #sigma_i_inv_vars += list(vars_R[r :])
        #sigma_i_inv_vars += [
                #vars_R[j] / vars_R[i] for j in range(r) if j != i]
        #sigma_i_inv = frac_R.hom(sigma_i_inv_vars, frac_R)
        sigma_i_inv = frac_R.hom(frac_R.gens(), frac_R)

        #print("sigma_i_inv check: ")

        comp = sigma_i_inv * sigma_i.extend_to_fraction_field()
        comp_vars = [comp(g) for g in frac_R.gens()]
        #print(frac_R.hom(comp_vars, frac_R).im_gens())
        comp = sigma_i.extend_to_fraction_field() * sigma_i_inv
        comp_vars = [comp(g) for g in frac_R.gens()]
        #print(frac_R.hom(comp_vars, frac_R).im_gens())

        #print("inclusion check: ")
        comp = inclusion * restriction 
        comp_vars = [comp(g) for g in frac_R.gens()]
        #print(frac_R.hom(comp_vars, frac_R).im_gens())
        comp = restriction * inclusion
        comp_vars = [comp(g) for g in frac_S.gens()]
        #print(frac_S.hom(comp_vars, frac_S).im_gens())

        if chart.tau is not None:
            f_i = strict_transform(sigma_i, chart.tau(R.ideal(f)), blowup_center)
            f_i = m2(f_i)
            l = f_i.mingens().flatten().entries()
            if len(l) == 1:
                f_i = l.first().sage()[0]
            else:
                raise ValueError("Blowup ideal is not principal in all charts")
        else:
            raise ValueError("Linear change of variables has not been computed")

        #print("f_i", f_i)
        #print("sigma_i", sigma_i)
        #print("id", id + (i,))
        #print("blowup_center", blowup_center)
        #print("inclusion", inclusion)
        #print("restriciton", restriction)
        #print("blowup_ideal", blowup_ideal)
        charts.append(AffineChart(f_i, sigma_i, id + (i,), blowup_center, inclusion, restriction, blowup_ideal, sigma_inv = sigma_i_inv))
    return charts
    

#Helper with gb elimination
def elimination_gb(ideal, dep_vars):
    R = ideal.ring()
    ind_vars = [v for v in R.gens() if v not in dep_vars]
    order = TermOrder('degrevlex', len(dep_vars)) + TermOrder("degrevlex", len(ind_vars))
    names = [str(v) for v in dep_vars] + [str(v) for v in ind_vars]
    E = sg.PolynomialRing(R.base_ring(), names, order=order)

    to_E = R.hom([E(str(v)) for v in R.gens()], E)
    from_E = E.hom([R.gens_dict()[str(v)] for v in E.gens()], R)

    J = E.ideal([to_E(g) for g in ideal.gens()])
    return [from_E(g) for g in J.groebner_basis()]

def graph_relations(ideal, dep_vars):
    dep_vars = list(dep_vars)
    relations = {}
    for g in elimination_gb(ideal, dep_vars):
        support = [c for c in dep_vars if g.degree(c) > 0]
        if len(support) != 1:
            return None
        c = support[0]
        if g.degree(c) != 1 or c in relations:
            return None
        coeff = g.coefficient({c : 1})
        if not coeff.is_constant():
            return None
        relations[c] = -(g - coeff * c) / coeff
    if len(relations) != len(dep_vars):
        return None
    return relations


    


#TO DO: Write using Plucker Type ideal
def blowup_affine_space_complete_intersection(chart : AffineChart, mingens, verify = True):
    R = chart.f.parent()
    r = len(mingens)
    f = chart.f
    id = chart.id
    vars_R = R.gens()
    blowup_center = R.ideal(mingens)
    n = len(vars_R)
    if r > n:
        raise ValueError("Number of generators is greater than dimension")
    elif r <= 1:
        raise ValueError("Attempted to blowup a principal ideal")
    k = R.base_ring()

    #Computes map sigma_i : A^n to A^n by the diagram in Hartshorne p.29
    S = sg.PolynomialRing(k, n + r, 'z')
    vars_S = S.gens()
    frac_R = R.fraction_field()
    frac_S = S.fraction_field()
    charts = []

    #pi : R \to S
    pi = R.hom([vars_S[j] for j in range(n)], S)
    mingens_S = [pi(g) for g in mingens]


    for i in range(r):
        #Construct blowup ideal
        blowup_ideal_gens = []
        for j in range(r):
            if i != j:
                blowup_ideal_gens.append(mingens_S[i] * vars_S[n + j] - mingens_S[j])
        blowup_ideal = S.ideal(blowup_ideal_gens)

        is_chart = False
        for dep_vars in combinations(vars_S[:n], r - 1):
            
            #Gets relations from blowup ideal using gb elimination
            relations = graph_relations(blowup_ideal, dep_vars)
            if relations is None:
                continue

            ind_idx = [j for j, v in enumerate(vars_S) if (v not in dep_vars) and (j != n + i)]

            is_graph = True


            phi_i_vars = []
            for j, v in enumerate(vars_S):
                if v in relations:
                    phi_i_vars.append(relations[v])
                elif n + i == j:
                    phi_i_vars.append(S(1))
                else:
                    phi_i_vars.append(v)

            phi_i = S.hom(phi_i_vars, S)

            chart_vars = [R(0)] * (n + r)
            for pos, j in enumerate(ind_idx):
                chart_vars[j] = vars_R[pos]
            chart_map = S.hom(chart_vars, R)

            sigma_i = chart_map * phi_i * pi

            frac_R = R.fraction_field()
            frac_S = S.fraction_field()

            #inclusion : S \to R (stored as fractional_map)
            #print("Inclusion", chart_map)
            inclusion_vars = (chart_map * phi_i).im_gens()
            inclusion = S.fraction_field().hom(inclusion_vars, R.fraction_field())

            frac_R = R.fraction_field()
            frac_S = S.fraction_field()


            #restriction : R \to S (fractional_map)
            restriction_vars = [frac_S(vars_S[j]) if j < n 
                                else frac_S(vars_S[j]) / frac_S(vars_S[n + i]) for j in ind_idx]
            restriction = frac_R.hom(restriction_vars, frac_S)

            sigma_i_inv_vars = []
            for j in ind_idx:
                if j < n:
                    sigma_i_inv_vars.append(frac_R(vars_R[j]))
                else:
                    m = j - n
                    sigma_i_inv_vars.append(frac_R(mingens[m]) / frac_R(mingens[i]))
            sigma_i_inv = frac_R.hom(sigma_i_inv_vars, frac_R)


            if verify:
                for x in vars_R:
                    assert sigma_i_inv(frac_R(sigma_i(x))) == frac_R(x)
                for w in frac_R.gens():
                    assert inclusion(restriction(w)) == w
            #print("sigma_i_inv check: ")

            #comp = sigma_i_inv * sigma_i.extend_to_fraction_field()
            #comp_vars = [comp(g) for g in frac_R.gens()]
            #print(frac_R.hom(comp_vars, frac_R).im_gens())
            #comp = sigma_i.extend_to_fraction_field() * sigma_i_inv
            #comp_vars = [comp(g) for g in frac_R.gens()]
            #print(frac_R.hom(comp_vars, frac_R).im_gens())

            #print("inclusion check: ")
            #comp = inclusion * restriction 
            #comp_vars = [comp(g) for g in frac_R.gens()]
            #print(frac_R.hom(comp_vars, frac_R).im_gens())
            #comp = restriction * inclusion
            #comp_vars = [comp(g) for g in frac_S.gens()]
            #print(frac_S.hom(comp_vars, frac_S).im_gens())

            if chart.tau is not None:
                f_i = strict_transform(sigma_i, chart.tau(R.ideal(f)), blowup_center)
                f_i = m2(f_i)
                l = f_i.mingens().flatten().entries()
                if len(l) == 1:
                    f_i = l.first().sage()[0]
                else:
                    raise ValueError("Blowup ideal is not principal in all charts")
            else:
                raise ValueError("Linear change of variables has not been computed")

        #print("f_i", f_i)
        #print("sigma_i", sigma_i)
        #print("id", id + (i,))
        #print("blowup_center", blowup_center)
        #print("inclusion", inclusion)
        #print("restriciton", restriction)
        #print("blowup_ideal", blowup_ideal)
            charts.append(AffineChart(f_i, sigma_i, id + (i,), blowup_center, inclusion, restriction, blowup_ideal, sigma_inv = sigma_i_inv))
            is_chart = True
            break

        if not is_chart:
            raise ValueError(f"Blowup of {mingens} cannot be expressed in charts")
    return charts

#Computes the affine chart tree corresponding to blowing up component n
def blowup_component(ledger : SingularLedger, i : int, is_intersection_point = False):
    if is_intersection_point:
        singular_component = ledger.initial_intersections[i]
        del ledger.initial_intersections[i]
    else:
        singular_component = ledger.components[i]
        del ledger.components[i]

    singular_component.blownup = True
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
    
            if active_chart.skipped == True:
                continue

            f = active_chart.f
            A = f.parent()
            if f == A(1):
                continue

            pullback_I = I_in_chart(leaf, I)

            if pullback_I == A.ideal(1):
                continue

            if not is_smooth_center(pullback_I):
                raise ValueError(f"Reached a non smooth center {pullback_I}")

            print(ideal_to_string(pullback_I))
            pullback_I_m2 = m2(pullback_I)
            generators = [g.sage() for g in pullback_I_m2.mingens().entries().flatten()]
            r = len(generators)
            
            #Extract generator information
            is_non_linear = False
            #if is_non_principal_linear_center(generators):
            if is_non_principal_linear_center(generators):
                active_chart.tau, active_chart.tau_inv = lcov_of_linear_center(generators, A)
            elif is_non_principal_complete_intersection(generators):
                active_chart.tau, active_chart.tau_inv = A.hom(A.gens(), A), A.hom(A.gens(), A)
                is_non_linear = True
            else:
                raise ValueError(f"Blowup center given by {generators} is non principle, singular, or not a complete intersection")


            if is_non_linear:
                #TO DO: Include exceptional divisor in the case of odd multiplicity
                try:
                    blowup_charts = blowup_affine_space_complete_intersection(active_chart, generators)
                except:
                    active_chart.skipped = True
                    continue
                if tree.mode == 'd' and order_along_ideal(f, pullback_I) % 2:
                    for idx, blowup_chart in enumerate(blowup_charts):
                        sigma_i = blowup_chart.sigma
                        blowup_chart.f *= sigma_i(generators[idx])
            else:
                blowup_charts = blowup_affine_space_linear(active_chart, r)
                if tree.mode == 'd' and order_along_ideal(f, pullback_I) % 2:
                    for blowup_chart in blowup_charts:
                        blowup_chart.f *= A.gen(0)

            for blowup_chart in blowup_charts:
                blowup_chart_node = Node(blowup_chart)
                tree.nodes[blowup_chart.id] = blowup_chart_node
                active_chart_node.children.append(blowup_chart_node)

            for component in ledger.components:
                if component.blownup:
                    continue

                for sub_component in component.sub_components:
                    J, sub_comp_origin, sub_comp_active_charts = (sub_component.defining_ideal, sub_component.origin, sub_component.active_charts)
    
                    if end_id not in sub_comp_active_charts:
                        continue
                
                    sub_component.active_charts.remove(end_id)
                    sub_leaf = tree.find_sub_path_from_ids(sub_comp_origin, end_id)
                    for blowup_chart_node in active_chart_node.children:
                        blowup_chart = blowup_chart_node.data
                        sigma_J = I_in_chart(sub_leaf + [blowup_chart_node], J)
                        B = sigma_J.ring()
                        if not sigma_J == B.ideal(1):
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
                        sigma_J = I_in_chart(sub_leaf + [blowup_chart_node], J)
                        B = sigma_J.ring()
                        if not sigma_J == B.ideal(1):
                            sub_component.active_charts.append(blowup_chart.id)

            print("----------")
            print(f"Blowup of {end_id}:")


            #Finds all new singular sub_componenets
            new_singular_sub_components = []

            for k, blowup_chart_node in enumerate(active_chart_node.children):
                blowup_chart = blowup_chart_node.data
                print(blowup_chart.f)
                Blf = blowup_chart.f
                B = Blf.parent()
                print("Chart", k)
                print("f =", Blf)
                if Blf == B(1):
                    continue

                E = exceptional_divisor(blowup_chart)
                singular_locus = m2(f"decompose ideal singularLocus ideal {m2(Blf).name()}")
                for sub_comp in singular_locus:
                    sub_comp = B.ideal([B(str(g)) for g in sub_comp.sage().gens()])

                    if E <= sub_comp:
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
                existing_comp = False
                if not new_singular_components:
                    chart = tree.nodes[sub_comp.origin].data
                    
                    if tree.mode == 'h':
                        new_singular_components.append(SingularComponent([sub_comp], discrepancy_hyperplane(chart.f, sub_comp.defining_ideal)))
                    else:
                        new_singular_components.append(SingularComponent([sub_comp], discrepancy_double_cover(chart.f, sub_comp.defining_ideal)))
                    continue

                for comp in new_singular_components:
                    old_sub_comp = comp.sub_components[0]
                    if same_component(sub_comp, old_sub_comp, tree):
                        comp.sub_components.append(sub_comp)
                        existing_comp = True
                        break

                if existing_comp:
                    continue

                chart = tree.nodes[sub_comp.origin].data
                if tree.mode == 'h':
                    new_singular_components.append(SingularComponent([sub_comp], discrepancy_hyperplane(chart.f, sub_comp.defining_ideal)))
                else:
                    new_singular_components.append(SingularComponent([sub_comp], discrepancy_double_cover(chart.f, sub_comp.defining_ideal)))
                
            ledger.components.extend(new_singular_components)










                            





                
        

        
        

