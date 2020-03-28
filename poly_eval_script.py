import collections
import random
import argparse
import math

class Node:
    def coeff_compatible(self, rhs):
        return self.sub_poly.coeff_compatible(rhs.sub_poly)

    @property
    def offset(self):
        return self.sub_poly.offset
    @property
    def coeffs(self):
        return self.sub_poly.coeffs



class Cst(Node):
    def __init__(self, coeff_index):
        super().__init__()
        self.sub_poly = SubPoly(coeff_index, 1 << coeff_index)
        self.coeff_index = coeff_index

    def op_count(self, processed_set=None):
        return 0
    @property
    def level(self):
        return 0

    def __str__(self):
        return "a_%d" % self.coeff_index

class Var(Node):
    arity = 0
    def __str__(self):
        return "x"

    def op_count(self, processed_set=None):
        return 0
    @property
    def level(self):
        return 0

class Op(Node):
    arity = None
    def __init__(self, sub_poly, *args):
        super().__init__()
        assert len(args) == self.arity
        self.sub_poly = sub_poly
        self.args = args

    def op_count(self, processed_set=None):
        if self in processed_set:
            # already counted
            return 0
        processed_set.add(self)
        return sum(op.op_count(processed_set) for op in self.args) + 1
    @property
    def level(self):
        return max(op.level for op in self.args) + 1


class FADD(Op):
    arity = 2
    def __str__(self):
        return "FADD({} + {})".format(self.args[0], self.args[1])

class FMUL(Op):
    arity = 2
    def __str__(self):
        return "FMUL({} * {})".format(self.args[0], self.args[1])

class FMA(Op):
    arity = 3
    def __str__(self):
        return "FMA({} + {} * {})".format(self.args[0], self.args[1], self.args[2])

class FDMA(Op):
    arity = 4
    def __str__(self):
        return "FDMA({} * {} + {} * {})".format(self.args[0], self.args[1], self.args[2], self.args[3])

class FDMDA(Op):
    arity = 5
    def __str__(self):
        return "FDMDA({} + {} * {} + {} * {})".format(self.args[0], self.args[1], self.args[2], self.args[3], self.args[4])




class SubPoly:
    def __init__(self, offset=0, coeffs=0):
        # offset: offset with respect to true final polynomial degree
        # coeffs: mask of coefficient present in this sub-polynomial
        # for example a_0 + a_2 x^2 is (offset=0, coeffs = 0x5)
        #             a_2 + a_3 x is (offset=2, coeffs = 0xc)
        self.offset = offset
        self.coeffs = coeffs

    def compatible(sub_poly_0, sub_poly_1):
        return (sub_poly_0.offset == sub_poly_1.offset) and sub_poly_0.compatible(sub_poly_1)

    def coeff_compatible(sub_poly_0, sub_poly_1):
        return (sub_poly_0.coeffs ^ sub_poly_1.coeffs == sub_poly_0.coeffs | sub_poly_1.coeffs)

    @property
    def key(self):
        return self.offset, self.coeffs


def random_compatible_mask(superset_mask):
    """ generate a random mask """
    pass

def generate_estrin_scheme(degree):
    terms = [Cst(index) for index in range(degree+1)]
    var = Var()
    while len(terms) > 1:
        new_terms = []
        while len(terms) >= 2:
            op0 = terms.pop(0)
            op1 = terms.pop(0)
            new_node = FMA(None, op0, op1, var)
            new_terms.append(new_node)
        if len(terms):
            # one last terme remaining
            new_terms.append(terms.pop(0))
        var = FMUL(None, var, var)
        terms = new_terms
    return terms[0]

def generate_horner_scheme(degree):
    terms = [Cst(index) for index in range(degree+1)]
    var = Var()
    op = terms.pop(-1)
    for term in terms[::-1]:
        op = FMA(None, term, var, op)
    return op

def generate_procedural_scheme(coeff_mask=None, offset=0):
    if coeff_mask is None:
        coeff_mask = 2**(degree+1) - 1

    min_index = min([index for index in range(degree + 1) if coeff_mask & (2**index)])
    mask_no_min = coeff_mask ^ (1 << min_index)
    even_mask = sum([2**index for index in range(0, degree +1 ,2) if coeff_mask & (2**index)])
    odd_mask = sum([2**index for index in range(1, degree +1 ,2) if coeff_mask & (2**index)])

    # NOTES/TODO: to be finished
    raise NotImplementedError



def generate_op_map(degree, coeff_mask=None, NUM_RANDOM_SAMPLE=100, operators=None):
    if coeff_mask is None:
        coeff_mask = 2**(degree+1) - 1

    # list of nodes
    op_map = []
    # map node -> level
    level_map = {}
    # map coeffs subset -> node set
    mask_map = collections.defaultdict(set)

    # initializing operation map with polynomial coefficients
    for index in range(degree+1):
        if (coeff_mask >> index) % 2 == 0:
            continue
        op = Cst(index)
        level_map[op] = 0
        op_map.append(op)

    # populating power map
    power_map = {}
    power_level = collections.defaultdict(lambda : None)
    power_level[0] = 0
    power_level[1] = 0
    power_level[2] = 1
    power_map = {0: None, 1: Var(), 2: FMUL(None, Var(), Var())}
    for index in range(degree+1):
        for m in range(1, index):
            lhs = m
            rhs = index - m
            level = max(power_level[lhs], power_level[rhs]) + 1
            if power_level[index] is None or power_level[index] > level:
                power_level[index] = level
                power_map[index] = FMUL(None, power_map[lhs], power_map[rhs])


    # Adding FMA like operation
    def generate_random_fma():
        lhs = random.choice(op_map)
        rhs = random.choice(op_map)
        if lhs.coeff_compatible(rhs):
            delta = lhs.offset - rhs.offset
            sub_poly = SubPoly(min(lhs.offset, rhs.offset), lhs.coeffs | rhs.coeffs)
            if delta == 0:
                node_level = max([level_map[lhs], level_map[rhs]]) + 1
                new_node = FADD(sub_poly, lhs, rhs)
            elif delta > 0:
                node_level = max([level_map[lhs], level_map[rhs], power_level[delta]]) + 1
                new_node = FMA(sub_poly, rhs, lhs, power_map[delta])
            else:
                node_level = max([level_map[lhs], level_map[rhs], power_level[-delta]]) + 1
                new_node = FMA(sub_poly, lhs, rhs, power_map[-delta])
            op_map.append(new_node)
            # computing level
            level_map[new_node] = node_level
            mask_map[new_node.coeffs].add(new_node)

    # Adding FDMA like operation
    def generate_random_fdma():
        lhs = random.choice(op_map)
        rhs = random.choice(op_map)
        if lhs.coeff_compatible(rhs):
            min_offset = min(lhs.offset, rhs.offset)
            if min_offset <= 0:
                # if min offset is zero or negative, no FDMDA is required
                # (a FMA or Add will work
                # and would have been catched (randomly) during previous stage
                return
            for offset in range(0, min_offset):
                # we compute sub_poly so that is offset is offset
                lhs_power = lhs.offset - offset
                rhs_power = rhs.offset - offset
                sub_poly = SubPoly(offset, lhs.coeffs | rhs.coeffs)
                new_node = FDMA(
                    sub_poly,
                    power_map[lhs_power],
                    lhs,
                    power_map[rhs_power],
                    rhs
                )
                node_level = max(power_level[lhs_power], level_map[lhs], power_level[rhs_power], level_map[rhs]) + 1 
                op_map.append(new_node)
                level_map[new_node] = node_level
                mask_map[new_node.coeffs].add(new_node)


    # Adding FDMDA like operation
    def generate_random_fdmda():
        op0 = random.choice(op_map)
        op1 = random.choice(op_map)
        op2 = random.choice(op_map)
        if not(op0.coeff_compatible(op1)) or not(op1.coeff_compatible(op2)) or not(op0.coeff_compatible(op2)):
            return
        min_offset = min(op0.offset, op1.offset, op2.offset)
        sub_poly = SubPoly(min_offset, op0.coeffs | op1.coeffs | op2.coeffs)
        op0, op1, op2 = sorted([op0, op1, op2], key=lambda op: op.offset)
        op1_power = op1.offset - min_offset
        op2_power = op2.offset - min_offset
        if op1_power == 0 or op2_power == 0:
            # if any power is 0, this is not really a FDMDA
            return
        new_node = FDMDA(
            sub_poly,
            op0,
            op1,
            power_map[op1_power],
            op2,
            power_map[op2_power]
        )
        node_level = max([level_map[op0], level_map[op1], level_map[op2], power_level[op1_power], power_level[op2_power]]) + 1
        op_map.append(new_node)
        level_map[new_node] = node_level
        mask_map[new_node.coeffs].add(new_node)

    for _ in range(NUM_RANDOM_SAMPLE):
        if "fma" in operators:
            generate_random_fma()
        if "fdma" in operators:
            generate_random_fdma()
        if "fdmda" in operators:
            generate_random_fdmda()

    if "fdmda" in operators:
        for _ in range(NUM_RANDOM_SAMPLE):
            generate_random_fdmda()

    print("len(op_map) is {} ({:.2f}% sample(s))".format(len(op_map), len(op_map) / NUM_RANDOM_SAMPLE * 100))

    # looking for complete candidate
    candidates = [node for node in op_map if node.coeffs == coeff_mask]
    also_candidates = list(mask_map[coeff_mask])
    print("# of candidates {}/{}".format(len(candidates), len(also_candidates)))
    min_depth_scheme = min(candidates, key=lambda node: (level_map[node], node.op_count(set())))
    min_num_op_scheme = min(candidates, key=lambda node: (node.op_count(set()), level_map[node]))
    print("min_depth_scheme is {} with level {} and {} op(s)".format(str(min_depth_scheme), level_map[min_depth_scheme], min_depth_scheme.op_count(set())))
    print("min_num_op_scheme is {} with level {} and {} op(s)".format(str(min_num_op_scheme), level_map[min_num_op_scheme], min_num_op_scheme.op_count(set())))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Polynomial scheme explorator')
    parser.add_argument('--num-steps', action='store',
                        default=100000, type=int,
                        help='number of random setps')
    parser.add_argument('--degree', action='store',
                        default=7, type=int,
                        help='polynomial degree (exclusive option with mask)')
    parser.add_argument('--coeff-mask', action='store',
                        default=None, type=lambda v: int(v, base=2),
                        help='coefficient mask (exclusive option with degree)')
    parser.add_argument('--operators', action='store',
                        default=["fma", "fadd", "fmul", "fdma", "fdmda"], type=lambda v: v.split(","),
                        help="',' separated list of allowed operators")
    args = parser.parse_args()
    if not args.coeff_mask is None:
        coeff_mask = args.coeff_mask
        degree = int(math.floor(math.log2(coeff_mask)))
        print("degree: {}".format(degree))
    else:
        coeff_mask = None
        degree = args.degree
    # standard scheme
    horner_scheme = generate_horner_scheme(degree)
    print("horner_scheme is {}, level={}, #ops={}".format(str(horner_scheme), horner_scheme.level, horner_scheme.op_count(set())))
    estrin_scheme = generate_estrin_scheme(degree)
    print("estrin_scheme is {}, level={}, #ops={}".format(str(estrin_scheme), estrin_scheme.level, estrin_scheme.op_count(set())))

    # random exploration
    generate_op_map(degree, coeff_mask=coeff_mask, NUM_RANDOM_SAMPLE=args.num_steps, operators=args.operators)

