import collections
import random

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

    @property
    def op_count(self):
        return 0

    def __str__(self):
        return "a_%d" % self.coeff_index

class Var(Node):
    arity = 0
    def __str__(self):
        return "x"

    @property
    def op_count(self):
        return 0

class Op(Node):
    arity = None
    def __init__(self, sub_poly, *args):
        super().__init__()
        assert len(args) == self.arity
        self.sub_poly = sub_poly
        self.args = args

    @property
    def op_count(self):
        return sum(op.op_count for op in self.args) + 1


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


def generate_op_map(degree, coeff_mask=None, NUM_RANDOM_SAMPLE=100):
    if coeff_mask is None:
        coeff_mask = 2**(degree+1) - 1
    # initializing operation map with polynomial coefficients
    op_map = []
    level_map = {}
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
                node_level = max(power_level[lhs_power], level_map[rhs], power_level[rhs_power], level_map[rhs]) + 1 
                op_map.append(new_node)
                level_map[new_node] = node_level


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

    for _ in range(NUM_RANDOM_SAMPLE):
        generate_random_fma()
        generate_random_fdma()
        generate_random_fdmda()

    print("len(op_map) is {}".format(len(op_map)))

    # looking for complete candidate
    candidates = [node for node in op_map if node.coeffs == coeff_mask]
    best_candidate = min(candidates, key=lambda node: (level_map[node], node.op_count))
    print("best_candidate is {} with level {} and {} op(s)".format(str(best_candidate), level_map[best_candidate], best_candidate.op_count))

generate_op_map(7, NUM_RANDOM_SAMPLE=10000000)



