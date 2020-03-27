
class Node:
    @property
    def offset(self):
        return self.sub_poly.offset
    @property
    def coeffs(self):
        return self.sub_poly.coeffs


class Cst(Node):
    def __init__(coeff_index):
        self.sub_poly = SubPoly(0, 1 << coeff_index)
        self.coeff_index = index

class Op(Node):
    arity = None
    def __init__(self, sub_poly, *args):
        assert len(args) == self.arity
        self.sub_poly = sub_poly
        self.args = args


class FADD(Op):
    arity = 2

class FMUL(Op):
    arity = 2

class FMA(Op):
    arity = 3

class FDMA(Op):
    arity = 4

class FDMDA(Op):
    arity = 5


power_map = {}
op_map


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


def generate_op_map(degree, coeff_mask=None, NUM_RAMDOM_SAMPLE=100):
    if coeff_mask is None:
        coeff_mask = 2**(degree+1) - 1
    # initializing operation map with polynomial coefficients
    op_map = set()
    level_map = {}
    for sub_poly in [[SubPoly(offset=index, coeffs=2**index) for index in range(degree+1) if (coeff_mask >> index) % 2 != 0]]:
        op = Cst(index)
        level_map[op] = 0
        op_map.add(Cst)

    # populating power map
    power_map = {}
    power_level = collection.defaultdict(lambda : None)
    power_level[0] = 0
    power_level[1] = 0
    power_level[2] = 1
    for index in range(degree+1):
        for m in range(index):
            if (index // m) * m == index:
                lhs = m
                rhs = index // m
                level = max(power_level[lhs], power_level[rhs]) + 1
                if power_level[index] is None or power_level[index] > level:
                    power_level[index] = level
                    power_map[index] = Mul(lhs, rhs)

    # Adding FMA like operation
    for _ in range(NUM_RANDOM_SAMPLE):
        lhs = random.choice(op_map)
        rhs = random.choice(op_map)
        if lhs.coeff_compatible(rhs):
            delta = lhs.offset - rhs.offset
            sub_poly = SubPoly(min(lhs.offset, rhs.offset), lhs.coeffs | rhs.coeffs)
            if delta == 0:
                new_node = FADD(sub_poly, lhs, rhs)
            elif delta > 0
                new_node = FMA(sub_poly, rhs, lhs, delta)
            else:
                new_node = FMA(sub_poly, hs, rhs, -delta)
            op_map.add(new_node)
            # computing level
            node_level = max(level_map[lhs], level_map[rhs]) + 1
            level_map[new_node] = node_level

    # Adding FDMA like operation
    for _ in range(NUM_RANDOM_SAMPLE):
        lhs = random.choice(op_map)
        rhs = random.choice(op_map)
        if lhs.coeff_compatible(rhs):
            min_offset = min(lhs.offset, rhs.offset)
            if min_offset <= 0:
                # if min offset is zero or negative, no FDMDA is required
                # (a FMA or Add will work
                # and would have been catched (randomly) during previous stage
                continue
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
                op_map.add(new_node)
                level_map[new_node] = node_level

    # Adding FDMDA like operation
    for _ in range(NUM_RANDOM_SAMPLE):
        op0 = random.choice(op_map)
        op1 = random.choice(op_map)
        op2 = random.choice(op_map)
        if not(op0.coeff_compatible(op1)) or not(op1.coeff_compatible(op2)) or not(op0.coeff_compatible(op2)):
            continue
        min_offset = min(op0.offset, op1.offset, op2.offset)
        sub_poly = SubPoly(min_offset, op0.coeffs | op1.coeffs | op2.coeffs)
        op0, op1, op2 = sorted([op0, op1, op2], key=lambda op: op.offset)
        op1_power = op1.offset - min_offset
        op2_power = op2.offset - min_offset
        new_node = FDMDA(
            sub_poly
            op0,
            op1,
            power_map[op1_power],
            op2,
            power_map[op2_power]
        )
        node_level = max([level_map[op0], level_map[op1], level_map[op2], power_level[op1_power], power_level[op2_power]) + 1
        op_map.add(new_node)
        level_map[new_node] = node_level

    # looking for complete candidate
    candidates = [node for node in op_map if node.coeffs == coeff_mask]
    best_candidate = min(candidates, key=lambda node: level_map[node])
    print("best_candidate is {} with level {}".format(best_candidate, level_map[best_candidate]))

generate_op_map(3)



