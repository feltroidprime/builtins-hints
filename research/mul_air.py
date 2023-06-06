import random
from research.tools import split, split_classic, NamedFelt, polynomial_multiplication_terms, EMULATED_PRIME, NATIVE_PRIME, DEGREE, BASE, N_LIMBS, N_LIMBS_UNREDUCED


# LOOKUP Columns:
# 1 column for selection between diff computation or Zero Poly Assertion
# 1 column for selection between add or sub in case of diff computation (not used in case of Zero Poly Assertion)
# In case of diff computation, N_LIMBS_UNREDUCED columns for which index to add/sub the value, 2 value columns for element product to add/sub
# In case of Zero Poly Assertion, N_LIMBS_UNREDUCED columns to select which diff_limb to check, and the 2 val columns are used for the flag value and the carry value
LOOKUP_N_COLUMNS = 1 + 1 + N_LIMBS_UNREDUCED
# Value Columns:
# 2 columns for the element product to add/sub or the flag value and the carry value. 
# All of these values come from the prover, except for a_limbs and b_limbs that comes from the Cairo memory. 
# N_LIMBS_UNREDUCED columns to store the result of the diff computation
VALUE_COLUMNS = 2 + N_LIMBS_UNREDUCED

AIR_N_COLUMNS = LOOKUP_N_COLUMNS + VALUE_COLUMNS
VAL_INDEX = LOOKUP_N_COLUMNS
# Index after lookup 
RES = LOOKUP_N_COLUMNS + 2


## --------------------Pre-proccessed columns (0/1) ----------------------------- ## -------------------- Other columns ------------------------  ##
## | Diff/ZeroPoly Gate | Add/Sub Gate | d0 | d1 | ... | d_{N_LIMBS_UNREDUCED-1} |||| val0 | val1 | res0 | res1 | ... | res_{N_LIMBS_UNREDUCED-1} |


# Reference : https://github.com/keep-starknet-strange/garaga/blob/7db40334d4b6a4afb3ee9bde90d80458a5512df6/src/bn254/fq.cairo#L231-L410
def mul_hint(a, b):
    def poly_mul(a:list, b:list,n=N_LIMBS) -> list:
        assert len(a) == len(b) == n
        result = [0] * N_LIMBS_UNREDUCED
        for i in range(n):
            for j in range(n):
                result[i+j] += a[i]*b[j]
        return result
    def poly_mul_plus_c(a:list, b:list, c:list, n=N_LIMBS) -> list:
        assert len(a) == len(b) == n
        result = [0] * N_LIMBS_UNREDUCED
        for i in range(n):
            for j in range(n):
                result[i+j] += a[i]*b[j]
        for i in range(n):
            result[i] += c[i]
        return result
    def poly_sub(a:list, b:list, n=N_LIMBS_UNREDUCED) -> list:
        assert len(a) == len(b) == n
        result = [0] * n
        for i in range(n):
            result[i] = a[i] - b[i]
        return result
    def abs_poly(x:list):
        result = [0] * len(x)
        for i in range(len(x)):
            result[i] = abs(x[i])
        return result
    def reduce_zero_poly(x:list):
        x = x.copy()
        carries = [0] * (len(x)-1)
        for i in range(0, len(x)-1):
            carries[i] = x[i] // BASE
            x[i] = x[i] % BASE
            assert x[i] == 0
            x[i+1] += carries[i]
        assert x[-1] == 0
        return x, carries
    
    q = a*b // EMULATED_PRIME
    r = a*b % EMULATED_PRIME
    a_limbs = split_classic(a, 'a')
    b_limbs = split_classic(b, 'b')
    q_limbs = split_classic(q, 'q')
    r_limbs = split_classic(r, 'r')
    p_limbs = split_classic(EMULATED_PRIME, 'P')
    val_limbs = poly_mul(a_limbs, b_limbs)
    q_P_plus_r_limbs = poly_mul_plus_c(q_limbs, p_limbs, r_limbs)
    diff_limbs = poly_sub(q_P_plus_r_limbs, val_limbs)
    _, carries = reduce_zero_poly(diff_limbs)
    carries = abs_poly(carries)
    flags =[0] * (N_LIMBS_UNREDUCED-1)
    for i in range(N_LIMBS_UNREDUCED-1):
        flags[i]=NamedFelt(1 if diff_limbs[i] >= 0 else 0, 'flag'+str(i))
        carries[i] = NamedFelt(carries[i], 'carry'+str(i))
    return flags, carries

def build_empty_row(size:int):
    return [NamedFelt(0, "0")]*size


def build_air(a_limbs, b_limbs, q_limbs, r_limbs, p_limbs, flags, carries):
    ROWS=[]
    n_terms = sum(polynomial_multiplication_terms(N_LIMBS))

    for i in range(n_terms + N_LIMBS):
        ROWS.append(build_empty_row(AIR_N_COLUMNS))
        ROWS[i][0] = NamedFelt(0, "0") # Diff computation
        ROWS[i][1] = NamedFelt(0, "0") # Add gate for q*P + r
    for i in range(n_terms+N_LIMBS, 2*n_terms + N_LIMBS):
        ROWS.append(build_empty_row(AIR_N_COLUMNS))
        ROWS[i][0] = NamedFelt(0, "0") # Diff computation
        ROWS[i][1] = NamedFelt(1, "1") # sub gate for - a*b
    for i in range(2*n_terms + N_LIMBS, 2*n_terms + N_LIMBS  + N_LIMBS_UNREDUCED + 1):
        ROWS.append(build_empty_row(AIR_N_COLUMNS))
        ROWS[i][0] = NamedFelt(1, "1") # Zero poly assertion

    # print("empty rows")
    # for r in ROWS:
    #     print(r)

    row_index=0

    # q*P
    for i in range(N_LIMBS):
        for j in range(N_LIMBS):
            index = i + j
            # Index in D
            d_index = index + 2

            ROWS[row_index][d_index] = NamedFelt(1, "1")
            ROWS[row_index][VAL_INDEX] = q_limbs[i]
            ROWS[row_index][VAL_INDEX + 1] = p_limbs[j]

            # Add term to next row:
            ROWS[row_index+1][RES+index] = ROWS[row_index][RES+index] + ROWS[row_index][VAL_INDEX] * ROWS[row_index][VAL_INDEX + 1]
            # Copy other terms:
            for k in range(N_LIMBS_UNREDUCED):
                if k != index:
                    ROWS[row_index+1][RES+k] = ROWS[row_index][RES+k]
            row_index += 1

    # + r 
    for i in range(N_LIMBS):
        index = i
        # Index in D
        d_index = index + 2

        ROWS[row_index][d_index] = NamedFelt(1, "1")
        ROWS[row_index][VAL_INDEX] = r_limbs[i]
        ROWS[row_index][VAL_INDEX + 1] = NamedFelt(1, "1")

        # Add term to next row:
        ROWS[row_index+1][RES+index] = ROWS[row_index][RES+index] + ROWS[row_index][VAL_INDEX] * ROWS[row_index][VAL_INDEX + 1]
        # Copy other terms:
        for k in range(N_LIMBS_UNREDUCED):
            if k != index:
                ROWS[row_index+1][RES+k] = ROWS[row_index][RES+k]
        row_index += 1

    # - a*b 
    for i in range(N_LIMBS):
        for j in range(N_LIMBS):
            index = i + j
            # Index in D
            d_index = index + 2

            ROWS[row_index][d_index] = NamedFelt(1, "1")
            ROWS[row_index][VAL_INDEX] = a_limbs[i]
            ROWS[row_index][VAL_INDEX + 1] = b_limbs[j]

            # Add term to next row:
            ROWS[row_index+1][RES+index] = ROWS[row_index][RES+index] - ROWS[row_index][VAL_INDEX] * ROWS[row_index][VAL_INDEX + 1]
            # Copy other terms:
            for k in range(N_LIMBS_UNREDUCED):
                if k != index:
                    ROWS[row_index+1][RES+k] = ROWS[row_index][RES+k]
            row_index += 1

    # Mocked 0 positive carry for the first limb
    ROWS[row_index][VAL_INDEX] = NamedFelt(1, "1") # Positive flag
    ROWS[row_index][VAL_INDEX + 1] = NamedFelt(0, "0") # Zero Carry 
    ROWS[row_index][2] = NamedFelt(1, "1")

    # Copy previous terms :
    for k in range(N_LIMBS_UNREDUCED):
        ROWS[row_index+1][RES+k] = ROWS[row_index][RES+k]
    row_index += 1

    for i in range(N_LIMBS_UNREDUCED-1):
        ROWS[row_index][VAL_INDEX] = flags[i]
        ROWS[row_index][VAL_INDEX + 1] = carries[i]
        d_index = i + 3
        ROWS[row_index][d_index] = NamedFelt(1, "1")
        # copy previous terms to next row:
        # if i != N_LIMBS_UNREDUCED-2: 
        for k in range(N_LIMBS_UNREDUCED):
            ROWS[row_index+1][RES+k] = ROWS[row_index][RES+k]
        row_index += 1

    return ROWS


def assert_row(rows:list, index:int):
    diff_or_zero_poly_gate = rows[index][0].value
    add_sub_gate = rows[index][1].value
    val0 = rows[index][VAL_INDEX].value
    val1 = rows[index][VAL_INDEX+1].value
    next_val0 = rows[index+1][VAL_INDEX].value
    next_val1 = rows[index+1][VAL_INDEX+1].value
    di = [rows[index][i].value for i in range(2, LOOKUP_N_COLUMNS)]
    res = [rows[index][i].value for i in range(RES, AIR_N_COLUMNS)]
    next_res = [rows[index+1][i].value for i in range(RES, AIR_N_COLUMNS)]
    add_diff_assertion = 0
    for i in range(N_LIMBS_UNREDUCED):
        add_diff_assertion += (res[i] + di[i] * val0 * val1 - next_res[i]) % NATIVE_PRIME
    sub_diff_assertion = 0
    for i in range(N_LIMBS_UNREDUCED):
        sub_diff_assertion += (res[i] - di[i] * val0 * val1 - next_res[i]) % NATIVE_PRIME

    add_diff_assertion = add_diff_assertion % NATIVE_PRIME
    sub_diff_assertion = sub_diff_assertion % NATIVE_PRIME

    # https://github.com/keep-starknet-strange/garaga/blob/7db40334d4b6a4afb3ee9bde90d80458a5512df6/src/bn254/fq.cairo#L374-L406
    zero_poly_assertion_postive_carry = 0
    for i in range(N_LIMBS_UNREDUCED-1):
        zero_poly_assertion_postive_carry += (di[i] * (val0*(res[i] + val1 - next_val1*BASE) + (1-val0)*(res[i] - val1 - next_val1*BASE))) % NATIVE_PRIME

    zero_poly_assertion_negative_carry = 0
    for i in range(N_LIMBS_UNREDUCED-1):
        zero_poly_assertion_negative_carry += (di[i] * (val0*(res[i] + val1 + next_val1*BASE) + (1-val0)*(res[i] - val1 + next_val1*BASE))) % NATIVE_PRIME

    diff_assertion = ((1-add_sub_gate) * add_diff_assertion + add_sub_gate * sub_diff_assertion) % NATIVE_PRIME
    zero_poly_assertion = (next_val0 * zero_poly_assertion_postive_carry + (1-next_val0) * zero_poly_assertion_negative_carry) % NATIVE_PRIME

    copy_constraint = 0
    for i in range(N_LIMBS_UNREDUCED):
        copy_constraint += (res[i] - next_res[i]) % NATIVE_PRIME

    final_assertion = ((1-diff_or_zero_poly_gate) * diff_assertion + diff_or_zero_poly_gate * (zero_poly_assertion + copy_constraint)) % NATIVE_PRIME
    assert final_assertion == 0

    
a=random.randint(0, EMULATED_PRIME-1)
b=random.randint(0, EMULATED_PRIME-1)
q = a*b // EMULATED_PRIME
r = a*b % EMULATED_PRIME
a_limbs = split(a, 'a')
b_limbs = split(b, 'b')
q_limbs = split(q, 'q')
r_limbs = split(r, 'r')
p_limbs = split(EMULATED_PRIME, 'P')
air=build_air(a_limbs, b_limbs, q_limbs, r_limbs, p_limbs, *mul_hint(a, b))


for i, row in enumerate(air):
    print(f"{i}:", row)

print(f"\nMul AIR Size NxM = {len(air)}x{len(air[0])} = {len(air)*len(air[0])} cells")

print(f"Estimated # step cost : (#cells) / 51 = {len(air)*len(air[0]) / 51} steps")

print(f"\n==> Asserting diff computation and Zero polynomial constraints for row 0 to N-1...")
for i in range(0, len(air)-1):
    # print(f"Asserting row {i}...", end='')
    assert_row(air, i)
print("==> Done!")
