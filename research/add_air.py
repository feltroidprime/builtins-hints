import random
from .tools import split, NamedFelt, polynomial_multiplication_terms, EMULATED_PRIME, NATIVE_PRIME, DEGREE, BASE, N_LIMBS, N_LIMBS_UNREDUCED
random.seed(42)


a=random.randint(0, EMULATED_PRIME-1)
b=random.randint(0, EMULATED_PRIME-1)
c = (a+b)%EMULATED_PRIME

a_limbs = split(a, 'a')
b_limbs = split(b, 'b')
c_limbs = split(c, 'c')
p_limbs = split(EMULATED_PRIME, 'P')