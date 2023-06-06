
BASE=2**86
DEGREE = 2
N_LIMBS = 3
N_LIMBS_UNREDUCED = 2*N_LIMBS - 1
EMULATED_PRIME = 21888242871839275222246405745257275088548364400416034343698204186575808495617 # bn254 
NATIVE_PRIME = 3618502788666131213697322783095070105623107215331596699973092056135872020481 # stark prime


class NamedFelt:
    def __init__(self, value, name):
        self.value = value
        self.name = name

    def __add__(self, other):
        if isinstance(other, NamedFelt):
            if other.value == 0:
                return NamedFelt(self.value, self.name)
            elif self.value == 0:
                return NamedFelt(other.value, other.name)
            else:
                return NamedFelt((self.value + other.value) % NATIVE_PRIME, f"{self.name}+{other.name}")
        else:
            raise TypeError("Unsupported operand type.")
    
    def __iadd__(self, other):
        if isinstance(other, NamedFelt):
            if other.value == 0:
                pass
            elif self.value == 0:
                self.value, self.name = other.value, other.name
            else:
                self.value = (self.value + other.value) % NATIVE_PRIME
                self.name = f"{self.name}+{other.name}"
            return self
        else:
            raise TypeError("Unsupported operand type.")

    def __sub__(self, other):
        if isinstance(other, NamedFelt):
            return NamedFelt((self.value - other.value) % NATIVE_PRIME, f"{self.name}-{other.name}")
        else:
            raise TypeError("Unsupported operand type.")

    def __mul__(self, other):
        if isinstance(other, NamedFelt):
            if other.value == 1:
                return NamedFelt(self.value, self.name)
            elif self.value == 1:
                return NamedFelt(other.value, other.name)
            else:
                return NamedFelt(self.value * other.value % NATIVE_PRIME, f"{self.name}*{other.name}")
        else:
            raise TypeError("Unsupported operand type.")

    def __repr__(self):
        return f"{self.name}"

    
def polynomial_multiplication_terms(n_limbs):
    result = []
    for i in range(2 * n_limbs - 1):
        if i < n_limbs:
            result.append(i + 1)
        else:
            result.append(2 * n_limbs - i - 1)
    return result

def split(x:int, name:str, degree=DEGREE, base=BASE):
    coeffs = []
    for n in range(degree, 0, -1):
        q, r = divmod(x, base ** n)
        coeffs.append(q)
        x = r
    coeffs.append(x)
    coeffs= coeffs[::-1]
    for i in range(len(coeffs)):
        coeffs[i]=NamedFelt(coeffs[i], f'{name}{i}')
    return coeffs

def split_classic(x:int, name:str, degree=DEGREE, base=BASE):
    coeffs = []
    for n in range(degree, 0, -1):
        q, r = divmod(x, base ** n)
        coeffs.append(q)
        x = r
    coeffs.append(x)
    coeffs= coeffs[::-1]
    return coeffs