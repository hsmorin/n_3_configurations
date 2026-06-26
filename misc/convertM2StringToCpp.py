import re
from sys import argv

def convert_powers_to_powermod(s: str, p: str = "p") -> str:
    pattern = r'\b(x\d+)\^(\d+)\b'
    replacement = lambda m: f"powerMod({m.group(1)}, {m.group(2)}, {p})"
    return re.sub(pattern, replacement, s)

if __name__ == "__main__":
    st = "x1^4*x2^2*x3^2-2*x1^2*x2^4*x3^2+x2^6*x3^2+4*x1^2*x2^3*x3^3-4*x2^5*x3^3-2*x1^2*x2^2*x3^4+6*x2^4*x3^4-4*x2^3*x3^5+x2^2*x3^6-4*x1^4*x2^2*x3*x4+6*x1^3*x2^3*x3*x4-2*x1^2*x2^4*x3*x4+2*x1*x2^5*x3*x4-2*x2^6*x3*x4-6*x1^3*x2^2*x3^2*x4-2*x1*x2^4*x3^2*x4+8*x2^5*x3^2*x4+2*x1^2*x2^2*x3^3*x4-2*x1*x2^3*x3^3*x4-12*x2^4*x3^3*x4+2*x1*x2^2*x3^4*x4+8*x2^3*x3^4*x4-2*x2^2*x3^5*x4+4*x1^4*x2^2*x4^2-8*x1^3*x2^3*x4^2+5*x1^2*x2^4*x4^2-2*x1*x2^5*x4^2+x2^6*x4^2+2*x1^4*x2*x3*x4^2+4*x1^3*x2^2*x3*x4^2-2*x1*x2^4*x3*x4^2-4*x2^5*x3*x4^2+4*x1^3*x2*x3^2*x4^2-7*x1^2*x2^2*x3^2*x4^2+10*x1*x2^3*x3^2*x4^2+6*x2^4*x3^2*x4^2+2*x1^2*x2*x3^3*x4^2-6*x1*x2^2*x3^3*x4^2-4*x2^3*x3^3*x4^2+x2^2*x3^4*x4^2-4*x1^4*x2*x4^3+6*x1^3*x2^2*x4^3-6*x1^2*x2^3*x4^3+4*x1*x2^4*x4^3-6*x1^3*x2*x3*x4^3+8*x1^2*x2^2*x3*x4^3-8*x1*x2^3*x3*x4^3-2*x1^2*x2*x3^2*x4^3+4*x1*x2^2*x3^2*x4^3+x1^4*x4^4"
    print(convert_powers_to_powermod(st))
