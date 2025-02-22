from py_ecc.bn128 import G1, G2, multiply, add, pairing, neg, curve_order, FQ, FQ2
from py_ecc.bn128.bn128_pairing import final_exponentiate
import random

def evaluate(poly, x):
    result = FQ(0)
    power = FQ(1)
    for coeff in poly:
        result += coeff * power
        power *= x
    return result

def div(numerator, denominator):
    quotient = []
    remainder = numerator[:]
    while len(remainder) >= len(denominator):
        lead_coeff = remainder[0] / denominator[0]
        quotient.append(lead_coeff)
        remainder = [a - b * lead_coeff for a, b in zip(remainder, denominator + [FQ(0)] * (len(remainder) - len(denominator)))]
        remainder.pop(0)
    return quotient

class KZG:
    def __init__(self, g1=G1, g2=G2, degree=2):
        self.g1 = g1
        self.g2 = g2
        self.degree = degree
        self.g2_tau = (0, 0)
        self.crs_g1 = []
        self.crs_g2 = []
    
    def setup(self, secret):
        self.crs_g1 = [multiply(self.g1, pow(secret, i, curve_order)) for i in range(self.degree + 1)]
        self.crs_g2 = [multiply(self.g2, pow(secret, i, curve_order)) for i in range(self.degree + 1)]
        self.g2_tau = multiply(self.g2, secret)
    
    def commit(self, poly):
        commitment = (0, 0)
        for i in range(self.degree + 1):
            commitment = add(commitment, multiply(self.crs_g1[i], poly[i]))
        return commitment
    
    def open(self, poly, point):
        value = evaluate(poly, point)
        numerator = [poly[0] - value] + poly[1:]
        denominator = [-point, FQ(1)]
        quotient = div(numerator, denominator)
        pi = (0, 0)
        for i in range(len(quotient)):
            pi = add(pi, multiply(self.crs_g1[i], quotient[i]))
        return pi
    
    def verify(self, point, value, commitment, pi):
        lhs = pairing(pi, add(self.g2_tau, multiply(self.g2, -point)))
        rhs = pairing(add(commitment, multiply(self.g1, -value)), self.g2)
        return final_exponentiate(lhs) == final_exponentiate(rhs)

# Implementation

def main():
# Setup Phase
    kzg = KZG(degree=3)  # Create KZG object with polynomial degree 3
    secret = random.randint(1, curve_order - 1)  # Trusted setup secret
    kzg.setup(secret)  # Generate CRS (Common Reference String)

    # Commitment Phase
    poly = [FQ(3), FQ(5), FQ(2)]  # 3 + 5x + 2x^2
    commitment = kzg.commit(poly)
    point = FQ(4)
    pi = kzg.open(poly, point)
    value = evaluate(poly, point)
    print("Verification:", kzg.verify(point, value, commitment, pi))

if __name__ == "__main__":
    main()