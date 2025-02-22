#
# This code is not audited. It is for research and education purposes only
# DO NOT USE IT IN PRODUCTION ENVIRONMENT
#

import galois
from py_ecc.bls12_381 import G1, G2, G12, add, multiply, curve_order, pairing, eq, neg, Z1
import random
from functools import reduce

class KZGCommitment:
    def __init__(self, curve_order, trusted_setup_degree = 10):
        print("Initializing Galois Field with BLS12_318 curve order. This may take a while.")
        self.GF = galois.GF(curve_order)
        print("Galois Field Initialized")
        (self.trusted_setup_g1, self.trusted_setup_g2) = self.generate_trusted_setup(trusted_setup_degree) # change this to the highest degree of the polynomial you'll use.

    # convert an integer to a finite field value
    def convert_int_to_gf(self, num):
        if num < 0:
            return -self.GF(-num)
        else:
            return self.GF(num)

    # interpolate the given points to return a galois polynomial
    def lagrange_interpolation_galois(self, xs, ys):
        assert len(xs) == len(ys), "Length mismatch"
        length = len(xs)

        def pi(j, y):
            p_result = self.GF(1)
            for k in range(length):
                if k == j:
                    continue
                q, _ = divmod(galois.Poly([self.GF(1), -xs[k]], self.GF), (xs[j] - xs[k]))
                p_result = p_result * q
            return y * p_result

        result = self.GF(0)
        for i in range(length):
            result = result + pi(i, ys[i])

        return result

    # lagrange interpolation to convert vector to a galois polynomial
    def vector_to_polynomial(self, vector):
        xs = [self.convert_int_to_gf(i) for i in range(len(vector))]
        vector_gf = [self.convert_int_to_gf(x) for x in vector]
        return self.lagrange_interpolation_galois(xs, vector_gf)

    # generate [s^i * G1] and [s^i * G2]
    def generate_trusted_setup(self, degree):
        print("Generating trusted setup")
        s = self.GF(random.randint(1, curve_order))
        return ([multiply(G1, int(s ** i)) for i in range(degree + 1)], [multiply(G2, int(s ** i)) for i in range(degree + 1)])

    def evaluate_at_trusted_setup(self, polynomial, trusted_setup):
        trusted_setup = trusted_setup[:polynomial.degree + 1]
        return reduce(add, (multiply(point, int(coeff)) for point, coeff in zip(trusted_setup, polynomial.coeffs[::-1])), Z1)

    def commit_polynomial(self, polynomial):
        return self.evaluate_at_trusted_setup(polynomial, self.trusted_setup_g1)

    # generate a proof for a given set of point
    def generate_proof(self, polynomial, points):
        x_s = [self.convert_int_to_gf(x) for x, y in points]
        y_s = [self.convert_int_to_gf(y) for x, y in points]

        points_polynomial = self.lagrange_interpolation_galois(x_s, y_s)
        numerator = polynomial - points_polynomial
        denominator = self.GF(1)
        for x in x_s:
            denominator = denominator * galois.Poly([self.GF(1), -x], self.GF)

        quotient, reminder = divmod(numerator, denominator)
        assert reminder == 0, "This is not a valid proof"

        return self.evaluate_at_trusted_setup(quotient, self.trusted_setup_g1)

    # Verify proof
    def verify_proof(self, commitment, points, proof):
        x_s = [self.convert_int_to_gf(x) for x, y in points]
        y_s = [self.convert_int_to_gf(y) for x, y in points]
        points_polynomial = self.lagrange_interpolation_galois(x_s, y_s)

        z = galois.Poly([1], self.GF)
        for x in x_s:
            z = z * galois.Poly([self.GF(1), -x], self.GF)
        z_s = self.evaluate_at_trusted_setup(z, self.trusted_setup_g2)
        lhs = pairing(z_s, proof)
        i_s = self.evaluate_at_trusted_setup(points_polynomial, self.trusted_setup_g1)
        rhs = pairing(G2, add(commitment, neg(i_s)))
        print("isProofValid: ", eq(lhs, rhs))
        assert eq(lhs, rhs), "The proof is invalid"


if __name__ == '__main__':
    kzg = KZGCommitment(curve_order)
    vector = [10, 20, 36, 50, 90]  # example vector
    print("Commiting Vector: ", vector)
    polynomial = kzg.vector_to_polynomial(vector)
    commitment = kzg.commit_polynomial(polynomial)
    points = [(0, 10), (1, 20)]  # element 10 at index 0 and element 20 at index 1
    proof = kzg.generate_proof(polynomial, points)
    print("Proof: ", proof)
    kzg.verify_proof(commitment, points, proof)
