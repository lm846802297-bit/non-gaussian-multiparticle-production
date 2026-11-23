from collections import Counter
from itertools import permutations
from math import factorial, sqrt

patterns = [
    (0, 0, 0, 0),
    (0, 0, 0, 1),
    (0, 0, 1, 1),
    (0, 1, 1, 2),
    (0, 1, 2, 3),
]

for pattern in patterns:
    occ = Counter(pattern)
    # Matrix element of (1/4!) sum_{p1...p4} a^dagger_p1 ... a^dagger_p4
    # between the vacuum and the normalized occupation state |{n_p}>.
    distinct_orderings = set(permutations(pattern))
    direct = 0.0
    for ordering in distinct_orderings:
        multiplicity = 1
        for mode, n in occ.items():
            multiplicity *= sqrt(factorial(n))
        direct += multiplicity / factorial(4)
    predicted = 1.0
    for n in occ.values():
        predicted /= sqrt(factorial(n))
    assert abs(direct - predicted) < 1e-14, (pattern, direct, predicted)
    print(pattern, f"direct={direct:.12g}", f"formula={predicted:.12g}", "PASS")
