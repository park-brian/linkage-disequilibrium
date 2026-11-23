# Fixture for src/pyfmt.js npStr4: numpy's str() of floats under ldsc's np.set_printoptions(precision=4), run in the LDlink
# oracle container (numpy 2.4.3).  usage (tools/ldlink-local): docker compose exec -T backend /usr/bin/python3.13 - \
#   < verify/reference/numpy_s4.py > test/fixtures/numpy_s4.json
import json, sys, random
import numpy as np
np.set_printoptions(linewidth=1000)
np.set_printoptions(precision=4)
s = lambda x: str(np.matrix(x)).replace('[', '').replace(']', '').strip()  # ldsc regressions.s + remove_brackets
random.seed(7)
vals = [0.03125, 0.15625, 0.34375, 1.03125, 0.00015625, 0.09375, 12.40625, 0.28125, 0.0, -0.0, 1.0, -1.0, 2.0, 10.0, 100.0, 0.5, 0.25, 0.1, 0.0001, 0.00009999, 0.00010001, 99999999.0, 100000000.0, 123456789.0,
        1e-20, 5e-05, 1.5e-05, 0.00012345, 0.2641, 0.26405722123648157, 0.1017613631357963, 1.2103, 0.9508933165947758, 0.047, -0.1715,
        0.00005, 0.12345, 0.12355, 0.99995, 9.99995, 1.00005, 0.000149995, 3.14159265, -2.71828, 1e8, 1e16, 1.234e-300, 12.5, 7.0,
        float('nan'), float('inf'), float('-inf')]
for _ in range(1500):
    e = random.uniform(-12, 10); v = random.choice([1, -1]) * random.random() * 10 ** e
    vals.append(v)
for _ in range(300):
    vals.append(round(random.uniform(-50, 50), random.randint(0, 6)))
out = [[repr(float(v)), s(v)] for v in vals]
json.dump({"numpy": np.__version__, "cases": out}, sys.stdout)
