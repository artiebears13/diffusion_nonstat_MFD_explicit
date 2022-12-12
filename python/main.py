import mfd_solver as solver
import numpy as np


def main():
    a = solver.Solver(21, 21, 21, int(1 / 0.002) + 1)
    our = a.solve()
    anal = a.get_analytic_solve()
    a.print()
    # print(a.u_actual.shape)
    print('mistake : ', solver.max_mistake(our, anal))


main()