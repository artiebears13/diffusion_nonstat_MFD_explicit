import mfd_solver as solver
import numpy as np


def main():
    a = solver.Solver(11, 11, 11, 1000)
    our = a.solve()
    anal = a.get_analytic_solve()
    # print(a.u_actual.shape)
    a.print()
    print('mistake : ', solver.max_mistake(our, anal))
    a.plot()
    a.plot_analytic()
    # print(a.u_actual.shape)



main()