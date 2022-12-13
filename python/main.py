import mfd_solver as solver
import numpy as np
from datetime import datetime
import time


def main():
    a = solver.Solver(11, 11, 11, 1000)
    start_time = datetime.now()
    our = a.solve()
    print(datetime.now() - start_time)
    anal = a.get_analytic_solve()
    # print(a.u_actual.shape)
    a.print()
    print('mistake : ', solver.max_mistake(our, anal))
    a.plot()
    a.plot_analytic()
    # print(a.u_actual.shape)



main()