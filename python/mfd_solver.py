import numpy as np

'''
Class Solver: solving equation:
∂U/∂t - div (D · grad U) = f(x, y, z, t) on [0; 1]³ and [0, T]

Using Diag tensor D: 
	    | d_x  0   0  |
	D = |  0  d_y  0  |
	    |  0   0  d_z |
where d_x = 0.25, d_y = 0.15, d_z = 0.1

Also we have 
g(x, y, z) = 0 - boundary condition
f(x, y, z) = (d_x+d_y+d_z)·π²·sin(πx)sin(πy)sin(πz) - right side
discretization:
n,m,k - number of nodes for x,y,z directions
our equation have analytic solve:
U_analityc = sin(πx)sin(πy)sin(πz)·(1 - exp(-(d_x+d_y+d_z)·π²·t))
Steps:
hx = 1/(Nx-1), hy = 1/(Ny-1), hz = 1/(Nz-1).

also dt is time step has condition:
dt < 0.5 / { d_x / (hx)² + d_y / (hy)² + d_z / (hz)² }

Solving:
Let us introduce discrete operators of secondary glass derivatives:
	Lx_ijk U^n = ( U_(i-1)jk^n - 2·U_ijk^n + U_(i+1)jk^n )/(hx)²
	Ly_ijk U^n = ( U_i(j-1)k^n - 2·U_ijk^n + U_i(j+1)k^n )/(hy)²
	Lz_ijk U^n = ( U_ij(k-1)^n - 2·U_ijk^n + U_ij(k+1)^n )/(hz)²
Also let: 
    f_ijk^n = f(i·hx, j·hy, k·hz, ndt) и g_ijk = g(i·hx, j·hy, k·hz)
    
And, finally our numerical scheme:
U_ijk^(n+1) = U_ijk^n + Δt·(f_ijk^n + d_x·Lx_ijk U^n + d_y·Ly_ijk U^n + d_z·Lz_ijk U^n), 
if 1 ⩽ i ⩽ Nx-2, 1 ⩽ j ⩽ Ny-2, 1 ⩽ k ⩽ Nz-2 (not boundary)
	U_ijk^(n+1) = g_ijk, if (i%(Nx-1)) · (j%(Ny-1)) · (k%(Nz-1)) = 0, (boundary)


'''


class Solver:
    D = [0.25, 0.15, 0.1]

    def __init__(self, n, m, k, t):
        self.n = n
        self.m = m
        self.k = k
        self.t = t
        self.dt = 1 / t
        self.u_actual = np.zeros((n, m, k))
        self.u_next = np.zeros((n, m, k))
        self.hx = 1 / (n - 1)
        self.hy = 1 / (m - 1)
        self.hz = 1 / (k - 1)

    def check_boundary(self, i, j, k):
        if i == 0 or i == self.n \
                or j == 0 or i == self.m \
                or k == 0 or k == self.k:
            return True
        else:
            return False

    def fill_boundaries(self):
        bi = [0, self.n - 1]
        bj = [0, self.m - 1]
        bk = [0, self.k - 1]
        for i in bi:
            for j in bj:
                for k in bk:
                    self.u_actual[i, j, k] = 0

    def Lx(self, i, j, k):
        return (self.u_actual[i - 1][j][k] - 2 * self.u_actual[i][j][k] + self.u_actual[i + 1][j][k]) \
            / (self.hx ** 2)

    def Ly(self, i, j, k):
        return (self.u_actual[i][j - 1][k] - 2 * self.u_actual[i][j][k] + self.u_actual[i][j + 1][k]) \
            / (self.hy ** 2)

    def Lz(self, i, j, k):
        return (self.u_actual[i][j][k - 1] - 2 * self.u_actual[i][j][k] + self.u_actual[i][j][k + 1]) \
            / (self.hz ** 2)

    # right side
    def f(self, x, y, z):
        return ((self.D[0] + self.D[1] + self.D[2])
                * (np.pi ** 2) * np.sin(np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z))

    # calculating u on next t step
    def get_next_u(self):
        for i in range(1, self.n - 1):
            for j in range(1, self.m - 1):
                for k in range(1, self.k - 1):
                    self.fill_boundaries()
                    # print('NOT boundary i={}, j={}, k={}'.format(i, j, k))
                    # print('F: {}'.format(self.f(i * self.hx, j * self.hy, k * self.hz)))
                    self.u_next[i, j, k] = self.u_actual[i, j, k] + \
                                           self.dt * (self.f(i * self.hx, j * self.hy, k * self.hz)
                                                      + self.D[0] * self.Lx(i, j, k)
                                                      + self.D[1] * self.Ly(i, j, k)
                                                      + self.D[2] * self.Lz(i, j, k))
        # print(self.u_next)

    # putting u on next t step in u_actual
    # REMEMBER: actual u in u_actual
    def calculate_next_u(self):
        self.get_next_u()
        # tmp = self.u_next
        # self.u_next = self.u_actual
        self.u_actual = self.u_next

    # solver function
    # call it to get answer
    def solve(self):
        for i in range(1, self.t):
            self.calculate_next_u()
        return self.u_actual

    def get_analytic_u(self, x, y, z):
        return (np.sin(np.pi * x)
                * np.sin(np.pi * y)
                * np.sin(np.pi * z)
                * (1 - np.exp(-((self.D[0] + self.D[1] + self.D[2]) * np.pi ** 2 * 1))))

    def get_analytic_solve(self):
        u = np.zeros((self.n, self.m, self.k))
        for i in range(1, self.n - 1):
            for j in range(1, self.m - 1):
                for k in range(1, self.k - 1):
                    self.fill_boundaries()
                    # print('NOT boundary i={}, j={}, k={}'.format(i, j, k))
                    # print('F: {}'.format(self.f(i * self.hx, j * self.hy, k * self.hz)))
                    u[i, j, k] = self.get_analytic_u(i * self.hx, j * self.hy, k * self.hz)
        return u



    def print(self):
        print('n = {}, m = {}, k = {}'.format(self.n, self.m, self.k))
        print('hx = {}, hy = {}, hz = {}'.format(self.hx, self.hy, self.hz))
        print('U: ')
        # print(self.u_actual)
        print('dt=',self.dt,' < ',0.5/(self.D[0]/(self.hx**2)+self.D[1]/(self.hy**2)+self.D[2]/(self.hz**2)))



def max_mistake(a,b):
    mist = 0
    our = a
    anal = b
    for i in range(len(a)):
        for j in range(len(a)):
            for k in range(len(a)):
                mist = max(np.abs(our[i, j, k] - anal[i, j, k]), mist)
    return mist


