# Task Formulation

We will deal with the solution of the non-stationary diffusion equation with the variable
$U = U(x, y, z, t):$
$∂U/∂t - div (D · grad U) = f(x, y, z, t), (x,y,z)$ belongs $Ω = [0; 1]^3$, t on the  $[0, T]$

$U(x, y, z, t) = g(x, y, z)$ on the edge of the region $∂Ω$

$U(x, y, z, 0) = 0$ at the initial time

We define the final moment of time as $T = 1$

In the problem, we will use the diagonal tensor $D:$

```math
D=\begin{pmatrix}
d_x & 0 & 0 \\
0 & d_y & 0 \\
0 & 0 & d_z 
\end{pmatrix}
```

где $d_x = 0.25, d_y = 0.15, d_z = 0.1$

Also set:

$g(x, y, z) = 0$

$f(x, y, z) = (d_x+d_y+d_z)·π²·sin(πx)sin(πy)sin(πz)$	 
This equation has analytic solution:
$U_analityc = sin(πx)sin(πy)sin(πz)·(1 - exp(-(d_x+d_y+d_z)·π²·t))$

# Sampling
Building a parallelepiped discretization of our region.

Set $Nx, Ny, Nz > 1$ - the number of nodes that will fit along the axis $Ox, Oy and Oz$, respectively.
Then we define the grid step $Δx = \frac{1}{(Nx-1)}, Δy = \frac{1}{(Ny-1)}, Δz = \frac{1}{(Nz-1)}. $
Also, define the time step $Δt$.
Define $V_{ijk}$ node of grid with coordinates $x_i = i·Δx, y_j = j·Δy, z_k = k·Δz.$

We will describe the discrete function $[U]^h$ at the time $nΔt$ by its degrees of freedom, which we will place at the nodes of the grid, and the degree of freedom at the node V_ijk will be denoted as

$U_{ijk}^n, 0 ⩽ i ⩽ Nx-1, 0 ⩽ j ⩽ Ny-1, 0 ⩽ k ⩽ Nz-1.$

We discretize our equation in space by the finite difference method, and in time by the explicit Euler scheme.

Note that for such a discretization, the time step must satisfy the Courant condition:

$Δt < \frac{0.5}{ \frac{d_x } {(Δx)^2} + \frac{d_y}{(Δy)^2} + \frac{d_z}{(Δz)^2} }$

Introduce discrete operators of the second spatial derivatives:

$Lx_{ijk} U^n = \frac{ U_{(i-1)jk}^n - 2·U_{ijk}^n + U_{(i+1)jk}^n }{Δx^2}$

$Ly_{ijk} U^n = \frac{ U_{i(j-1)k}^n - 2·U_{ijk}^n + U_{i(j+1)k}^n }{Δy)^2}$

$Lz_{ijk} U^n = \frac{ U_{ij(k-1)}^n - 2·U_{ijk}^n + U_{ij(k+1)}^n }{Δz^2}$

Also, let

$f_{ijk}^n = f(i·Δx, j·Δy, k·Δz, nΔt)$ и $g_{ijk} = g(i·Δx, j·Δy, k·Δz)$

Then we have the following numerical scheme:
$U_{ijk}^{(n+1)} = U_{ijk}^n + Δt·(f_{ijk}^n + d_x·Lx_{ijk} U^n + d_y·Ly_{ijk} U^n + d_z·Lz_{ijk} U^n), 


if


1 ⩽ i ⩽ Nx-2, 1 ⩽ j ⩽ Ny-2, 1 ⩽ k ⩽ Nz-2$

$U_{ijk}^{(n+1)} = g_{ijk}$, if $(i%(Nx-1)) · (j%(Ny-1)) · (k%(Nz-1)) = 0$,
где $x%y$ - operation of taking the remainder of a division $x$ на $y$.

# How to open results in ParaView
1. Generate ```.vtk``` file [guide](https://github.com/artiebears13/diffusion_nonstat_MFD_explicit/blob/main/CPP/ReadMe.md)
2. Install [ParaView](https://www.paraview.org/)
3. In ParaView open file ```numeric.vtk```
4. Choose u:

![choose u](https://github.com/artiebears13/diffusion_nonstat_MFD_explicit/blob/main/examples/choose_u.png)

6. Choose clip and cut 

![clip](https://github.com/artiebears13/diffusion_nonstat_MFD_explicit/blob/main/examples/choose_clip.png)

# Examples
![example1](https://github.com/artiebears13/diffusion_nonstat_MFD_explicit/blob/main/examples/clip_perp_surface.png)

![example2](https://github.com/artiebears13/diffusion_nonstat_MFD_explicit/blob/main/examples/clip_surface_edges.png)

 
