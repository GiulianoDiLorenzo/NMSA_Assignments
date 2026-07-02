# Numerical Modeling and Simulation for Acoustics - Assignments

Course: Numerical Modeling and Simulation for Acoustics

Programme: MSc in Music and Acoustic Engineering, Politecnico di Milano

Work Team: Di Lorenzo Giuliano, Ouali Ernest

## Overview

Three MATLAB assignments covering the theoretical and numerical analysis of acoustic problems: from 1D finite element method (FEM), to leap-frog scheme for implicit equation, and study of flows in finite volume method (FVM).

## HW1 - Finite Element Method for a 1D Helmholtz-type problem

Solves the frequency-domain wave equation

$$
\begin{cases}
\rho(x)\omega u(x) + \dfrac{d}{dx}(\mu(x) \frac{du}{dx})(x) = f(x), \\
u(0) = g_D, \\
u(L) - \dfrac{i}{\rho(x)\omega} \dfrac{du}{dx}(L) = g_A.
\end{cases} \qquad
x \in (0, L)
$$

with variable stiffness $\mu(x)$ and density $\rho(x)$, and a Robin (impedance) condition at $x = L$.

- Weak formulation and Galerkin discretization with piecewise-linear basis functions, leading to the system $\omega M \cdot U - A \cdot U = F$.
- Mass and stiffness matrix assembly via Gauss-Legendre quadrature.
- Convergence verified against the exact solution $u_\text{ex}(x) = \sin(2\pi x)$: measured $L^2$ convergence order $\approx 2$, consistent with linear FEM.
- Sensitivity of accuracy to piecewise-constant $\mu(x)$ and $\rho(x)$ (material discontinuities), showing negligible influence of the exact material values on the error when the exact solution is fixed by construction.
- Frequency sweep over $\omega \in (0, 20)$ with 10 points-per-wavelength mesh sizing, gain/phase spectra, and identification of resonances for constant and piecewise-constant wave speed $c(x)$.
- 2D space-time visualization of the harmonic solution $v(x,t) = e^{i\omega t}u(x)$.

## HW2 - Webster's equation for 1D acoustic tubes

Solves Webster's horn equation for a 1D acoustic tube of variable cross-section $S(x)$

$$
S(x) \frac{\partial^2 φ}{\partial t^2} = \gamma ^2 \frac{\partial}{\partial x}(S(x) \frac{\partial φ}{\partial x}) + f, \quad
\gamma  = c/L
$$

with excitation at $x = 0$ and radiation at $x = L$.

- Explicit finite difference scheme, reduces to the standard leap-frog scheme for constant cross-section, extended to variable $S(x)$ via a second-order local polynomial approximation of the profile.
- Verified against a closed-form exact solution for both constant and variable cross-section, with $L^2$-error convergence: linear in time ($O(\Delta t)$), and space-refinement study showing the error plateaus unless the time step is refined first - a coupling effect between space and time discretization worth noting for anyone reusing this scheme.
- CFL stability condition derived and demonstrated empirically (unstable configuration blows up almost immediately).
- Applied to a **vocal tract simulation**: variable cross-section profiles for vowels `e` and `a`, driven by a glottal-pulse-like input, using a staggered finite difference scheme. Output shows pitch-periodic behavior consistent with vowel sounds.

## HW3 - Traffic flow as a scalar conservation law

Treats the LWR traffic flow model

$$
\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x} = 0
$$

as a scalar hyperbolic conservation law, using the analogy between vehicle density and 1D wave propagation.

- Two flux models compared: the standard linear speed-density relation $u(\rho) = u_\text{max}(1 - \frac{\rho}{\rho_\text{max}})$, and a logarithmic flux $f(\rho) = \rho \cdot log(\frac{\rho_\text{max}}{\rho})$ that better matches real traffic data.
- Finite volume discretization with the **Godunov method**: first-order (constant reconstruction) and second-order (linear reconstruction with a Monotonized Central slope limiter) schemes.
- Three canonical scenarios - Traffic Jam, Green Light, Traffic Flow - each with physically motivated boundary conditions, comparing shock formation/propagation and rarefaction wave behavior between flux models and scheme orders.
- Discussion of the accuracy/sharpness trade-off: the second-order scheme preserves discontinuities more faithfully, while the first-order scheme is more diffusive.
