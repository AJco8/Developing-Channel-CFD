# Developing-Channel-CFD

## Abstract
A MATLAB program simulating a steady, two-dimensional, incompressible flow with constant properties by applying the finite volume method. The program uses a Mash Convergence study to find the optimal size for the Tri-diagonal Matrix Algorithm to discretize a non-dimensional model. 

## Results
![](./Figures/Developing_Flow.jpg)
![](./Figures/Vectors.jpg)

## Verification
This is verified by comparing the results to numeric predictions, a polynomial approximation of the boundary layer, and the commencement flow in a circular pipe. 

## Math
Continuity Equation
$$ \frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0 $$

Conservation of Momentum
$$ \rho \left[u \frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}\right]=-\frac{\partial p}{\partial x}+μ\left(\frac{\partial ^2 u}{\partial y^2}+\frac{\partial ^2 u}{\partial x^2}\right) $$

Conservation of Mass
$$\int_0^1u^* dy^*=1$$

Fully Developed Flow
$$ u=\frac{1}{2\mu} \frac{\partial P}{\partial x}\left(y^2-Hy\right) $$

### Boundary Layer Polynomial Approximation
$$ u^* = 2\eta - \eta^2 $$

$$ \eta = \frac{y}{\delta} $$

$$ \delta = \frac{3}{2} \left(1 - \frac{U_i}{U}\right) $$

### Commencement Flow in Circular Pipe
$$ u(t,r) = -\frac{R^2}{4\mu} \frac{dp}{dz} \left[1 - \left(\frac{r}{R}\right)^2 - \sum_{n=1}^{\infty} \frac{8J_0\left(\frac{\alpha_n r}{R}\right)}{\alpha_n^3 J_1(\alpha_n)} e^{-\left(\frac{\alpha_n^2 \nu t}{R^2}\right)} \right] $$
