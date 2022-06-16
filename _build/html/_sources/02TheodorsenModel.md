---
myst:
  enable_extensions: ["colon_fence","amsmath"]
  dmath_allow_labels: True
---


# The Theodorsen Model

## Relevant mathematics

### Theodorsen 

As shown in the previous report, in {cite}`VonKarman1938` Von Kármán arrives at the following form for the usteady lift of an airfoil undergoing a general motion:

:::{math}
:label: VonKarmanModel
 L = -\rho\frac{d}{dt}\int_{-1}^1\gamma_0(x)xdx + \rho U \Gamma_0 + \rho U\int_1^{\infty}\frac{\gamma(\xi)d\xi}{\sqrt{\xi^2-1}} 
:::

where we can identify the following terms:
* $L_1 = -\rho (d/dt)\int_{-1}^1\gamma_0(x)xdx$ is the added mass term.
* $L_0 = \rho U \Gamma_0$ is the quasi-steady lift.
* $L_2 = \rho U\int_1^{\infty}\gamma(\xi)d\xi/\sqrt{\xi^2-1}$ is the contribution of the wake to the overall lift.

The coordiante system used in the formula is shown in {numref}`VonKarmanCoordinates` 
```{figure} images/VonKarmanCoordinates.png
---
height: 200px
name: VonKarmanCoordinates
---
Coordiante system for the Von Kármán model {eq}`VonKarmanModel`
```

In equation {eq}`VonKarmanModel` $\Gamma_0$ is the quasi-steady circulation, $\gamma_0(x)$ the vorticity distribution on the surface of the airfoil and $\gamma(\xi)$ the vorticity distribution of the wake.

In the case of a purely sinusoidal motion with frequency $k$ it is possible to write a simpler formula for the lift coefficient:
:::{math}
	C_L = C_1 \left(\ddot{h}+\dot{\alpha}-a\ddot{\alpha} \right) + C_2 \left(\alpha + \dot{h} + \dot{\alpha}\left(\frac{1}{2}-a \right) \right)C(k) = C_L^{AM} + C_L^{QS}C(k)
:::

where $C(k)$ is the Theodorsen function, defined as:
:::{math}
C(k) = \frac{H_1^{(2)}(k)}{H_1^{(2)}(k) + jH_0^{(2)}(k)}
:::

where $H_\nu^{(2)}(k)$ are Hankel functions, defined by an expression using Bessel functions, $H_\nu^{(2)}(k) = J_\nu + jY_\nu$. The variables of the system are the angle of attack $\alpha$ and the vertical velocity $\dot h$ with their respective derivatives $\dot \alpha$, $\ddot \alpha$ and $\ddot h$. The pitch rotation centre is called $a$, and $C_1$ and $C_2$ are real coefficients for the added mass ($C_L^{AM}$) and the quasi-steady vortical lift ($C_L^{QS}$) respectively. For the classical model $C_1 = \pi$ and $C_2 = 2\pi$.

The effect of the Theodorsen function on the vortical lift is of delay and scaling on changes in angle of attack. In particular, it is common to consider the vertical velocity $\dot h$ and $\alpha$ together as an effective angle of attack $\alpha_e$. 

As linear approximations of the Theodorsen function in the time domain already exist in literature {cite}`brunton2013empirical`, we have 
decided to use them as a first approach to generate the data necessary for the system identification. Such representation allow to construct a linear Multi-Input-Single-Output system (MISO), which has been implemented in a pythpn script presented in <span style="color: red;">section</span>. That data has then been used to extract a first data-driven approximation of the $C(k)$ function in <span style="color: red;">section</span> using the algorithm presented in <span style="color: red;">section</span>.

### The Wagner function 

The Wagner function can be seen as the time domain representation of the Theodorsen function. In particular it can be defined as the time-response of an airfoil undergoing a sudden change in angle of attack for time $t=0$, i.e. the step response of the system to a change in $\alpha_e$ {cite}`dawson2022improved`.

Calling $\phi(t)$ the Wagner function, the analytical computation can be carried out as:
:::{math}
\begin{equation}
        \phi(t) = \mathcal{L}^{-1}\left(\frac{C(s)}{s}\right)
\end{equation}
:::

where $\mathcal{L}^{-1}$ is the inverse Laplace transform, and $C(s)$ is the Theodorsen function. Defining the real and imaginary part of the Theodorsen function as $G_1(k)$ and $G_2(k)$ respectively, it is possible to compute the Wagner function as {cite}`peters2008two`:

:::{math}
\begin{eqnarray}
  \phi(t) &=& \frac{1}{2} + \frac{2}{\pi}\int_0^\infty \frac{1}{k}\left(G_1(k) - \frac{1}{2}\right)\sin(kt)dk :label: small_t\\
  &=& 1 + \frac{2}{\pi}\int_0^\infty \frac{1}{k}G_2(k)\cos(kt)dk\label{eqn:int_large_t} :label: large_t
\end{eqnarray}
:::

The first equation has better numerical properties when integrated for small time steps, while the second is more stable for large ones.

The computation of the Wagner function in this fashion is rather cumbersome, as it requires to integrate numerically the Theodorsen function. In particular for large times the trigonometric terms in the above equations make the integrals slowly convergent and require a very fine discretization of the domain, making the computational time quite elevated. As a first step in our project we ahave then decided to test the selected system identification algorithm, presented in <span style="color: red;">section</span>, to reproduce the Wagner function, in an attempt to replicate the study conducted in {cite}`dawson2022improved` and gain confidence with the numerical methods involved. The results are presented in <span style="color: red;">section</span>.

