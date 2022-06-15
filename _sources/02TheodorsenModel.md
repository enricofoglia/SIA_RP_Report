---
myst:
  enable_extensions: ["colon_fence","amsmath"]
  dmath_allow_labels: True
---


# The Theodorsen Model

## Relevant mathematics

As shown in the previous report, Von Kármán arrives at the following form for the usteady lift of an airfoil undergoing a general motion:

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

In equation {eq}`VonKarmanModel` $\Gamma_0$ is the quas-steady circulation, $\gamma_0(x)$ the vorticity distribution on the surface of the airfoil and $\gamma(\xi)$ the vorticity distribution of the wake.

In the case of a purely sinusoidal motion with frequency $k$ it is possible to write a simpler formula for the lift coefficient:
:::{math}
	C_L = C_1 \left(\ddot{h}+\dot{\alpha}-a\ddot{\alpha} \right) + C_2 \left(\alpha + \dot{h} + \dot{\alpha}\left(\frac{1}{2}-a \right) \right)C(k) = C_L^{AM} + C_L^{QS}C(k)
:::

where $C(k)$ is the Theodorsen function, defined as:
:::{math}
C(k) = \frac{H_1^{(2)}(k)}{H_1^{(2)}(k) + jH_0^{(2)}(k)}
:::

where $H_\nu^{(2)}(k)$ are Hankel functions, defined by an expression using Bessel functions, $H_\nu^{(2)}(k) - J_\nu + jY_\nu$. The variables of the system are the angle of attack $\alpha$ and the vertical velocity $\dot h$ with their respective derivatives $\dot \alpha$, $\ddot \alpha$ and $\ddot h$. The pitch rotation centre is called $a$, and $C_1$ and $C_2$ are real coefficients for the added mass ($C_L^{AM}$) and the quasi-steady vortical lift ($C_L^{QS}$) respectively. For the classical model $C_1 = \pi$ and $C_2 = 2\pi$.