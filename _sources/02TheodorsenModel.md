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
height: 150px
name: VonKarmanCoordinates
---
Coordiante system for the Von Kármán model {eq}`VonKarmanModel`
```