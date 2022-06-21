---
myst:
  enable_extensions: ["colon_fence","amsmath"]
  dmath_allow_labels: True
---


# The Theodorsen Model

The Theodorsen model is a classical frequency-domain model of unsteady lift, derived analytically using unsteady potential flow theory. Because of this assumption, it is not expected to be accurate for low Reynolds numbers flows. Nevertheless, it is highly important to the project, because:

1. It is a well established model, often used in research problems related to unsteady lift <span style="color: red;">citation</span>.
2. It is formulated in a very physically meaningful way, with different terms directly related to specific flow phenomena.
3. It can be used to generate reliable unsteady lift data for high Reynolds number flows.
4. It can be used as a test problem for the system identification methodology developed for the project.

## Formulation

As shown in the previous report, in {cite}`VonKarman1938` Von Kármán arrives at the following form for the usteady lift of an airfoil undergoing a general motion:

:::{math}
:label: VonKarmanModel
 L = -\rho\frac{d}{dt}\int_{-1}^1\gamma_0(x)xdx + \rho U \Gamma_0 + \rho U\int_1^{\infty}\frac{\gamma(\xi)d\xi}{\sqrt{\xi^2-1}} 
:::

where we can identify the following terms:
* $L_1 = -\rho (d/dt)\int_{-1}^1\gamma_0(x)xdx$ is the added mass term.
* $L_0 = \rho U \Gamma_0$ is the quasi-steady lift.
* $L_2 = \rho U\int_1^{\infty}\gamma(\xi)d\xi/\sqrt{\xi^2-1}$ is the contribution of the wake to the overall lift.

The coordinate system used in the formula is shown in {numref}`VonKarmanCoordinates` 
```{figure} images/VonKarmanCoordinates.png
---
height: 200px
name: VonKarmanCoordinates
---
Coordinate system used in the Von Kármán model {eq}`VonKarmanModel` (and also in the Theodorsen model <span style="color: red;">citation</span>).
```

In equation {eq}`VonKarmanModel` $\Gamma_0$ is the quasi-steady circulation, $\gamma_0(x)$ the vorticity distribution on the surface of the airfoil and $\gamma(\xi)$ the vorticity distribution of the wake.

The Theodorsen model applies this framework to a thin airfoil of chord $c$ (and half-chord $b = c/2$) in a potential flow of free-stream velocity $U_\infty$ undergoing a pitching (change of angle of attack) and plunging (normal to t) motion. For purely sinusoidal motion with nondimensionalised frequency $k = \omega b / U_\infty$, the resulting lift can be formulated as <span style="color: red;">citation</span>:
:::{math}
	C_L = C_1 \left(\ddot{h}+\dot{\alpha}-a\ddot{\alpha} \right) + C_2 \left(\alpha + \dot{h} + \dot{\alpha}\left(\frac{1}{2}-a \right) \right)C(k)
:::

where:
* $\alpha$ is the angle of attack of the airfoil (along with its derivatives $\dot \alpha$ and $\ddot \alpha$)
* $\dot h$ is the vertical (normal to the free-stream flow) velocity of the airfoil (with its derivative $\ddot h$)
* $a$ is the pitch rotation centre, i.e. the point around which the pitching movement is performed, expresed in the frame shown in {numref}`VonKarmanCoordinates`
* $C_1$ is a coefficient related to added-mass effects, equal to $\pi$ for thin airfoil theory
* $C_2$ is a coefficient related to quasi-steady lift, equal to $2\pi$ for thin airfoil theory
* $C(k)$ is the frequency-domain Theodorsen function

In the expression for $C_L$, two distinct terms can be identified:
:::{math}
	C_L = C_L^{AM} + C_L^{QS}C(k)
:::
where $C_L^{AM}$ is the added-mass lift, corresponding to inertial forces due to air being displaced by the airfoil's movement, while $C_L^{QS}$ is the quasi-steady lift, corresponding to a steady-state flow condition around the airfoil. Contrary to the added-mass effects, which are instant, the quasi-steady lift experiences an attenuation and phase shift due to unsteady flow effects, which can be related to the vorticity of the wake behind the airfoil. The Theodorsen function was derived to model this exact phenomenon. It is analytically formulated as:
:::{math}
C(k) = \frac{H_1^{(2)}(k)}{H_1^{(2)}(k) + jH_0^{(2)}(k)}
:::

where $H_\nu^{(2)}(k)$ are Hankel functions, defined by an expression using Bessel functions, $H_\nu^{(2)}(k) = J_\nu + jY_\nu$. Moreover, the Theodorsen fuction can be generalised to a Laplace-domain transfer function of the nondimensionalised Laplace variable $\overline{s} = sb/U_\infty$  <span style="color: red;">citation</span>. Thanks to this, the Theodorsen model can be used to compute unsteady lift responses to arbitrary pitching and plunging motion <span style="color: red;">citation</span>.

## The model as a dynamical system

The effect of the Theodorsen function on the vortical lift is of delay and scaling on changes in angle of attack. In particular, it is common to consider the vertical velocity $\dot h$ and $\alpha$ together as an effective angle of attack $\alpha_e$. 

## Summary

As linear approximations of the Theodorsen function in the time domain already exist in literature {cite}`brunton2013empirical`, we have 
decided to use them as a first approach to generate the data necessary for the system identification. Such representation allow to construct a linear Multi-Input-Single-Output system (MISO), which has been implemented in a python script presented in <span style="color: red;">section</span>. That data has then been used to extract a first data-driven approximation of the $C(k)$ function in <span style="color: red;">section</span> using the algorithm presented in <span style="color: red;">section</span>.