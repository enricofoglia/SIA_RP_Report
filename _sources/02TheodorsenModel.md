---
myst:
  enable_extensions: ["colon_fence","amsmath"]
  dmath_allow_labels: True
---

(section:theodorsen_model)=
# The Theodorsen Model

The Theodorsen model is a classical frequency-domain model of unsteady lift, derived analytically using unsteady potential flow theory. Because of this assumption, it is not expected to be accurate for low Reynolds numbers flows. Nevertheless, it is highly important to the project, because:

1. It is a well established model, often used in research problems related to unsteady lift {cite}`theodorsen1979` {cite}`brunton2013empirical` {cite}`brunton_reduced-order_2013`.
2. It is formulated in a very physically meaningful way, with different terms directly related to specific flow phenomena.
3. It can be used to generate reliable unsteady lift data for high Reynolds number flows much more efficiently that CFD simulations.
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
Coordinate system used in the Von Kármán model {eq}`VonKarmanModel` (and also in the Theodorsen model {cite}`brunton2013empirical`).
```

In equation {eq}`VonKarmanModel` $\Gamma_0$ is the quasi-steady circulation, $\gamma_0(x)$ the vorticity distribution on the surface of the airfoil and $\gamma(\xi)$ the vorticity distribution of the wake.

```{figure} images/unsteady_airfoil.png
---
width: 400px
name: unsteady_airfoil
---
Visualisation of the basic elements of the unsteady airfoil problem.
```

The Theodorsen model applies this framework to a thin airfoil of chord $c$ (and half-chord $b = c/2$) in a potential flow of free-stream velocity $U_\infty$ undergoing a pitching (change of angle of attack) and plunging (normal to t) motion, as shown in {numref}`unsteady_airfoil`. The lift can be nondimensionalised using the notion of the lift coefficient:

:::{math}
	C_L = \frac{2 L}{\rho U_\infty^2 c}
:::

Note that in the Von Kármán framework, the lift is expressed with respect to unit span (since the model is 2D), therefore the lift is divided only by the chord $c$ instead of the reference area.

 For purely sinusoidal motion with nondimensionalised frequency $k = \omega b / U_\infty$, the resulting lift can be formulated as {cite}`brunton2013empirical`:

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

where $H_\nu^{(2)}(k)$ are Hankel functions, defined by an expression using Bessel functions, $H_\nu^{(2)}(k) = J_\nu + jY_\nu$. Moreover, the Theodorsen fuction can be generalised to a Laplace-domain transfer function of the nondimensionalised Laplace variable $\overline{s} = sb/U_\infty$  <span style="color: red;">citation</span>. Thanks to this, the Theodorsen model can be used to compute unsteady lift responses to arbitrary pitching and plunging motion {cite}`brunton2013empirical`.

## The model as a dynamical system

The Theodorsen model can be interpreted as a dynamical system, accepting $\ddot \alpha$ and $\ddot h$ as input and outputting $C_L$, with an internal state evolving in response to the input and itself. As explained before, 3 subsystems can be distinguished, as visualised in {numref}`Theodorsen_diagram`:

1. The added-mass lift, accepting $\ddot \alpha$ and $\ddot h$ as input and outputting $C_L^{AM}$.
2. The quasi-static lift, accepting $\ddot \alpha$ and $\ddot h$ as input and outputting $C_L^{QS}$.
3. The Thedorsen function, accepting $C_L^{QS}$ as input and outputting the attenuated circulatory lift $C_L^{circ}$.

```{figure} images/Theodorsen_diagram.png
---
width: 500px
name: Theodorsen_diagram
---
Diagram of the subsystems of the dynamical system defined by the Theodorsen model (modified from {cite}`brunton2013empirical`).
```

The input was defined using the second derivatives of $\alpha$ and $h$, because a dynamical system that need to internally differentiate its input has undesireable properties, both analytically and numerically. The choice of second derivatives is also often useful in practical applications, since physically speaking, the movement of the airfoil with respect to the free-stream flow is related to pitching moments (which cause an angular acceleration that translates to $\ddot{\alpha}$) and normal forces (which cause the vertical acceleration $\ddot{h}$).

After splitting the system into 2 single-input, single-output (SISO) systems, Bode plots can be created for both, as shown in {numref}`Theodorsen_bode`.

```{figure} images/Theodorsen_bode.png
---
width: 500px
name: Theodorsen_bode
---
Bode plot of the Theodorsen model with respect to both $\ddot h$ and $\ddot \alpha$ inputs, the latter depending also on the pitching center location (copied from {cite}`brunton2013empirical`).
```

The Bode plot provides insight into the dynamics of the system, as described in {cite}`brunton_reduced-order_2013` and {cite}`brunton2013empirical`. For the $\ddot{\alpha} \rightarrow C_L$ system, in the low-frequency limit, the dominating effects are those of $\alpha$, since $\dot{\alpha}$ and $\ddot{\alpha}$ are significantly smaller for low frequencies. The system simplifies to the classical steady linear model of lift:

:::{math}
C_L \approx C_2 \alpha C(k) \approx C_2 \alpha
:::

Since the input to the system is $\ddot{\alpha}$, the system dynamics consists mainly of integrating the input twice (which corresponds to the magnitude slope of -40dB/decade, since integrating the input always reduces the slope by 20dB/decade) and applying a gain of $C_2$. The phase of $C_L$ is $180 [deg]$ behind $\alpha$, because each integration of a sine wave moves its phase back by $90[deg]$. For the $\ddot(h) \rightarrow C_L$ subsystem, the low frequency limit is:

:::{math}
C_L \approx C_2 \dot{h} C(k) \approx C_2 \dot{h}
:::

This shows that the plunging velocity has the same effect on the lift as the angle of attack, because of its direct influence on the apparent velocity of the incoming flow (this remains true for other frequencies as well). Because of this, an equivalent angle of attack is sometimes used as an alternative state for $\alpha$ and $\dot{h}$ (which has the benefit of being observable, contrary to  $\alpha$ and $\dot{h}$ whose effects are not distinguishable):

:::{math}
\alpha_e = \alpha + \dot{h}
:::

For the low frequency $\ddot(h) \rightarrow C_L$ system, since the dominant action is integration (but only once) of the input, the slope of the Bode plot is -20dB/decade, and the phase is $-90[deg]$, due to considerations analogous to the $\ddot{\alpha} \rightarrow C_L$ system.

In the high frequency limit, the $\ddot{\alpha} \rightarrow C_L$ system is dominated by the $\ddot{\alpha}$ term of the added-mass force, resulting in the form:

:::{math}
C_L \approx - C_1 a \ddot{\alpha} C(k) \approx - C_1 a  \ddot{\alpha}/4
:::

Depending on the location of the pitching centre $a$, 3 dynamics are possible:

1. If the pitch centre is before the half-chord, i.e. $a<0$, the system becomes $C_L \approx - C_1 a \ddot{\alpha}/4 = C  \ddot{\alpha}$, where $C > 0$. This means that $C_L$ is proportional to $\ddot{\alpha}$, with a gain asymptotically approaching $C = -C_1 a /4 > 0$). Because $C>0$, C_L is in phase with $\ddot{\alpha}$ (hence the phase of the Bode plot approaching $0[deg]$), meaning that it is also in opposite phase to $\alpha$. This means that during a rapid pitch up movement, the initial resulting lift is actually directed downwards.
2. If the pitch centre is behind the half-chord, i.e. $a>0$, the system becomes $C_L \approx - C_1 a \ddot{\alpha}/4 = C  \ddot{\alpha}$, where $C < 0$. This means that $C_L$ is proportional to $\ddot{\alpha}$, with a gain asymptotically approaching $C = -C_1 a /4 < 0$). Because $C<0$, C_L is in opposite phase to $\ddot{\alpha}$ (hence the phase of the Bode plot approaching $-180[deg]$), meaning that it is also in phase with $\alpha$. This means that during a rapid pitch up movement, the initial resulting lift is directed upwards.
3. If the pitch centre is exactly at the half-chord, i.e. $a=0$, the $\ddot{\alpha}$ term gets cancelled, meaning that $\dot{\alpha}$ terms are dominant instead:
:::{math}
C_L \approx C_1 \dot{\alpha} + C_2 \dot{\alpha} (0.5-a) C(k) \approx C_1 \dot{\alpha} + C_2 \dot{\alpha}/8 = (C_1 + C_2/8) \dot{\alpha}
:::
This means that the dynamics of $C_L$ consist of a single integration of the input $\ddot{\alpha}$ and a constant gain of $C = (C_1 + C_2/8)$. Consistent with this, the Bode plot has a gain slope of -20dB/decade an phase asymptotically approaching $-90[deg]$.

The $\ddot{h} \rightarrow C_L$ system in the high frequency limit is dominated by the $\ddot{h}$ term of the added-mass force, resulting in the form:

:::{math}
C_L \approx C_1 \ddot{h}
:::

Analogously to the pitching dynamics with $a<0$, this leads to an asymptotic gain of $C_1$ and phase of $0[deg]$.

## Summary

The Thedorsen model is a useful tool for understanding and predicting unsteady lift at high Reynolds number flows. While it is difficult to implement in an exact manner, precise linear state-space approximations are available in the literature — one of them, originating from {cite}`brunton2013empirical`, is implemented and demonstrated in {ref}`section:theodorsen_implementation`, and later used throughout the report as a reference model for the system identification procedures in {ref}`section:results_theodorsen`.