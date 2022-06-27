(section:identification_lift)=
# System identification for unsteady lift

## Literature considerations

This subsection highlights and discusses some works with elements highly relevant to our project.

In {cite}`brunton2013empirical` a an approximate state-space representation of the Theodorsen has been obtained by combining the already linear terms of the model with a balanced linear model of the Theodorsen function, as shown in {ref}`section:theodorsen_implementation`. It was computed by evaluating the Theodorsen function at 1200 logarithmically spaced reduced frequencies, and fitted to a state-space model of order 11 using black box frequency domain methods implemented in the MATLAB _fitfrd_ function. Finally, balanced truncation {cite}`brunton2019data` was used to reduce the order of the model to 4 (models of higher order did not present significantly lower error).

In {cite}`dawson2022improved`, it was proposed to apporoximate the Wagner function by the output of a nonlinear dynamic system with $C_L$ (i.e. the output itself) as the only state. The system had no input (only an initial condition). The base of the nonlinear model were just polynomials of increasing order, resulting in a model of the form $\dot{C}_L = \sum C_L^n$. Importantly, only terms up to the order 6 contributed to an increase in accuracy. Another approach presented added $\dot{C}_L$ as a state, thus increasing the order of the nonlinear model to 2. Overall, the SINDy-based models were shown to yield better results than many previous approximations.

3 successive works by Brunton et al.: {cite}`brunton_unsteady_2012`, {cite}`brunton_reduced-order_2013` and {cite}`brunton_state-space_2014`, focus on creating a linear state-space model of unsteady lift, especially for low Reynolds number cases. The obtained models were linearised around a given angle of attack, leading to maintaining accuracy only in a certain range of angles of attack, most notably limited by stall phenomoena. Their general form was simillar to that of the LTI approximation to the Theodorsen model shown in {cite}`brunton2013empirical` and described in {ref}`section:theodorsen_implementation`. The inputs to the system are the angle of attack $\alpha$ and vertical position $h$, corresponding to pitching and plunging motion, resulting in a generalised input $u = [\alpha, \ h]^T$. Just like in {cite}`brunton2013empirical`, it is more practical to treat $\ddot{u}$ as the input to avoid differentiation in the model. Due to this, the states of the system are the integrals of the input, $\dot{u}$ and $u$, as well as a generalised state of the unsteady flow, $\tilde{x}$. In particular, $\tilde{x}$ can represent either the vorticity or velocity field of the flow (a derivation of this kind of representation of unsteady lift model with vorticity-based states can be found in {cite}`peters2008two`). The output of the system is simply the lift coefficient, $y = C_L$.  The resulting linear model has the form:

:::{math}
	\frac{d}{dt} \begin{bmatrix} \tilde{x} \\ u \\ \dot{u} \end{bmatrix} = \begin{bmatrix} A & 0 & 0\\ 0 & 0 & I\\ 0 & 0 & 0 \end{bmatrix} \begin{bmatrix} \tilde{x} \\ u \\ \dot{u} \end{bmatrix} + \begin{bmatrix} B_{\ddot{u}} \\ 0 \\ I \end{bmatrix} \ddot{u}
:::  

:::{math}
	y = \begin{bmatrix} C & C_u & C_{\dot{u}} \end{bmatrix} \begin{bmatrix} \tilde{x} \\ u \\ \dot{u} \end{bmatrix} + C_{\ddot{u}} \ddot{u}
:::  

In this framework, the added-mass and quasi-steady forces are captured by the $C_u$, $C_\dot{u}$ and coefficients $C_\ddot{u}$, while the matrices $A$, $B_{\ddot{u}}$ and $C$ constitute a model of the transient flow effects on the lift. Importantly, while this is a very high-dimensional system, dimensionality reduction can be used to find a balanced reduced order model. In the case of {cite}`brunton_reduced-order_2013`, the eigensystem realisation algorithm {cite}`juang_eigensystem_1985` was used. Moreover, it has been suggested that the quasi-static and added-mass effects can be identified separately from the transient ones, e.g. by studying the impulse response in $\dot{u}$ or $\ddot{u}$ of the full system. It is also noted that those coefficients ensure the correct asymptotic behaviour of the model for both low and high frequencies. After subtracting the effect of the previously computed coefficients the transient part of the model can be identified using balanced methods. In {cite}`brunton_reduced-order_2013`, the reduced order model of the transient flow was of order 7. It is also worth noting that the identification (with 3 different methods shown) was based on rather simple input signals (impulse/step on one of the inputs, ramp followed by a step), and that it is possible to parametrise the model using the position of the pitching axis $a$.

Given the review of the state of the art above, one can conclude that the unsteady lift model developed in this project can improve on the current ones (most notably the ones from {cite}`brunton_reduced-order_2013`) in 2 aspects:

1. Retaining accuracy in a broader range of angles of attack.
2. Reducing the number of states (without sacrificing accuracy).

Up to this point of this report, the focus of the work has been the reduction of states, since considerations on the accuracy of the viscous models can only be reliable once they start being identified and analysed, which will happen at a later phase of the project.

## Applied approach

The general approach currently used in this phase of the project is to try to fit a model using data generated by the code presented in {ref}`section:theodorsen_implementation`. Consistent with previous works in the field, the inputs to the system are defined to be $\ddot h$ and $\ddot \alpha$, by extension meaning that their integrals: $\dot{h}$, $\dot{\alpha}$ and $\alpha$ are available to the identified model as states. The exact choice of states is discussed in later sections. As explained in {ref}`section:identification_algorithms`, the state evolution of the system is being modelled by performing sparse regression over a nonlinear library of functions $\Theta(x, u)$ of the form:

:::{math}
	\dot{X} = \Xi \Theta^T (X, U)
:::  

where:

* $X$ is the vector of $\mathbb R^n$ snapshots of the states of the system, i.e. $x = [x(t_0),\; x(t_1),\;x(t_2),\;\dots,\,x(t^{N-1})]^T$
* $\dot{X}$ is the vector of $\mathbb R^n$ snapshots of the state derivatives of the system, i.e. $\dot{X} = [\dot{x}(t_0),\; \dot{x}(t_1),\;\dot{x}(t_2),\;\dots,\,\dot{x}(t^{N-1})]^T$
* $U$ is the vector of the history of the inputs, i.e. $U = [u(t_0),\; u(t_1),\;u(t_2),\;\dots,\,u(t^{N-1})]^T$
* $\Theta \in \mathbb R^{N,M}$ contains the time series for $M$ possible base functions
* $\Xi \in \mathbb R^M$ contains the coefficients for all of them

In general, the matrix $\Theta$ can contain as columns any function of the states of the systems. In fact, a careful choice of the base fuctions may help with the accuracy and sparsity of the resulting model, as each term will carry more information about the system. A common library are polynomials up to a maximum degree $d$. An example of how such a basis could be constructed (for a state of dimension 2, and disregarding the control inputs for clarity) is:

:::{math}
 \Theta = [1,\;x_1,\;x_2,\;x_1 x_2,\; x_1^2,\;x_2^2\;\dots,\;x_1^d,\;x_2^d]
:::

Even though this choice seems like one of the most natural, it can in fact present some disadvantages. In particular, in general it cannot assured that it is orthogonal. Orthogonal bases are to be preferred, because they allow for a simpler generalisation of the results, since new terms added to the model don't influence the already found parameters. A visual example might be useful to understand the importance of such a base.

Let the vector $v = [1.5,\; 0.8]^T$ be the alias for the function we are trying to retrieve (see {numref}`OrthogonalBasis`). When projecting to the first canonical basis vector $e_1 = [1,\; 0]^T$ the result is just the first component of $v$. When adding the second canonical base vector $e_2 = [0,\; 1]^T$ the result obtained is the vector $v$ as defined before: the addition of a second base vector didn't modify the first component that we had calculated. If instead of $e_2$ we use the as the second base vector $u = [0.5,\; 1]$ the components of $v$ in this base will be $[1.1,\; 0.8]$: the addition of a second non orthogonal vector modified the already computed value for the first component.  

```{figure} images/orthogonal_basis.png
---
height: 300px
name: OrthogonalBasis
---
Representation of the projection of the vector $v$ on two different sets of basis vectors, $\mathcal C = \{e_1,\;e_2\}$ and $\mathcal B = \{e_1,\; u\}$
```

The fact that newly added parameters don't influence with already computed ones id of great importance if we set as an objective the interpretability of the results: if any modification of the library changes all the parameters, it will be difficult, if not impossible, to identify those that represent a real physical phenomenon and those that are just a byproduct of the fitting technique used. The following 2 sections show a derivation and implementation of a multivariate polynomial basis, which is hoped to provide better models in the next steps of the project.