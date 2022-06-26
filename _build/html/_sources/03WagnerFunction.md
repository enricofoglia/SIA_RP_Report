(section:wagner)=
# The Wagner function 

The Wagner function can be seen as the time domain representation of the Theodorsen function. In particular, it can be defined as the lift time response of an airfoil undergoing a sudden change in angle of attack for time $t=0$, i.e. the step response of the system to a change in $\alpha_e$ {cite}`dawson2022improved`.

Calling $\phi(t)$ the Wagner function, the analytical computation can be carried out as:
:::{math}
\begin{equation}
        \phi(t) = \mathcal{L}^{-1}\left(\frac{C(s)}{s}\right)
\end{equation}
:::

where $\mathcal{L}^{-1}$ is the inverse Laplace transform, and $C(s)$ is the Theodorsen function. Defining the real and imaginary part of the Theodorsen function as $G_1(k)$ and $G_2(k)$ respectively, it is possible to compute the Wagner function as {cite}`peters2008two`:

:::{math}
:label: Wagner_analytical
\begin{eqnarray}
  \phi(t) &=& \frac{1}{2} + \frac{2}{\pi}\int_0^\infty \frac{1}{k}\left(G_1(k) - \frac{1}{2}\right)\sin(kt)dk \\
  &=& 1 + \frac{2}{\pi}\int_0^\infty \frac{1}{k}G_2(k)\cos(kt)dk\label{eqn:int_large_t} 
\end{eqnarray}
:::

The first equation has better numerical properties when integrated for small time steps, while the second is more stable for large ones {cite}`dawson2022improved`.

The computation of the Wagner function in this fashion is rather cumbersome, as it requires to integrate numerically the Theodorsen function. In particular for large times the trigonometric terms in the above equations make the integrals slowly convergent and require a very fine discretization of the domain, leading to a long computation time. Many analytical approximations to the function exist {cite}`brunton2013empirical`, but usually their accuracy is limited to a certain time scale {cite}`dawson2022improved`.

As a first step in the project, the SINDy algorith presented in {ref}`section:identification_algorithms` has been used to identify a dynamical system generating the Wagner function, in an attempt to replicate the study conducted in {cite}`dawson2022improved` and gain confidence in the numerical methods involved. The results are presented in {ref}`section:results_wagner`.

The Wagner function, even in an approximated form, can be the base of useful unsteady lift model, since the response to arbitrary input can be computed from it using a convolution integral. However, models based on this framework have the 2 main drawbacks of not being easy to use for many control synthesis problems (state-space representation are often more useful in this regard) and not being easy to parametrise, since the entire dynamics of unsteady lift are compressed into a single function {cite}`brunton_reduced-order_2013`.