# System identification algorithms

## General considerations

The central idea of system identification is creating a model of a dynamic system by observing its behaviour. Such a model can have many engineering applications, principally the simulation of the system, controller synthesis and gaining insight into its dynamics. These needs also appear in other sciences — dynamical models of various phenomena, both physical, biological and social, can be used to predict and understand how they change over time.

### Black box methods

For some practical application, like controller synthesis, the accuracy and robustness of the model is usually more important than interpretability. The internal dynamics do not really matter; the system is treated like a "black box" accepting certain inputs and producing corresponding outputs. These problems give rise to black box system identification methods, which strive to reproduce the I/O characteristics of the system without any consideration to the internal structure of the model or its interpretability. This limits generalisability, but is very efficient for single-case applications. There are a few popular black box identification methods (usually formulated for linear discrete-time systems), such as subspace-based methods {cite}`mckelvey_subspace-based_1996` used in the MATLAB System Identification Toolbox. However, because the goal of the presented project is to find interpretable models, black box system identification is not of significant interest.

### Interpretable model discovery

For more general and scientific purposes, the interpretability of the identified models is of higher importance. This tends to lead to a focus on identifying the underlying physical dynamic processes. In these cases, the general form of the model (or at least of its elements) can be derived from first principles, e.g. Newton's laws of mechanics for mechanical systems. The models tend to consist of (or at least contain) terms of known form, but with uncertain parameters. This class of problems leads to the domain of grey box system identification, where prior insight into the dynamics of the system is leveraged to identify more interpretable (and often better) models. If some of the terms of the model can't be established theoretically a priori, data-driven methods can help to identify them. One of such methods is symbolic regression {cite}`schmidt_distilling_2009` {cite}`bongard_automated_2007`, which fits symbolic equations to observation data using genetic programming methods. However, symbolic regression tends to be computationally expensive and does not scale well for larger systems {cite}`brunton2016_SINDy`. Our project utilises a new, simpler method of nonlinear model discovery, called Sparse Identification of Nonlinear Dynamics (SINDy) {cite}`brunton2016_SINDy`, described briefly in the following subsection.

### Reduced order modeling

Another challenge in data-driven model discovery is the high dimensionality of the problems involved. In a fluid flow problem, for example, the state of the system, e.g. defined as the velocity field for potential flow, is infite, since the field is continuous in space. Even if it is discretised, the state can consist of thousands or millions of parameters. The modelling of a system this big is possible (since this is exactly what CFD simulations do), but also very costly in terms of time and computational power. However, even in high-dimensional physical systems, there are usually global, low order patterns playing a significant role in their dynamics. For instance, a low Reynolds number flow over a cylinder or airfoil is significantly marked by a periodic von Karman vortex street. To approximate the global pattern of the flow, it is not necessary to exactly know the state of the velocity field —  simply knowing the current positions of the vortices can be enough. This procedure of finding low-order strucutres that dominate (e.g. carry a significant part of the system's mechanical energy) the dynamics of the system, and subsequent creation of lower-order methods approximating it, is called reduced order modeling (ROM) {cite}`brunton2019data`. One of the main ROM algorithms is the proper orthogonal decomposition, which uses a singular value decomposition to find global patterns in the evolution of the system and sort them by their significance, which allows to discard terms (i.e. combinations of states) of little importance to the overall dynamics {cite}`brunton2019data`. Therefore, ROMs allow to approximate the behaviour of complex systems by using a smaller number of states, therefore significantly reducing the computational cost.

### Balanced models

In control applications it is neither possible nor necessary to know and/or control the full state of such a system, even approximated, but rather some local or global property (which can be defined as an output of the system), like pressure at a certain point of interest or the total lift acting on an object. Models designed for such applications should strive not only to reduce the number of states, like in reduced-order modelling. Since the goal of a controller is to control the output, the model it employs should prioritise the ability to affect the output using the available inputs. This amounts to seeking models with states with two prominent properties from control theory:
1. Observability, i.e. having a measureable influence on the input of the system: if a state does not affect the output in any way, it is useless from the control perspective.
2. Controllability, i.e. the ability to be influenced by the input: if a state cannot be controlled, from the control perspective it is not only useless, but even detrimental if it is observable at the same time, since it constitutes a perturbance in the output that cannot be directly countered. 
Reduced order models that use states which balance observability and controllability are called balanced models {cite}`brunton2019data`. The balanced truncation method has been used in {cite}`brunton2013empirical` to find an accurate and balanced linear model of the Theodorsen function, which is being used in this project.

## SINDy

{cite}`brunton2016_SINDy`

## Input signals

For an accurate identification of a dynamic model, the training data should contain input signals with a wide frequency spectrum, to make sure that all of the system dynamics are present in the data. It is also useful to have signals representing a typical physical scenario, to verify whether the identified model behaves in an intuitive way. The following points describe input signals which have been found useful for the identification of the Theodorsen model (and, by extension, the low Reynolds number models to be developed later).

### Square wave

A square wave signal periodically switches between a high and low state with a given half-period $T$. Its spectrum has significant peaks at the harmonics of its principal frequency, given by the equation $f_n = 1/2T + n/T [Hz], n=0,1,2,...$. However, it is a quite useful signal because of its simplicity, enabling to observe the system's response to steps and periodicity in the inputs at the same time.

```{figure} images/signal_square_wave.png
---
width: 600px
name: signal_square_wave
---
Square wave signal with half-period $T=1$ and amplitude of 1. The spectrum has very pronounced peaks at frequencies $f_n = 1/2T + n/T [Hz], n=0,1,2,...$.
```

### Chirp

A linear chirp is a sine wave whose frequency changes linearly with time from an initial value $\omega_0$ to a final $\omega_f$ over a time $T$:

:::{math}
	\omega(t) = \frac{\omega_f - \omega_0}{T}
:::

:::{math}
	u(t) = A \cdot \sin{(\omega(t) + \phi)}
:::

This structures ensures that every frequency in the range $f \in (\omega_f/2\pi, \omega_0/2\pi) [Hz]$ is excited in a uniform way, as shown in {numref}`signal_chirp`. Another advantage of the signal is its continuity and smoothness.

```{figure} images/signal_chirp.png
---
width: 600px
name: signal_chirp
---
Linear chirp signal with initial $\omega_0$ of $10 [rad/s]$, final $\omega_f$ of $0.1 [rad/s]$ and amplitude of 1, spanning over $20 [s]$. Note that the spectrum is flat in the range $f \in (\omega_f/2\pi, \omega_0/2\pi) [Hz]$.
```

### Pseudo-Random Binary Signal (PRBS)

A pseudo-random binary signal (PRBS) works like a square wave, swithing between a high and low state, except that the switching is random in nature. In the present implementation, after every base time step $\Delta t$ (higher or equal to the sampling step), the signal switches states with a probability of 0.5 based on a pseudo-random binary sequence generated by a linear-feedback shift register (LFSR). The randomness of the swithching behaviour ensures a broad spectrum without significant peaks, up to the frequency related to the base time step $\Delta t$, when the spectrum starts to diminish. The singal and its spectrum is visualised in {numref}`signal_prbs`

```{figure} images/signal_prbs.png
---
width: 600px
name: signal_prbs
---
Pseudo-random binary signal with a base time step of $\Delta t = 0.2$. Note that the spectrum begins diminishing around the maximum switching frequency $f=1/\Delta t [Hz]$. The signal can have better spectral properties for lower base steps $\Delta t$, but this would make this visualisation less clear.
```

### White noise

White noise, by definition, has a flat spectrum, and therefore constitutes a very useful signal for system identification. For the project, it has been generated by sampling a normal distribution of a given mean and standard deviation $\sigma$. 

```{figure} images/signal_white.png
---
width: 600px
name: signal_white
---
White noise with a mean of 0 and standard deviation of 1, sampled in $N=4000$ uniform points over a timespan of $20 [s]$.
```

Since the only limit to the frequency content is the sampling time, the signal can exhibit very erratic behaviour. If this poses a problem to the identification method or simulation model, the frequency spectrum can be reduced using a moving average filtration. Given an integer averaging radius $n$, every point of the signal is transformed into the average of itself and its $n$ neighbours from both sides, therefore limiting the frequency content above the frequency $f = 1/(2 n \Delta t)$.

```{figure} images/signal_white_averaged.png
---
width: 600px
name: signal_white_averaged
---
White noise with a mean of 0 and standard deviation of 1, sampled every $\Delta t = 0.005 [s]$ and filtered by a moving average of radius $n=5$. Note that The spectrum begins diminishing at the frequency $f = 1/(2 n \Delta t)$.
```