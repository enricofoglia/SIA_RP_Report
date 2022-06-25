# Introduction

## Scope of the work

The goal of the project here presented is to identify sparse and reliable models for the unsteady lift of an airfoil undergoing a general motion (with a focus on pitching and plunging) in a low Reynolds number flow, using data from data captured from simulations of increasing complexity. 

Classical analytical models already exist, based on unsteady potential flow theory. In particular, Theodorsen derived a frequency domain function that models a phase shift and gain of unsteady lift in response to sinusoidal variations of its vertical position and angle of attack {cite}`theodorsen1979`. Using the same assumptions, Wagner {cite}`dawson2021` computed the time domain lift response to a step in angle of attack, which was generalised by Von Kármán for arbitrary motion {cite}`VonKarman1938`. However, both time-domain approaches involve integration of the complex Theodorsen function or of the distribution of vorticity in the wake, and can be thus impractical to compute in many real time applications. Moreover, their true limit lies in the assumption of potential flow, discarding viscous effects which can be very pronounced in low Reynolds number flows. 

Nonetheless, improvements to the classical models have been developed in recent years, providing linear state-space representations {cite}`brunton2013empirical`, adapting them to wings of finite span {cite}`boutet2018`, viscous flow with low frequecy oscillations {cite}`liu2015`, or even generalised movement under viscous flow linearised at a certain angle of attack {cite}`brunton_reduced-order_2013`. In fact, it has been proved by Limacher that the classical separation in added mass and circulatory forces already identified by Theodorsen and Von Kármán is indeed more general and applies beyond potential flow theory {cite}`Limacher2018`; this fact has already contributed to the low Reynolds number system identification framework used in {cite}`brunton_reduced-order_2013`. The limitation of the models in {cite}`brunton_reduced-order_2013` is, however, their linearity, leading to accuracy only at a certain range of angles of attack, despite using a reduced order linear model of the unsteady flow and its effects on lift. One of the prinicipal goals of this project is to try to improve on these results by identifying nonlinear models of unsteady lift, which would be either accurate in a broader range of conditions, simpler (e.g. containing less states) or more interpretable. In order to do so, the Sparse Identification for Non-linear Dynamics (SINDy) {cite}`brunton2016_SINDy` machine learning algorithm is being used. The attractiveness of this particular algorithm is due to the fact that it has a mature Python implementation {cite}`kaptanoglu2021_PySINDy` and is designed to retrieve sparse models, i.e. with as little terms as possible. Hopefully, this will result in a nonlinear model that is both accurate and interpretable.

The system identification performed in the project will be based on 3 types of simulation data:

1. Data generated using the Theodorsen model approximation found in {cite}`brunton_reduced-order_2013`.
2. Inviscid CFD simulations.
3. Low Reynolds number viscous simulations.

The rationale is to initially develop the identification methodology based on data that easy to generate, and that is known to be correct (because of using a good model from the literature). Also the identification itself is easier in this phase, because the obtained models can be explicitly compared to the reference one. An intermediate step after that is the repetition of the identification procedure on inviscid CFD data, which should yield a very similar model, followed by attempts to identify a model using the viscous CFD data. The latter will be a long, iterative process, with gradually adapted methodolody and analysis of the results.

## Applications

The project here presented is not just a thought exercise in order to go beyond classical analytical models with modern machine learning techniques. The need for precise and interpretable models for airfoils in unsteady flow conditions is felt in diverse fields of the aerospace industry.

The first analytical models by Wagner and Theodorsen can be applied to the prediction of the onset of flutter and the structural response of a wing to a gust of various shapes {cite}`CHIOCCHIA`.

There are many modern developments in small and agile drones, either in a rotorcraft, fixed-wing or bio-inspired flapping-wing configurations. Given their small sizes and resulting low velocities, accurate models of low Reynolds number unsteady lift can help in their design, both in terms of structure and control algorithm synthesis {cite}`brunton_unsteady_2012`.

The forward flight of an helicopter induces on the rotor a cyclic load due to the fact that the blade encounters a flow with a changing velocity due to the summation of the linear and rotating motion. The prediction of such loads, and thus of the lift force generated by the rotor, can only be made considering the time-varying nature of the flow {cite}`johnson2013`. 

Another emerging field where unsteady airfoil aerodynamics comes into play is Boundary Layer Ingestion engines. Such systems are projected to have up to a 19\% efficiency increment over conventional configurations (considering the overall system improvement {cite}`uranga2018`) and are thus studied with great attention. However, due to the presence at the inlet of a varying velocity field, the fan blades undergo time varying loading that can result in flutter or fatigue damage {cite}`bakhle2018`.

The intrinsic unsteadiness of the flow generates acoustic waves that propagate from the airfoil. This noise can, in fact, be of considerable intensity when it come to the biggest rotors, and can be a deal braker in situations where a low acoustic footprint is paramount, like in the case of urban air mobility using drone-like vehicles. 

## Outline of the report

The first part of the work introduces the topic of unsteady lift by describing in detail the classical Theodorsen model and presenting a Python implementation of its approximate linear state-space representation, reproduced from {cite}`brunton2013empirical`. The Wagner function, related to the Theodorsen model, is also briefly introduced.

The second part introduces system identification concepts relevant to the project. Some very relevant works on unsteady lift system identification are briefly discussed, followed by an overview of the key aspects of the currently used identification methodology (which is subject to evolved over time).

Afterwards, two main pieces of results achieved so far are presented. The first one aims at using SINDy to discover a dynamic model outputting the Wagner function, replicating the results found in {cite}`dawson2022improved`. The second part aims at reproducing the linearised Theodorsen model found in {cite}`brunton2013empirical`, using a modified approach involving less states than the reference model. Finally, the current achievements and future considerations are summarised.