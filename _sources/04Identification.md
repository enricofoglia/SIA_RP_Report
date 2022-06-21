# System identification algorithms

## General considerations

## SINDy

## Dimensionality reduction

## Input signals

For an accurate identification of a dynamic model, the training data should contain input signals with a wide frequency spectrum, to make sure that all of the system dynamics are present in the data. It is also useful to have signals representing a typical physical scenario, to verify whether the identified model behaves in an intuitive way. The following points described input signals which have been identified as useful for the identification of the Theodorsen model (and, by extension, the low Reynolds number models to be developed later).

### Square wave

A square wave signal periodically switches between a high and low state with a given half-period $T$. Its spectrum has significant peaks at the harmonics of its principal frequency. However, it is a quite useful signal because of its simplicity, enabling to observe the system's response to steps and periodicity in the inputs at the same time.

```{figure} images/signal_square_wave.png
---
width: 600px
name: signal_square_wave
---
Square wave signal with half-period $T=1$ and amplitude of 1. The spectrum has very pronounced peaks at frequencies $f_n = 1/2T + i/T$ [Hz].
```

### Chirp

```{figure} images/signal_chirp.png
---
width: 600px
name: signal_chirp
---
Linear chirp signal with initial $\omega_0$ of $10 [rad/s]$, final $\omega_f$ of $0.1 [rad/s]$ and amplitude of 1, spanning over $20 [s]$. Note that the spectrum is flat in the range $f \in (\omega_f/2\pi, \omega_0/2\pi) [Hz]$.
```

### Pseudo-Random Binary Signal (PRBS)

```{figure} images/signal_prbs.png
---
width: 600px
name: signal_prbs
---
Pseudo-random binary signal with a base time step of 0.2. Note that the spectrum begins diminishing around the maximum switching frequency $1/0.2=5 [Hz]$.
```

### White noise

```{figure} images/signal_white.png
---
width: 600px
name: signal_white
---
White noise with a mean of 0 and standard deviation of 1, sampled in $N=4000$ uniform points.
```

```{figure} images/signal_white_averaged.png
---
width: 600px
name: signal_white_averaged
---
White noise with a mean of 0 and standard deviation of 1, sampled every $\Delta t = 0.005 [s]$ and filtered by a moving average of radius $n=5$. Note that The spectrum begins diminishing at the frequency $f = 1/(2 n \Delta t)$
```