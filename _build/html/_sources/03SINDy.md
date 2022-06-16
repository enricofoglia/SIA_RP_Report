# System Identification Algorithms

## Generalities

The general approach that will has been used in this work is to try to fit a model from data produced by the codes presented in <span style="color: red;">section</span> and <span style="color: red;">section</span>. In general we assume the inputs to the system to be $\ddot h$ and $\ddot \alpha$, which means that by analytical integration we have access to all the physical states present in the classical definition of the unsteady lift coefficients for sinusoidal inputs presented in <span style="color: red;">section</span>. The input-output response of the system will be then approximated by a regression model over a library of functions $\Theta(\alpha, \dot \alpha, \ddot \alpha, \dot h, \ddot h)$, in the form:

:::{math}
	C_L = \Theta\xi
:::  
