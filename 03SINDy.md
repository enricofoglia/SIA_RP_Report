# System Identification Algorithms

## SI: context and state of the art



## Our approach

The general approach that will has been used in this work is to try to fit a model from data produced by the codes presented in <span style="color: red;">section</span> and <span style="color: red;">section</span>. In general we assume the inputs to the system to be $\ddot h$ and $\ddot \alpha$, which means that by analytical integration we have access to all the physical states present in the classical definition of the unsteady lift coefficients for sinusoidal inputs presented in <span style="color: red;">section</span>. The input-output response of the system will be then approximated by a regression model over a library of functions $\Theta(\alpha, \dot \alpha, \ddot \alpha, \dot h, \ddot h)$, in the form:

:::{math}
	C_L = \Theta\xi
:::  

Where $C_L$ is the vector of $\mathbb R^n$ of the collected temporal data of the output of the symulations, so that $C_L = [C_L(0),\; C_L(t_1),\;C_L(t_2),\;\dots,\, C_L(t^{n-1})]^T$, the matrix $\Theta \in \mathbb R^{n,M}$ contains the time series for $M$ possible base functions and $\xi\in \mathbb R^M$ contains the coefficients for all of them. In particular we would hope that $\xi$ is sparse, i.e. that most of its entries are 0, since this would allow some more theoretical analysis of the results.

In general the matrix $\Theta$ can contain as columns any function of the states of the systems. In fact a careful choice of the base fuctions on which the lift coefficient will be projected may help with the accuracy and sparsity of the resulting model, as each term will carry more information about the system. A common library consist in polynomials up to a maximum degree $d$. An example of how such a basis could be constructed is:

:::{math}
 \Theta = [1,\;\alpha_e,\;\dot \alpha,\;\ddot \alpha,\; \ddot h,\;\alpha_e^2,\;\alpha_e\ddot\alpha,\;\dots,\;\ddot h ^d]
:::

Even though this choice seems like one of the most natural it can in fact present some disadvantages. In particular, in general we cannot assure that it is othogonal. Orthogonal basis are to be preferred since they allow for a simpler generalisation of the results, since new terms added to the model don't influence already found parameters. A visual example might be useful to understand the importance of such a base.

Let the vector $v = [1.5,\; 0.8]^T$ be the alias for the function we are trying to retrieve (see {numref}`OrthogonalBasis`). When projecting to the first canonical basis vector $e_1 = [1,\; 0]^T$ the result is just the first component of $v$. When adding the second canonical base vector $e_2 = [0,\; 1]^T$ the result obtained is the vector $v$ as defined before: the addition of a second base vector didn't modify the first component that we had calculated. If instead of $e_2$ we use the as the second base vector $u = [0.5,\; 1]$ the components of $v$ in this base will be $[1.1,\; 0.8]$: the addition of a second non orthogonal vector modified the already computed value for the first component.  

```{figure} images/orthogonal_basis.png
---
height: 200px
name: OrthogonalBasis
---
Representation of the projection of the vector $v$ on two different sets of basis vectors, $\mathcal C = \{e_1,\;e_2\}$ and $\mathcal B = \{e_1,\; u\}$
```