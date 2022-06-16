# Orthogonal basis and Polynomial Chaos

## Mathematical Background

Let $x$ be a random variable with support $\mathcal X$. As a first instance $x$ will be one dimensional, but it will be showed later that this hypothesis doens't limit our treatment of multiple input systems. Let $w: \mathcal X \rightarrow \mathbb R$ be the probability distribution of $x$. If $f(x)$ and $g(x)$ are two functions of $x$ we can define their inner product as:

:::{math}
:label: innerProduct
 \langle f(x), g(x)\rangle_w = \int_{\mathcal X}f(x)g(x)w(x)dx
:::

Thus it is possible to define the notion of "orthogonality" between functions in the same way as we would for vectors: two functions $f$ and $g$ are orthogonal if and only if their inner product is 0, or:

:::{math}
 \langle f(x), g(x)\rangle_w = \int_{\mathcal X}f(x)g(x)w(x)dx = 0
:::

As stated in {cite}`augustin2008polynomial`, if a function $f(x)$ has finite variance (\langle f(x), f(x)\rangle_w < \infty) then we "expand" it in a series of orthogonal polynomials of $x$ as:

:::{math}
	f(x) = \sum_{i=0}^{\infty} \beta_i\Phi_i(x)
:::

where $\beta_i$ are real coefficients and ${\Phi_i(x)}_ {i\in\mathbb N}$ is the basis of orthogonal polynomials of $x$. Such an expansion can be truncated to an order $d$, thus obtaining an approximation of the function $f(x)$ in a set of orthogonal polynomials in a similar fashion as what is commonly done with Fourier series.

