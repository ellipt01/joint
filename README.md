# Joint Inversion

Joint inversion of magnetic and gravity data using group lasso regularization.

---
## 1. Overview

This program performs a joint inversion of magnetic and gravity anomaly data to reveal subsurface structures, specifically the distribution of magnetization and density.

Geophysical inverse problems are typically **ill-posed**, meaning multiple models can fit the observed data equally well. This ambiguity can lead to non-unique solutions. To address this, joint inversion applies a **structural coupling constraint**. This constraint encourages similarities between the magnetic and gravity structures, providing a more robust and geologically plausible solution than inverting each dataset independently.

Our approach uses **group lasso regularization** to enforce this structural coupling. Group lasso is a sparse regularization method that performs variable selection on predefined groups of parameters. It encourages all parameters within a group to be either all zero or all non-zero simultaneously.

In this code, we model the subsurface by dividing it into a grid of small cells. For each cell $j$, we group its magnetization ($\beta_j$) and density ($\rho_j$) parameters. By applying group lasso to these ($\beta_j, \rho_j$) pairs, we enforce structural similarity: either a cell has **both** non-zero magnetization and density, or it has **neither**. This directly couples the two physical models, leading to a more correlated and reliable result.

---
### Objective Function

Let $\boldsymbol{\beta}$ be the subsurface magnetization distribution, $\boldsymbol{\rho}$ be the density distribution, and $\mathbf{f}$ and $\mathbf{g}$ be the observed magnetic and gravity anomalies, respectively. The objective function to be minimized is:

$$
L(\boldsymbol{\beta},\boldsymbol{\rho};\lambda_1, \lambda_2)=
\frac{1}{2}\left\|
	\mathbf{f}-\mathbf{K}\boldsymbol{\beta}
\right\|_2^2
+\frac{1}{2}\left\|
	\mathbf{g}-\mathbf{G}\boldsymbol{\rho}
\right\|_2^2
+\lambda_2\left(
	\frac{1}{2}\left\|\boldsymbol{\beta}\right\|_2^2
	+\frac{1}{2}\left\|\boldsymbol{\rho}\right\|_2^2
 \right)
+\lambda_1P_{\text{group}}(\boldsymbol{\beta},\boldsymbol{\rho})
$$

where $\mathbf{K}$ and $\mathbf{G}$ are the kernel matrices for magnetic and gravity data. The $P_{\text{group}}(\boldsymbol{\beta},\boldsymbol{\rho})$ term is the group lasso penalty:

$$P_{\text{group}}(\boldsymbol{\beta},\boldsymbol{\rho})=\sum_{j=1}^M\sqrt{\beta_j^2+\rho_j^2}$$

Therefore, this code solves a joint inversion problem with a mixed penalty: a standard $L_2$ Tikhonov regularization term ($\lambda_2$) and the group lasso penalty ($\lambda_1$). This problem is solved using an iterative algorithm called the **Alternating Direction Method of Multipliers (ADMM)**.

---
## 2. Compilation

To compile the program, edit the `make.config` file to match your system's environment and then run the `make` command.

You must configure the `BLAS_LIB` and `BLAS_CFLAGS` variables to correctly link against your system's **BLAS library**. The sample `make.config` provided in this repository is configured for the Intel oneAPI toolkit (2024), using the `icx` and `icpx` compilers and the MKL BLAS library.

```bash
# 1. Edit the configuration file
nano make.config

# 2. Compile the source code
make

# 3. Install the executable to the ./bin directory
make install
