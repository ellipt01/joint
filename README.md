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

We define a Cartesian coordinate system where the $x$-axis is **positive to the East**, the $y$-axis is **positive to the North**, and the $z$-axis is **positive upward**.

Let $\boldsymbol{\beta}$ be the subsurface magnetization distribution, $\boldsymbol{\rho}$ be the density distribution, and $\mathbf{f}$ and $\mathbf{g}$ be the observed magnetic and gravity anomalies, respectively. The objective function to be minimized is:

$$
L(\boldmath{\beta}},\boldsymbol{\rho};\lambda_1, \lambda_2)=
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
```


## 3. Usage
After installation, the executable `jinv` will be available in the `./bin` directory.

### Command-Line Options
```
USAGE: jinv [options]

Required:
  -f <path>    Path to the magnetic anomaly data file.
  -g <path>    Path to the gravity anomaly data file.
  -l <l1:l2>   Regularization parameters as log10(λ₁) and log10(λ₂).
               (e.g., -l -3:-2)
  -a <a:l>     Regularization parameters as α and log10(λ).
               (e.g., -a 0.5:-3)
<br>
Optional:
  -t <path>    Path to the gridded terrain elevation file.
  -s <path>    Path to the settings file (default: settings.par).
  -v           Enable verbose mode for detailed progress output.
  -h           Show this help message.
```

### Regularization Parameters

The regularization parameters $\lambda_1$ and $\lambda_2$ can be set in one of two mutually exclusive ways:

1.  **Directly with the `-l` option:**
    * `-l <log10(λ₁):log10(λ₂)>`

2.  **Using a mixing parameter `α` with the `-a` option:**
    * `-a <α:log10(λ)>`
    * This sets $\lambda_1 = \alpha \cdot \lambda$ and $\lambda_2 = (1-\alpha) \cdot \lambda$.

**Note:** The `-l` and `-a` options cannot be used at the same time.

---
### Input File Formats

#### Anomaly Data (`-f`, `-g`)
The anomaly data files must be **tab-separated** text files with four columns:
> `x_obs(km)   y_obs(km)   z_obs(km)   anomaly(A/m or mgal)`

#### Terrain Data (`-t`)
The optional terrain file defines the surface topography as a grid of elevation points. The format is a text file with three **space or tab-separated** columns:
> `x(km)   y(km)   z_elevation(km)`

If no terrain file is provided, the surface is assumed to be a flat plane at z=0. The grid layout (x and y coordinates) of the terrain file must match the subsurface grid defined in the settings file.

---
### Output Files

The program generates the following output files in the working directory:

* `model.data`: The final inverted model.
    > Format: `x(km)   y(km)   z(km)   magnetization(A/m)   density(g/cc)`
* `recover_mag.data`: The recovered magnetic anomaly data calculated from the final model.
* `recover_grv.data`: The recovered gravity anomaly data calculated from the final model.

---
### Settings File (`-s`)

The settings file configures the model space, geophysical parameters, and solver settings. Lines starting with `#` are ignored as comments.

```
# example settings file
1. nx, ny, nz:		50, 50, 25
2. x, y, zrange (km):	-2., 2., -2., 2., 0., -2.
3. exf_inc, exf_dec, mgz_inc, mgz_dec(deg.):	45., -7., 45., -7.
4. tol, maxiter:	1.e-5, 100000
5. mu:			1.0
6. nu, beta0, rho0:	1.0, 0., 0.
```

| ID | Description | Values | Example |
| :--- | :--- | :--- | :--- |
| **1** | Number of grid cells | `nx, ny, nz` | `1. nx, ny, nz: 50, 50, 25` |
| **2** | Model space area (km) | `x_min(west), x_max(east),`<br>`y_min(south), y_max(north),`<br> `z_top, z_bottom` | `2. x, y, zrange (km): -2., 2., -2., 2., 0., -2.` |
| **3** | Field/Mag directions (deg) | `exf_inc, exf_dec, mgz_inc, mgz_dec` | `3. ...inc, dec(deg.): 45., -7., 45., -7.` |
| **4** | Solver settings | `tolerance, max_iterations` | `4. tol, maxiter: 1.e-5, 100000` |
| **5** | ADMM penalty parameter | `mu (μ)` | `5. mu: 1.0` |
| **6** | Lower bound constraints | `nu (ν), beta_min, rho_min` | `6. nu, beta0, rho0: 1.0, 0., 0.` |

#### Notes on Settings:

* **Model Space (ID 2):**
    * If no terrain file is used, `z_top` and `z_bottom` define the absolute vertical range. In the example, the space is $z \in [-2, 0]$ km.
    * If a terrain file *is* specified, `z_top` and `z_bottom` are interpreted as offsets **relative to the surface elevation**. In the example, the space would be from `surface - 2km` to `surface + 0km`.
* **Lower Bounds (ID 6):**
    * `nu` ($\nu$) is the ADMM penalty parameter for the lower-bound constraint. If `nu` is zero or negative, the constraint is disabled.
* **Solver (ID 4):**
    * The ADMM solver will stop when the residual falls below the `tolerance` or when 

