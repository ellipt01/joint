# joint

## Code for magnetic and gravity joint inversion based on the group lasso regularization.

## 1. overview
This code performs the inversion using magnetic and gravity anomaly data jointly, and reveals the subsurface magnetic and gravity structures, i.e., the subsurface distribution of magnetization and density. It is well known that since the magnetic and gravity inverse problems are ill-posed problems, there are a number of models that equivalently recover the observed data, which may lead to arbitrariness in the derived model. Therefore, previous studies have jointly used magnetic and gravity data, and performed joint inversion by applying structural coupling, which forces to derive highly correlated magnetic and gravity structure, to compensate for the weak constraint of each data on the model.

In our code, we use the group lasso regularization. The group lasso is a type of sparse regularization method that achieves variable selection based on the groups of model elements. This method can be used when the model elements have inherent groups or clusters, and the model elements of each group together take zero or nonzero values through this regularization.

Generally, to derive the subsurface magnetic and gravity structure, the subsurface space is divided into small grid cells, and the magnetization and density are assigned to each cell, and an attempt is made to estimate these parameters to recover the observed anomalies. In this code, the magnetization and density of each cell are grouped and a group lasso is applied. In this way, the magnetization and density of each cell together take on a zero or non-zero value. As a result, the magnetization and density structure become similar, and a highly correlated model is expected to be derived.

#### group lasso
Now, denote the subsurface magnetization distribution as $\boldsymbol{\beta}$, that of density as $\boldsymbol{\rho}$, and observed magnetic and gravity anomalies as $\mathbf{f}$ and $\mathbf{g}$, respecively.
The objective function to be minimized is

$$\displaystyle
L(\boldsymbol{\beta},\boldsymbol{\rho};\lambda_1, \lambda_2)=
\frac{1}{2}\left\|
	\mathbf{f}-\mathbf{K}\cdot\boldsymbol{\beta}
\right\|^2
+\frac{1}{2}\left\|
	\mathbf{g}-\mathbf{G}\cdot\boldsymbol{\rho}
\right\|^2
+\lambda_2\left(
	\frac{1}{2}\left\|\boldsymbol{\beta}\right\|^2
	+\frac{1}{2}\left\|\boldsymbol{\rho}\right\|^2
 \right)
+\lambda_1P_{group}(\boldsymbol{\beta},\boldsymbol{\rho}),$$

where $\mathbf{K}$ and $\mathbf{G}$ are the kernel matrix for magnetic and gravity anomaly, respectively.
The function $P_{group}(\boldsymbol{\beta},\boldsymbol{\rho})$ is the following penalty function of group lasso:

$$\displaystyle
P_{group}(\boldsymbol{\beta},\boldsymbol{\rho})=\sum_{j=1}^M\sqrt{\beta_j^2+\rho_j^2}.$$

Thus, the problem treated by this code is a mixed $L_2$ norm and group lasso regularized inversion for the magnetic and gravity anomalies.

## 2. compilation
To compile this program, edit "make.config" and run make.

make.config specifies some compilation options. You will need to edit and modify this file according to the C and C++ compiler you are using. In particular, you will need to modify the BLAS_LIB and BLAS_CFLAGS flags, which specify the BLAS library and its compiler options according to the BLAS library you are using.
The sample make.config included in this repository assumes the use of Intel OneApi 2024 (icx and icpx compilers and BLAS from MKL).

## 3. usage of jinv
After running "make", "make install" creates the inversion program "jinv", the $L_2$ norm and group lasso regularized joint inversion program, in the ./bin directory.

    USAGE: jinv
           -f <magnetic anomaly filename>
           -g <gravity anomaly filename>
           -l <log10(lambda1):log10(lambda2)>
           -a <alpha:log10(lambda)>
    [optional]
           -t <terrain filename>
           -s <setting filename:default is settings.par>
           -v (verbos mode)
           -h (show this message)

The -l and -a options specify the regularization parameters $\lambda_1$ and $\lambda_2$.

The -l option specifies

-l $\log_{10}(\lambda_1):\log_{10}(\lambda_2)$,

while -a specifies

-a $\alpha:\log_{10}(\lambda)$

where $\lambda$ is the regularization parameter and $\alpha$ is the mixing ratio of the $L_2$ norm and the group lasso penalty, that spesify $\lambda_1$ and $\lambda_2$ as follows:

$\lambda_1=\alpha\cdot\lambda,\quad \lambda_2=(1-\alpha)\cdot\lambda$

These options -l and -a cannot be used together.

The format of the input magnetic and gravity anomaly data file is

    xobs(km)  yobs(km)  zobs(km)  anomaly(A/m or mgal)

where the delimiter must be Tab (\t).
The terrain file specified by the -t option is the gridded terrain elevation file. The format is

    x(km)  y(km)  z(elevation, km)
If no terrain file is specified, the surface topography of the study area is assumed to be a flat plane. The terrain grid must be the same as that of the subsurface space specified in the settings file described below.



The output is the following files:

    model.data: derived model. format is x(km)  y(km)  z(km)  beta(A/m)  rho(g/cc)
    recover_mag.data: recovered magnetic anomaly
    recover_rho.data: recovered gravity anomaly

The file specified by the -p option is the settings file for the joint inversion, which must contain the following 6 fields:

    1. number of grid cells: nx, ny, nz
    2. area of analysis: west-edge, east-edge, south-edge, north-edge, z-top, z-bottom (km)
    3. inclination and declination of geomagnetic field and magnetization vector: exf_inc, exf_dec, mgz_inc, mgz_dec (deg)
    4. converging tolerance and maximum number of iterations: tol, maxiter
    5. penaltyn parameter: mu
    6. lower bounds: penalty parameter nu, magnetization lower (A/m), density lower (g/cc)

The format of the settings file is

    identifier(1-6). explanations:	VALUES, VALUES,...
    
    ex.
    1. nx, ny, nz:		50, 50, 25
    2. x, y, zrange (km):	-2., 2., -2., 2., 0., -2.
    3. exf_inc, exf_dec, mgz_inc, mgz_dec(deg.):	45., -7., 45., -7.
    4. tol, maxiter:	1.e-5, 100000
    5. mu:			1.0
    6. nu, beta0, rho0:	1.0, 0., 0.

Lines starting with # are considered comments.
In the case of the above example, the subsurface space $x\in$ [-2., 2. (km)], $y\in$ [-2., 2. (km)], and $z\in$ [-2., 0. (km)] is divided into nx=50, ny=50, and nz=25 grid cells, and assign the magnetization $\beta_j$ and density $\rho_j$ ($j=1,2,\cdots,$ nx $\times$ ny $\times$ nz) to each. nu ($\nu$) is a penalty parameter for the lower-bound constraint; if nu is zero or negative, the lower-bound constraint is not applied.

The $L_2$ norm-group lasso regularized problem is solved by an iterative method called the Alternating Direction Method of Multipliers (ADMM),
and the convergence torelance is, in the above example, 1.e-5, and maximum number of iteration is 10000.


If terrain file is spesicied, zrange of 2., z-top, z-bottom (km) indicates that, the vertical range of the study area is from (surface (topography) + z-bottom) km to (surface + z-top) km.
In the case of the above example, $z\in$ [surface - 2, surface (km) ].
