# joint

## Code for magnetic and gravity joint inversion based on the group lasso regularization.

## 1. overview
This code performs the inversion using magnetic and gravity anomaly data jointly, and reveals the subsurface magnetic and gravity structures, i.e., the subsurface distribution of magnetization and density. It is well known that since the magnetic and gravity inverse problems are ill-posed problems, there are a number of models that equivalently recover the observed data, which may lead to arbitrariness in the derived model. Therefore, previous studies have jointly used magnetic and gravity data, and perform joint inversion by applying the structural coupling, which force to derive highly correlated magnetic and gravity structure, to compensate for the weak constraint of each data on the model.

In our code, we use group lasso regularization. The group lasso is a type of sparse regularization method that achieves variable selection based on the groups of model elements. This method can be used when model elements have inherent groups or clusters, and the model elements of each group together take zero or nonzero values through this regularization.

Generally, to derive the subsurface magnetic and gravity structure, the subsurface space is divided into small grid cells, and the magnetization and density are assigned to each cell, and an attempt is made to estimate these parameters to recover the observed anomalies. In this code, the magnetization and density of each cell are grouped and a group lasso is applied. In this way, the magnetization and density of each cell together take on a zero or non-zero value. As a result, the magnetization and density structure become similar, and a highly correlated model is expected to be derived.

#### group lasso
Now, denote the subsurface magnetization distribution as $\boldsymbol{\beta}$, that of density as $\boldsymbol{\rho}$, and observed magnetic and gravity anomalies as $\mathbf{f}$ and $\mathbf{g}$, respecively.
The objective function to be minimized is
$$\displaystyle
L(\boldsymbol{\beta},\boldsymbol{\rho};\lambda, \alpha)=
\frac{1}{2}\left\|
	\mathbf{f}-\mathbf{K}\cdot\boldsymbol{\beta}
\right\|^2
+\frac{1}{2}\left\|
	\mathbf{g}-\mathbf{G}\cdot\boldsymbol{\rho}
\right\|^2
+\lambda\alpha\left(
	\frac{1}{2}\left\|\boldsymbol{\beta}\right\|^2
	+\frac{1}{2}\left\|\boldsymbol{\rho}\right\|^2
 \right)
+\lambda(1-\alpha)P_{group}(\boldsymbol{\beta},\boldsymbol{\rho}),$$
where $\mathbf{K}$ and $\mathbf{G}$ are the kernel matrix for magnetic and gravity anomaly, respectively.
The function $P_{group}(\boldsymbol{\beta},\boldsymbol{\rho})$ is the following penalty function of group lasso:
$$\displaystyle
P_{group}(\boldsymbol{\beta},\boldsymbol{\rho})=\sum_{j=1}^M\sqrt{\beta_j^2+\rho_j^2}.$$
Thus, the problem treated by this code is a mixed $L_2$ norm and group lasso regularized inversion for the magnetic and gravity anomalies, and $\alpha$ controls the mixing ratio of these two regularizations.

## 2. compilaton
To compile this program, edit "make.config" and run make.

make.config specifies some compilation options. You will need to edit and modify this file according to the C and C++ compiler you are using. In particular, you will need to modify the BLAS_LIB and BLAS_CFLAGS flags, which specify the BLAS library and its compiler options according to the BLAS library you are using.
The sample make.config included in this repository assumes the use of Intel OneApi 2024 (icx and icpx compilers and BLAS from MKL).

## 3. usage of jinv
After running make, make install creates the core inversion program jinv, the $L_2$ norm and group lasso regularized joint inversion program, in the ./bin directory.
    USAGE: jinv
           -f <magnetic anomaly filename>
           -g <gravitic anomaly filename>
           -l <log10(lambda): regularization parameter>
           -a <alpha: mixing ratio of L2 (alpha) and group lasso (1-alpha)>
    [optional]
           -t <terrain filename>
           -p <setting filename:default is settings.par>
           -x (output kernel matrices)
           -h (show this message)
The format of the input magnetic and gravity anomaly data file is

    xobs(km)  yobs(km)  zobs(km)  anomaly(A/m or mgal)

where the delimiter must be Tab (\t).
The terrain file specified by the -t option is the gridded terrain elevation file. The format is

    x(km)  y(km)  z(elevation, km)
If no terrain file is specified, the surface topography of the study area is assumed to be a flat plane. The grid of the terrain must be the same as that of the subsurface space specified in the settings file described below.

If the -x option is specified, the kernel matrices $\mathbf{K}$ and $\mathbf{G}$ are written to the files "K.mat" and "G.mat", respectively, in MatrixMarket format.
https://math.nist.gov/MatrixMarket/


The output is the following files:

    model.data: derived model. format is x(km)  y(km)  z(km)  beta(A/m)  rho(g/cc)
    recover_mag.data: recovered magnetic anomaly
    recover_rho.data: recovered gravity anomaly

The file specified by the -p option is the settings file for the joint inversion, which must contain the following 6 fields:

    1. number of grid cells: nx, ny, nz
    2. area of analysis: west-edge, east-edge, south-edge, north-edge, z-top, z-bottom (km)
    3. altitude of observations: zobs (km)
    4. inclination and declination of geomagnetic field and magnetization vector: exf_inc, exf_dec, mgz_inc, mgz_dec (deg)
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
