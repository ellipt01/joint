# joint

## Code for magnetic and gravity joint inversion based on the group lasso regularization.

### 1. Description
This code performs the inversion using magnetic and gravity anomaly data jointly, and reveals the magnetic and gravity structure, which is the subsurface distribution of magnetization and density. It is well known that the magnetic and gravity inverse problems are ill-posed problems, and therefore there are a number of models that equivalently recover the observed data, which can lead to arbitrariness in the derived model. Therefore, previous studies have jointly used magnetic and gravity data, and perform joint inversion by applying the structural coupling, which force to derive highly correlated magnetic and gravity structure, to compensate for the weak constraint of each data on the model.

In our code, we use group lasso regularization. The group lasso is a type of sparse regularization method that achieves variable selection based on the groups of model elements. This method can be used when model elements have inherent groups or clusters, and the model elements of each group together take zero or non-zero values by this regularization.

Generally, to derive the subsurface magnetic and gravity structure, the subsurface space is divided into small grid cells, and the magnetization and density are assigned to each cell, and tried to estimate these parameters to recover the observed anomalies. In this code, the magnetization and density of each cell are grouped and a group lasso is applied. Thus, the magnetization and density of each cell together take on a zero or non-zero value. As a result, the magnetization and density structure become similar, and a highly correlated model is expected to be derived.

### 2. group lasso
Now, denote the subsurface magnetization distribution as $\boldsymbol{\beta}$, that of density as $\boldsymbol{\rho}$, and observed magnetic and gravity anomalies as $\mathbf{f}$ and $\mathbf{g}$, respecively.
The objective function tobe minimized is
$$L(\boldsymbol{\beta},\boldsymbol{\rho};\lambda, \alpha)
# joint

## Code for magnetic and gravity joint inversion based on the group lasso regularization.

### 1. Description
This code performs the inversion using magnetic and gravity anomaly data jointly, and reveals the magnetic and gravity structure, which is the subsurface distribution of magnetization and density. It is well known that the magnetic and gravity inverse problems are ill-posed problems, and therefore there are a number of models that equivalently recover the observed data, which can lead to arbitrariness in the derived model. Therefore, previous studies have jointly used magnetic and gravity data, and perform joint inversion by applying the structural coupling, which force to derive highly correlated magnetic and gravity structure, to compensate for the weak constraint of each data on the model.

In our code, we use group lasso regularization. The group lasso is a type of sparse regularization method that achieves variable selection based on the groups of model elements. This method can be used when model elements have inherent groups or clusters, and the model elements of each group together take zero or non-zero values by this regularization.

Generally, to derive the subsurface magnetic and gravity structure, the subsurface space is divided into small grid cells, and the magnetization and density are assigned to each cell, and tried to estimate these parameters to recover the observed anomalies. In this code, the magnetization and density of each cell are grouped and a group lasso is applied. Thus, the magnetization and density of each cell together take on a zero or non-zero value. As a result, the magnetization and density structure become similar, and a highly correlated model is expected to be derived.

### 2. group lasso
Now, denote the subsurface magnetization distribution as $\boldsymbol{\beta}$, that of density as $\boldsymbol{\rho}$, and observed magnetic and gravity anomalies as $\mathbf{f}$ and $\mathbf{g}$, respecively.
The objective function tobe minimized is
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
Thus, the problem treated by this code is a mixed $L_2$ and group lasso regularized inversion for the magnetic and gravity fields, and $\alpha$ controls the mixing ratio of these two regularizations.

## 3. usage of jinv, $L_2$-group lasso regularized inversion program

    USAGE: jinv
           -f <magnetic anomaly filename>
           -g <gravitic anomaly filename>
           -l <log10(lambda)>
           -a <alpha>
    [optional]
           -t <terrain file>
           -p <parameter filename:default is settings.par>
           -x (output kernel matrices)
           -h (show this message)
The output is

    model.data: derived model. format is x(km)  y(km)  z(km)  beta(A/m)  rho(g/cc)
    recover_mag.data: recovered magnetic anomaly
    recover_rho.data: recovered gravity anomaly

The file specified by the -p option is the settings file for the joint inversion, which must contain the following 6 fields:

    1. number of grid cells: nx, ny, nz
    2. area of analysis: xleft, xright, yleft, yright, ztop, zbottom (km)
    3. altitude of observations: zobs (km)
    4. magnetic inclination and declination (degree): inc, dec
    5. penaltyn parameter: mu
    6. inversion parameters: tolerance, num. of maxiter

    ex.
    1. nx, ny, nz:		50, 50, 25
    2. x, y, zrange (km):	-2., 2., -2., 2., 0., -2.
    3. exf_inc, exf_dec, mgz_inc, mgz_dec(deg.):	45., -7., 45., -7.
    4. tol, maxiter:	1.e-5, 100000
    5. mu:			1.0
    6. nu, beta0, rho0:	1.0, 0., 0.

The format of the settings file is

    identifier(1-6). explanations:	VALUES, VALUES,...


The "terrain file" specified by the -t option is the terrain elevation file. The format is x(km)  y(km)  z(elevation, km).




