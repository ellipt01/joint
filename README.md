# joint

## 1. Description
Code for magnetic and gravity joint inversion based on the group lasso regularization.

This code performs the inversion using magnetic and gravity anomaly data jointly, and reviel the magnetic and gravity structre that is the subsurface distribution of the magnetization and density.

It is well known that the magnetic and gravity inverse problems are ill-posed problems, and therefore there are a number of models that equivalently recover the observed data, which can lead to arbitrariness in the derived model. Therefore, previous studies have jointly used magnetic and gravity data to compensate for the weak constraint of each data on the model, and have attempted to derive a model of subsurface structure that reduces ambiguity.

To this end, we use group lasso regularization. Group lasso is a type of sparse regularization method that achieves variable selection based on the groups of model elements. This method can be used when model elements have inherent groups or clusters. The model elements of each group together take zero or non-zero values by this regularization.
