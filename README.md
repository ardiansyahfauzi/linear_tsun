# linear_tsun

This is a linear tsunami propagation model in spherical coordinate system with Perfectly Matched Layer (PML) boundary condition.
This model is an extension of original tsunami propagation model in cartesian coordinate system written by Maeda (2016).
Parallel programming is also applied in this model.

# References
Maeda, T., H. Tsushima, and T. Furumura (2016) An effective absorbing boundary condition for linear long-wave and linear dispersive wave tsunami simulations, Earth, Planets, and Space, 68, 63

# How to use
Just compile llw1.f90 for linear long wave or ldw1.f90 for linear dispersive wave model. 
If parallel computation is needed, add parallelization command such as "OMP_NUM_THREADS".
Example will be uploaded in the near future.
