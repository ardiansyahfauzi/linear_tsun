# linear_tsun

This is a linear tsunami propagation model in spherical coordinate system with Perfectly Matched Layer (PML) boundary condition.
This model is an extension of original tsunami propagation model in cartesian coordinate system written by Maeda (2016).
Parallelization (OpenMP) is also applied in this model.

# References
Maeda, T., H. Tsushima, and T. Furumura (2016) An effective absorbing boundary condition for linear long-wave and linear dispersive wave tsunami simulations, Earth, Planets, and Space, 68, 63

# How to use
Compile main_llw.f90 for linear long wave or main_ldw.f90 for linear dispersive wave model.
It will include:
- sub_param.f90      (parameters for computation)
- sub_bathymetry.f90 (load bathymetry/topography file)
- sub_initheight.f90 (load initial water surface file)
- sub_station.f90    (load wave gauges file)

If parallel computation is needed, add parallelization command such as "OMP_NUM_THREADS".

.$ OMP_NUM_THREADS=8
.$ gfortran -o llw -fopenmp main_llw.f90
       
Example will be uploaded in the near future.
