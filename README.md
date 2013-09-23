QPM_MOCKS
=========

This is part of the mockFactory software package.

This code does two separate things:
 1. creates mock halo file from local particle densities
 2. populates halos with mock galaxies 

If you use this code in research that results in publications, please cite the following paper: 

    White, M., Tinker, J., McBride, C.K., 2013, MNRAS, submitted.
    "Mock galaxy catalogs using the quick particle mesh method"


GENERAL USAGE
-------------

Options are specified in the parameter file. There is an anotated 
parameter file at: `examples/annotated_qpm.bat.template_0.6452`


DEPENDENCIES
------------

This code requires some numerical recipes (Press et al), which is not free to 
distribute. For private versions, this resides in the nrlib/ subdirectory. 
For public code, please find the approriate code and make sure the Makefile
points to that directory. 

