% change path to eigen3 as necessary
mex -I. -I/opt/local/include/eigen3 sturm_mex.cpp
mex -I. -I/opt/local/include/eigen3 roots_mex.cpp