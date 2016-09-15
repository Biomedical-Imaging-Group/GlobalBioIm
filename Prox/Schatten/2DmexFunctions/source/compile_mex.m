function compile_mex 

if ispc
  mex -v projectSpMat2x2.cpp CFLAGS="\$CFLAGS /openmp" LDFLAGS="\$LDFLAGS /openmp" -largeArrayDims  -I../headers/
else
  mex -v projectSpMat2x2.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims  -I../headers/
end