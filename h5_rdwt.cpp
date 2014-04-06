//#include <complex>
#include "stdio.h"
#include "stdlib.h"
#include "hdf5.h"
#include <iostream>
#define FILE "092238.h5"

typedef struct{
  double complex[2];
}tuple;


int main() {
  tuple *complex;
  hid_t file_id, dataset_id, complex_id, complex_array_id, space_id;  /* identifiers */
  hsize_t dims[2], datasize_1d[1];
  herr_t status;

  int i, j, k, rank, ndims,
    fitorder[1], dset_data[1], length[1];
  int * l_value;
  double *fit;
  double *mpole;
  double **mpoles;

  file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset_id = H5Dopen(file_id, "/isotop/fitorder", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,fitorder);
  printf("fitorder:%2d\n", fitorder[0]);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotop/length", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,length);
  printf("length:%2d\n", length[0]);
  status = H5Dclose(dataset_id);

  l_value = (int*)malloc(length[0]*sizeof(int));
  dataset_id = H5Dopen(file_id, "/isotop/l_value", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,l_value);
  printf("1st l value:%d\n", l_value[0]);
  free(l_value);
  status = H5Dclose(dataset_id);
  
  fit = (double*)malloc(3*2000*(fitorder[0]+1)*sizeof(double));
  dataset_id = H5Dopen(file_id, "/isotop/fit", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fit);
  printf("fit[0][0][0]:%g\n", fit[11]);
  free(fit);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotop/mpdata", H5P_DEFAULT);
  
  space_id = H5Dget_space(dataset_id);
  ndims = H5Sget_simple_extent_dims(space_id, dims, NULL);
  printf("number of dims:%d, length of dim1:%d,of dim2:%d\n",(int)ndims,(int)dims[0],(int)dims[1]);
  //  dims[0] = 4;
  //  dims[1] = length[0];
  rank = 1; 
  datasize_1d[0]=2;
  complex_array_id = H5Tarray_create(H5T_NATIVE_DOUBLE, rank, datasize_1d);
  complex_id = H5Tcreate(H5T_COMPOUND, sizeof(tuple));
  status = H5Tinsert(complex_id, "point", HOFFSET(tuple, complex), complex_array_id);
  complex = (tuple*)malloc(dims[0]*dims[1]*sizeof(tuple));
  status = H5Dread(dataset_id, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, complex);
  printf("%g+%gi\n",complex[4].complex[0],complex[4].complex[1]);
  
  status = H5Fclose(file_id);
}

