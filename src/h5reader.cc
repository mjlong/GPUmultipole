#include "stdlib.h"
#include "hdf5.h"

void readh5_(char* filename, float* x, float* y, float* z){

  hid_t       file_id, dataset_id;  /* identifiers */
  herr_t      status;

  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset_id = H5Dopen(file_id, "x", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (float*)x);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "y", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (float*)y);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "z", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (float*)z);
  status = H5Dclose(dataset_id);

  status = H5Fclose(file_id);
}


void readh5_(char* filename, int* cnt){

  hid_t       file_id, dataset_id;  /* identifiers */
  herr_t      status;

  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset_id = H5Dopen(file_id, "batch_cnt", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (int*)cnt);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
}


void readh5_(char* filename, int* gridsize, int* nbat, 
	     int* meshes, double* width, 
	     double* sigt, double* pf, double* pc){

  hid_t       file_id, dataset_id;  /* identifiers */
  herr_t      status;

  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset_id = H5Dopen(file_id, "num_history", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (int*)gridsize);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "num_batch", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (int*)nbat);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "num_cells", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (int*)meshes);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "width", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double*)width);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "sigma", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double*)sigt);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "pf", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double*)pf);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "pc", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double*)pc);
  status = H5Dclose(dataset_id);

  status = H5Fclose(file_id);
}
