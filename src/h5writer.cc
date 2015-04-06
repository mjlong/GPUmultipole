#include "hdf5.h"
/*
Author: Jilang Miao <jlmiao@mit.edu>
Date  : Feb 3, 2014

writeh5_nxm_ writes "num_vec" vectors of length "length" to h5 file "filename"
the length by num_vec data MUST be allocated continuously like 
length X num_vec array or
long vector with size length X num_vec;

allocating "num_vec" vectors with length length will NOT work
 */
void createmptyh5(char *filename){
  hid_t file /*,group*/;
  herr_t status;
  file = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //group = H5Gcreate (file, "", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //status = H5Gclose (group);
  status = H5Fclose (file);
}

void writeh5_nxm_(char *filename, char *dsetname, int *vec1, int *num_vec, int *length){

  hid_t       file,group, dataset;         /* file and dataset handles */
  hid_t       datatype, dataspace;   /* handles */
  hsize_t     dimsf[2];              /* dataset dimensions */
  herr_t      status;                             

  /*
   * Create a new file using H5F_ACC_TRUNC access,
   * default file creation properties, and default file
   * access properties.
   */
  //file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);


  /* Open an existing group of the specified file. */
  group = H5Gopen(file, "/", H5P_DEFAULT);

  /*
   * Describe the size of the array and create the data space for fixed
   * size dataset. 
   */
  dimsf[0] = *num_vec;
  dimsf[1] = *length;
  dataspace = H5Screate_simple(2, dimsf, NULL); 

  /* 
   * Define datatype for the data in the file.
   * We will store little endian INT numbers.
   */
  datatype = H5Tcopy(H5T_NATIVE_INT);
  status = H5Tset_order(datatype, H5T_ORDER_LE);

  /*
   * Create a new dataset within the file using defined dataspace and
   * datatype and default dataset creation properties.
   */
  dataset = H5Dcreate(group, dsetname, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /*
   * Write the data to the dataset using default transfer properties.
   */
  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double*)vec1);

  /*
   * Close/release resources.
   */
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Gclose(group);
  H5Fclose(file);

  return;
}


void writeh5_nxm_(char *filename, char *dsetname, double *vec1, int *num_vec, int *length){

  hid_t       file,group, dataset;         /* file and dataset handles */
  hid_t       datatype, dataspace;   /* handles */
  hsize_t     dimsf[2];              /* dataset dimensions */
  herr_t      status;                             

  /*
   * Create a new file using H5F_ACC_TRUNC access,
   * default file creation properties, and default file
   * access properties.
   */
  //file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);


  /* Open an existing group of the specified file. */
  group = H5Gopen(file, "/", H5P_DEFAULT);

  /*
   * Describe the size of the array and create the data space for fixed
   * size dataset. 
   */
  dimsf[0] = *num_vec;
  dimsf[1] = *length;
  dataspace = H5Screate_simple(2, dimsf, NULL); 

  /* 
   * Define datatype for the data in the file.
   * We will store little endian INT numbers.
   */
  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
  status = H5Tset_order(datatype, H5T_ORDER_LE);

  /*
   * Create a new dataset within the file using defined dataspace and
   * datatype and default dataset creation properties.
   */
  dataset = H5Dcreate(group, dsetname, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /*
   * Write the data to the dataset using default transfer properties.
   */
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double*)vec1);

  /*
   * Close/release resources.
   */
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
  H5Gclose(group);
  H5Fclose(file);

  return;
}
