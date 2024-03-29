#include "hdf5IO.h"

void h5read(char filename[]) {
}

/*
  //  tuple *complext;
  tuple z;
  hid_t file_id, 
    dataset_id, 
    complex_id, complex_array_id, 
    space_id;  
  hsize_t dims[2], datasize_1d[1];
  herr_t status;

  int i, j, k, ndims,
    iW, cnt,maxwindow=0;
  int ivalue[1];
  CMPTYPE dvalue[1];

  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset_id = H5Dopen(file_id, "/isotope/fissionable", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ivalue);
  pole.fissionable = ivalue[0];
  //printf("fissionable:%2d\n", pole.fissionable);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotope/mode", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ivalue);
  pole.mode = ivalue[0];
  //printf("mode:%2d\n", pole.mode);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotope/windows", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ivalue);
  pole.windows = ivalue[0];
  //printf("windows:%2d\n", pole.windows);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotope/sqrtAWR", H5P_DEFAULT);
#if defined(__CFLOAT)
  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,dvalue);
#else
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,dvalue);
#endif
  pole.sqrtAWR = dvalue[0];
  //printf("sqrtAWR:%g\n", pole.sqrtAWR);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotope/startE", H5P_DEFAULT);
#if defined(__CFLOAT)
  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,dvalue);
#else
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,dvalue);
#endif
  pole.startE = dvalue[0];
  //printf("startE:%g\n", pole.startE);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotope/endE", H5P_DEFAULT);
#if defined(__CFLOAT)
  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,dvalue);
#else
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,dvalue);
#endif
  pole.endE = dvalue[0];
  //printf("endE:%g\n", pole.endE);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotope/spacing", H5P_DEFAULT);
#if defined(__CFLOAT)
  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,dvalue);
#else
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,dvalue);
#endif
  pole.spacing = dvalue[0];
  //printf("spacing:%g\n", pole.spacing);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotope/fitorder", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ivalue);
  pole.fitorder = ivalue[0];
  //printf("fitorder:%2d\n", pole.fitorder);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotope/length", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ivalue);
  pole.length = ivalue[0];
  //printf("length:%2d\n", pole.length);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotope/numL", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,ivalue);
  pole.numL = ivalue[0];
  //printf("numL:%2d\n", pole.numL);
  status = H5Dclose(dataset_id);

  pole.pseudo_rho = (CMPTYPE*)malloc(pole.numL*sizeof(CMPTYPE));
  dataset_id = H5Dopen(file_id, "/isotope/pK0RS", H5P_DEFAULT);
#if defined(__CFLOAT)
  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pole.pseudo_rho);
#else
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pole.pseudo_rho);
#endif
  //printf("pK0RS:%g,%g\n", pole.pseudo_rho[0],pole.pseudo_rho[1]);
  status = H5Dclose(dataset_id);

  pole.l_value = (unsigned*)malloc(pole.length*sizeof(unsigned));
  dataset_id = H5Dopen(file_id, "/isotope/l_value", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,pole.l_value);
  //printf("4th l value:%d\n", pole.l_value[3]);
  status = H5Dclose(dataset_id);
  
  pole.fit = (CMPTYPE*)malloc(3*pole.windows*(pole.fitorder+1)*sizeof(CMPTYPE));
  dataset_id = H5Dopen(file_id, "/isotope/fit", H5P_DEFAULT);
#if defined(__CFLOAT)
  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pole.fit);
#else
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pole.fit);
#endif
  //printf("fit[0][0][0]:%g\n", pole.fit[11]);
  status = H5Dclose(dataset_id);

  pole.w_start = (int*)malloc(pole.windows*sizeof(int));
  dataset_id = H5Dopen(file_id, "/isotope/wstart", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pole.w_start);
  //printf("w_start:%2d\n", pole.w_start[11]);
  status = H5Dclose(dataset_id);

  pole.w_end = (int*)malloc(pole.windows*sizeof(int));
  dataset_id = H5Dopen(file_id, "/isotope/wend", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pole.w_end);
  //printf("w_end:%2d\n", pole.w_end[11]);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/isotope/mpdata", H5P_DEFAULT);
  space_id = H5Dget_space(dataset_id);
  ndims = H5Sget_simple_extent_dims(space_id, dims, NULL);
  //printf("number of dims:%d, length of dim1:%d,of dim2:%d\n",(int)ndims,(int)dims[0],(int)dims[1]);
  datasize_1d[0]=2;
#if defined(__CFLOAT)
  complex_array_id = H5Tarray_create(H5T_NATIVE_FLOAT, 1, datasize_1d);
#else
  complex_array_id = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, datasize_1d);
#endif
  complex_id = H5Tcreate(H5T_COMPOUND, sizeof(CMPTYPE)*2);
  status = H5Tinsert(complex_id, "point", HOFFSET(tuple, complex), complex_array_id);
  pole.mpdata = (CPUComplex<CMPTYPE>*)malloc(dims[0]*dims[1]*sizeof(tuple));
  status = H5Dread(dataset_id, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, pole.mpdata);
  //printf("%g+%gi\n",real(pole.mpdata[4]),imag(pole.mpdata[4]));
  status = H5Fclose(file_id);
}

*/
