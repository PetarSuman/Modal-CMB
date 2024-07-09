#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include "global.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>


int main( int argc, char *argv[] ){
    
    char inifile[MAXLEN];
	strcpy(inifile, argv[1]);
	initilise(inifile);


    // B = \sum \alpha Q




}


void reconstruct_bispectrum(){
    double *** reconstructed = create_3Darray(alpha_points,alpha_points,alpha_points);
    for (int i = 0; i < alpha_points; i++){
        for (int j = 0; j < alpha_points, j++){
            for (int k = 0; k < alpha_points; k++){
                double sum = 0e0;
                for (int n = 0; n < alpha_max; n++){
                    sum += get_eigen(n) * pijk(r,i,j,k);
                }
                reconstructed[i][j][k] = sum;
            }
        }
    }

    hid_t file_id, dataset_id, dataspace_id, plist_id;
    hsize_t dims[3] = {DIM0, DIM1, DIM2};
    int data[DIM0][DIM1][DIM2];
    int i, j, k;


    // Create the data to be written
    for (i = 0; i < DIM0; i++) {
        for (j = 0; j < DIM1; j++) {
            for (k = 0; k < DIM2; k++) {
                data[i][j][k] = i + j + k;
            }
        }
    }

    // Create a new file using the default properties
    file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Create dataspace
    dataspace_id = H5Screate_simple(3, dims, NULL);

    // Create the dataset creation property list
    plist_id = H5Pcreate(H5P_DATASET_CREATE);

    // Create the dataset
    dataset_id = H5Dcreate(file_id, DATASET, H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);

    // Write the data to the dataset
    H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    // Close the dataset
    H5Dclose(dataset_id);

    // Close the dataspace
    H5Sclose(dataspace_id);

    // Close the file
    H5Fclose(file_id);

    printf("3D array is written to HDF5 file successfully!\n");

    return 0;
}