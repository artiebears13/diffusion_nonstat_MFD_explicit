//
// Created by artem on 12/17/22.
//
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>

#define file_name "out.txt"

void Count_XYZ(int n, int* arr, int NX, int NY){

    arr[2] = n / (NY*NX);
    arr[1] = (n - arr[2] * (NY*NX)) / NX;
    arr[0] = n - arr[2] * (NY*NX) - arr[1] * NX;
    return;
};

void write_to_file(double *massiv, int NX, int NY, int NZ){
    try{
        std::ofstream MyFile(file_name);
        int *xyz = new int[3];
        // Write to the file
        MyFile << NX << ' ' << NY << ' ' << NZ <<"\n"; // NX NY NZ

        for(int n = 0; n < NX*NY*NZ; n++){

            Count_XYZ(n, xyz, NX,NY);
            MyFile << massiv[n] << ' ' << xyz[0] << ' ' << xyz[1] << ' ' << xyz[2] <<"\n"; // U x y z
        }

        // Close the file
        MyFile.close();
    }
    catch(int cod_error){std::cout << "error write to " << file_name << "with error code" << cod_error;}
};

