//ver2.9
//general head file
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include<stdio.h>
#include<stdlib.h>
#include<limits>
#include<omp.h>

//project head file
#include"particles.h"
#include"bfield.h"
#include"diagnose.h"
#include"plasma.h"
#include"input.h"

using namespace std;

int main()
{
    double start_time = omp_get_wtime();
    //creat plasma
    PlasmaSystem plasma;

    plasma.Run();
    double stop_time = omp_get_wtime();
    cout << "  Using Time: " << stop_time - start_time << endl;
    return 0;
}
