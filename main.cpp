//ver2.9
//general head file
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include<stdio.h>
#include<stdlib.h>
#include<limits>
#include<chrono>

//project head file
#include"particles.h"
#include"bfield.h"
#include"diagnose.h"
#include"plasma.h"
#include"input.h"

using namespace std;

int main()
{
    auto start_time = chrono::high_resolution_clock::now();
    //creat plasma
    PlasmaSystem plasma;

    plasma.Run();
    auto stop_time = chrono::high_resolution_clock::now();
    auto using_time = chrono::duration_cast<chrono::seconds>(stop_time - start_time).count();
    cout << "  Using Time: " << using_time << "s" << endl;
    cout << "--------------------------------------------" << endl;
    return 0;
}
