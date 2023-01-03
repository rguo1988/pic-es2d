//general head file
#include <iostream>
#include<chrono>

//project head file
#include"plasma.h"

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
