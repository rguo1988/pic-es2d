#ifndef _plasma_h
#define _plasma_h

#include<gsl/gsl_rng.h>
#include<sys/timeb.h>
#include<eigen3/Eigen/Core>

#include"particles.h"
#include"poisson_solver_2d.h"
#include"bfield.h"
#include"input.h"

using namespace std;
using namespace Eigen;

class PlasmaSystem: public Input
{
  public:
    PoissonSolver2D_DirichletBC poisson_solver;
    ConstMagneticField B;
    MatrixXd charge;

    vector<double> Ek;
    vector<double> Ep;
    vector<double> Et;
    vector<double> T;
    vector<double> num_in_grid;

    PlasmaSystem();

    //print important parameters in simulation
    void PrintParameters() const;

    //trace energy during simulation
    void CalculateE();
    //calculate the temperature in every cell
    void CalculateT();

    //Setup Charge
    void SetupSpeciesChargeOnGrids();
    void SetupBackgroundChargeOnGrids();

    //BorisPusher
    void PushOneStep(int if_init);
    void Run();

};
#endif
