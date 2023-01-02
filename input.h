#ifndef _input_h
#define _input_h
#include<cmath>
#include"particles.h"
#include<string>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<gsl/gsl_rng.h>
#include<sys/timeb.h>
using namespace std;

class Input
{
  public:
    //title
    const string title = "Landau damping of 2D Langmuir waves";
    //if continue from data
    const bool if_continued = 0;

    //simulation box
    static constexpr double kx = 0.5;
    static constexpr double ky = 0.5;
    static constexpr double Lx = 2.0 * M_PI / kx;
    static constexpr double Ly = 2.0 * M_PI / ky;
    const double v_max = 5.0;
    const double vx_width = 2.0 * v_max;
    const double vy_width = 2.0 * v_max;

    static const int nx = 96;
    static const int ny = 96;
    static const int nx_grids = nx - 1;
    static const int ny_grids = ny - 1;
    const double dx = Lx / nx_grids;
    const double dy = Ly / ny_grids;

    //speices
    static constexpr double m_e = 1.0;
    static constexpr double NePerCell = 200;
    static constexpr double T_e = 1;
    const double N_e = NePerCell * nx_grids * ny_grids;
    const double n_e_aver = N_e / Lx / Ly;
    const double q_e = -sqrt(1.0 / n_e_aver);
    const double w_pe = sqrt(1.0 / m_e);//set nq^2=1
    const double lambda_De = sqrt(T_e);

    const double w_p = sqrt(w_pe * w_pe);
    const double lambda_D = sqrt(T_e);

    //time parameters
    const int maxsteps = 200;
    const int time_ran = 0;
    const double timestep_condition = 0.1;
    //const double dt = timestep_condition / w_p;
    const double dt = 0.1;

    //special settings
    static constexpr double d = 0.1;

    //data path
    const string data_path = "./data/";
    const int data_steps = maxsteps;
    const int data_num = maxsteps / data_steps + 1;

    vector<Particles> species;
    static double GetElecInitDistrib(double x, double y, double vx, double vy)
    {
        double ue = 1.0 + d * cos(kx * x);
        double f = sqrt( m_e / (2 * M_PI * T_e) ) * exp(-0.5 * m_e * (vx * vx + vy * vy) / T_e);
        return ue / Lx / Ly * f;
    }
    static double GetElecXDistrib(double x, double y)
    {
        double ue = 1.0 + d * cos(kx * x);
        return ue;
    }
    static double GetElecVDistrib(double vx, double vy)
    {
        double f = sqrt( m_e / (2 * M_PI * T_e) ) * exp(-0.5 * m_e * (vx * vx + vy * vy) / T_e);
        return f;
    }

    //electrons as background
    double GetBackgroundDensity(double x, double y)
    {
        return 1.0;
    }

    void Initialize()
    {
        Particles electrons(N_e, q_e, m_e, "electrons");
        electrons.InitializeXV_Random(GetElecInitDistrib, v_max, Lx, Ly);
        species.push_back(electrons);
    }
    void PrintSpecialInformation()
    {
        cout << "  " << title << endl;
        cout << "  " << "d = " << d << endl;
        cout << "--------------------------------------------" << endl;
    }
};

#endif
