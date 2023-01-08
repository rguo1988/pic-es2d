#ifndef _input_h
#define _input_h
#include<cmath>
#include"particles.h"
#include<string>
#include<iostream>
using namespace std;

class Input
{
  public:
    //title
    const string title = "Landau damping of 2D Ion-Acoustic waves (X-direction)";
    //if continue from data
    const bool if_continued = 0;

    //simulation box
    static constexpr double kx = 1.1;
    static constexpr double ky = 1.0;
    static constexpr double Lx = 2.0 * M_PI / kx;
    static constexpr double Ly = 2.0 * M_PI / ky;
    const double v_max = 5.0;
    const double vx_width = 2.0 * v_max;
    const double vy_width = 2.0 * v_max;

    static const int nx = 256;
    static const int ny = 64;
    static const int nx_grids = nx - 1;
    static const int ny_grids = ny - 1;
    const double dx = Lx / nx_grids;
    const double dy = Ly / ny_grids;

    //speices
    static constexpr double m_e = 1.0;
    static constexpr double NePerCell = 500;
    static constexpr double T_e = 1;
    const double N_e = NePerCell * nx_grids * ny_grids;
    const double n_e_aver = N_e / Lx / Ly;
    const double q_e = -sqrt(1.0 / n_e_aver);
    const double w_pe = sqrt(1.0 / m_e);//set nq^2=1
    const double lambda_De = sqrt(T_e);

    static constexpr double m_i = 100.0;
    static constexpr double NiPerCell = 500;
    static constexpr double T_i = 1;
    const double N_i = NiPerCell * nx_grids * ny_grids;
    const double n_i_aver = N_i / Lx / Ly;
    const double q_i = sqrt(1.0 / n_i_aver);
    const double w_pi = sqrt(1.0 / m_i);//set nq^2=1
    const double lambda_Di = sqrt(T_i);

    //const double w_p = sqrt(w_pe * w_pe);
    //const double lambda_D = sqrt(T_e);
    const double w_p = sqrt(w_pe * w_pe + w_pi*w_pi);
    const double lambda_D = 1.0 / sqrt(1.0 / T_e + 1.0 / T_i);

    //time parameters
    const int maxsteps = 600;
    const int time_ran = 0;
    const double timestep_condition = 0.1;
    //const double dt = timestep_condition / w_p;
    const double dt = 0.1;

    //special settings
    static constexpr double d = 0.3;

    //data path
    const string data_path = "./data/";
    const int data_steps = maxsteps;
    const int data_num = maxsteps / data_steps + 1;

    vector<Particles> species;
    static double GetElecInitDistrib(double x, double y, double vx, double vy)
    {
        //double ue = 1.0 + d * cos(kx * x);
        double ue = 1.0;
        double f = sqrt( m_e / (2 * M_PI * T_e) ) * exp(-0.5 * m_e * (vx * vx + vy * vy) / T_e);
        return ue / Lx / Ly * f;
    }
    static double GetIonInitDistrib(double x, double y, double vx, double vy)
    {
        double ui = 1.0 + d * cos(kx * x);
        //double ui = 1.0;
        double f = sqrt( m_i / (2 * M_PI * T_i) ) * exp(-0.5 * m_i * (vx * vx + vy * vy) / T_i);
        return ui / Lx / Ly * f;
    }

    //electrons as background
    double GetBackgroundDensity(double x, double y)
    {
        return 1.0;
    }

    void Initialize()
    {
        Particles electrons(N_e, q_e, m_e, "electrons");
        Particles ions(N_i, q_i, m_i, "ions");
        electrons.InitializeXV_Random(GetElecInitDistrib, v_max, Lx, Ly);
        ions.InitializeXV_Random(GetIonInitDistrib, v_max, Lx, Ly);
        species.push_back(electrons);
        species.push_back(ions);
    }
    void PrintSpecialInformation()
    {
        cout << "  " << title << endl;
        cout << "  " << "d = " << d << endl;
        cout << "--------------------------------------------" << endl;
    }
};

#endif
