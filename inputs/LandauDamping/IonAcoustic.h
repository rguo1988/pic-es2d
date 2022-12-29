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
    const string title = "Landau damping of Ion-acoustic waves (QuietStart)";
    //if continue from data
    const bool if_continued = 0;

    //simulation box
    static constexpr double k = 1.0;
    static constexpr double L = 2.0 * M_PI / k;
    const double v_max = 5.0;
    const double vx_width = 2.0 * v_max;

    static const int nx = 201;
    static const int nx_grids = nx - 1;
    const double dx = L / nx_grids;

    //speices
    static constexpr double m_e = 1.0;
    static constexpr double NePerCell = 20000;
    static constexpr double T_e = 1;
    const double N_e = NePerCell * nx_grids;
    const double n_e_aver = N_e / L;
    const double q_e = -sqrt(1.0 / n_e_aver);
    const double w_pe = sqrt(1.0 / m_e);//set nq^2=1
    const double lambda_De = sqrt(T_e);

    static constexpr double m_i = 100.0;
    static constexpr double NiPerCell = 20000;
    static constexpr double T_i = 1;
    const double N_i = NiPerCell * nx_grids;
    const double n_i_aver = N_i / L;
    const double q_i = sqrt(1.0 / n_i_aver);
    const double w_pi = sqrt(1.0 / m_i);//set nq^2=1
    const double lambda_Di = sqrt(T_i);

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

    static double GetElecInitDistrib(double x, double v)
    {
        double ue = 1.0;
        double f = sqrt( m_e / (2 * M_PI * T_e) ) * exp(-0.5 * m_e * v * v / T_e);
        return ue / L * f;
    }
    static double GetElecXDistrib(double x)
    {
        return 1.0;
    }
    static double GetElecVDistrib(double v)
    {
        double f = sqrt( m_e / (2 * M_PI * T_e) ) * exp(-0.5 * m_e * v * v / T_e);
        return f;
    }
    static double GetIonInitDistrib(double x, double v)
    {
        double ui = 1.0 + d * cos(k * x);
        double f = sqrt( m_i / (2 * M_PI * T_i) ) * exp(-0.5 * m_i * v * v / T_i);
        return ui / L * f;
    }
    static double GetIonXDistrib(double x)
    {
        double ui = 1.0 + d * cos(k * x);
        return ui;
    }
    static double GetIonVDistrib(double v)
    {
        double f = sqrt( m_i / (2 * M_PI * T_i) ) * exp(-0.5 * m_i * v * v / T_i);
        return f;
    }

    //electrons as background
    double GetBackgroundDensity(double x)
    {
        return 1.0;
    }

    void Initialize()
    {
        Particles electrons(N_e, q_e, m_e, "electrons");
        Particles ions(N_i, q_i, m_i, "ions");
        //electrons.InitializeXV_Random(GetElecInitDistrib, v_max, L);
        electrons.InitializeXV_Quiet(GetElecXDistrib, GetElecVDistrib, L, dx, nx_grids);
        //ions.InitializeXV_Random(GetIonInitDistrib, v_max, L);
        ions.InitializeXV_Quiet(GetIonXDistrib, GetIonVDistrib, L, dx, nx_grids);
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
