#ifndef _input_h
#define _input_h
#include<string>
#include<cmath>
#include"particles.h"

class UniversalParameters
{
  protected:
    //if the simulation continued from last calculation
    const bool if_continued = 0;

    //configuration space from x_min to x_max, from vx_min to vx_max
    const double x_min = 0.0;
    const double x_max = 10.0;

    const double vx_min = -10.0;
    const double vx_max = 10.0;

    const double L = fabs(x_max - x_min);
    const double vx_width = fabs(vx_max - vx_min);

    //grids number
    const int grids_num = 200;
    const double grid_width = L / grids_num;

    //simulated steps & dt
    const int maxsteps = 1000;
    const int time_ran = 0.0;
    const double timestep_condition = 0.1;

    //data path
    const string data_path = "./data/";
};

class Input: public UniversalParameters
{
  public:
    //electron mass
    static constexpr double m_e = 1.0;
    //electron number
    static constexpr double N_e = 1000;
    //electron temperature
    static constexpr double T_e = 0.1;

    //special parameters
    //two stream bulk speed
    const double u1 = 1.6;
    const double u2 = -1.4;

    //electron/ion average number
    const double n_e_aver = N_e / L;

    //calculated parameters
    const double q_e = -1.0 * sqrt(1.0 / n_e_aver);
    const double w_pe = sqrt(1.0 / m_e);//set nq^2=1
    const double w_p = w_pe;
    const double dt = timestep_condition / w_p;
    const double lambda_D = sqrt(T_e / n_e_aver / q_e / q_e);

    //define species
    vector<Particles> species;

    double max_probability_density = GetStreamOneInitDistrib(0.0, u1);
    double GetStreamOneInitDistrib(double x, double v);
    double GetStreamTwoInitDistrib(double x, double v);
    double GetBackgroundIonDensity(double x);
    void Initialize();

};
#endif
