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
    const double x_max = 60.0;

    const double vx_min = -100.0;
    const double vx_max = 100.0;

    const double L = fabs(x_max - x_min);
    const double vx_width = fabs(vx_max - vx_min);
    //wave number
    const double k = 1 * 2.0 * M_PI / L;

    //grids number
    static const int grids_num = 100;
    const double grid_width = L / grids_num;

    //simulated steps & dt
    const int maxsteps = 10;
    const int time_ran = 0;
    const double timestep_condition = 0.1;

    //data path
    const string data_path = "./data/";
    const int data_steps = 10;
    const int data_num = maxsteps / data_steps + 1;
};

class Input: public UniversalParameters
{
  public:
    //electron mass
    static constexpr double m_e = 1.0;
    //electron number
    static constexpr double N_e = 1000000;
    //electron temperature
    static constexpr double T_e = 1;

    //electron/ion average number
    const double n_e_aver = N_e / L;

    //calculated parameters
    const double q_e = -1.0 * sqrt(1.0 / n_e_aver);
    const double w_pe = sqrt(1.0 / m_e);//set nq^2=1
    const double w_p = w_pe;
    const double dt = timestep_condition / w_p;
    const double lambda_D = sqrt(T_e / n_e_aver / q_e / q_e);

    //special settings
    const double kappa = 0.5;
    const double B = 0.5;

    //caculating the normalization
    double A = L * 1.0 / GetNormalization();
    double d = B * (k * k * lambda_D * lambda_D * kappa );

    //define species
    vector<Particles> species;

    double max_probability_density = GetElecInitDistrib(0.0, 0.0);
    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double fx = A * pow(1.0 - B * cos(k * x), -kappa - 1);
        double Tx = T_e * pow(fx, -1.0 / (kappa + 1.0));
        double fv = sqrt(1.0 / (2 * M_PI * Tx * kappa)) * tgamma(kappa + 1.5) / tgamma(kappa + 1.0) * pow(1 + pow(v, 2) / (2 * Tx * kappa), -kappa - 1.5);
        return fx * fv;
    }
    double GetBackgroundIonDensity(double x)
    {
        double fx = A * pow(1.0 - B * cos(k * x), -kappa - 1) + d * cos(k * x);
        return fx;
    }
    double GetNormalization()
    {
        double sum = 0.0;
        double dx = 0.05;
        int nx = L / dx;
        for (int i = 0; i < nx; i++)
        {
            sum += pow(1.0 - B * cos(k * i * dx), -kappa - 1);
        }
        return sum * dx;
    }
    void Initialize();
    void PrintSpecialInformation();
};
#endif
