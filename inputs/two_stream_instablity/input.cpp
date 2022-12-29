#include<cmath>
#include<string>
#include<iostream>
#include<iomanip>
#include "input.h"
#include "particles.h"
#include<gsl/gsl_rng.h>
#include<sys/timeb.h>
using namespace std;

double Input::GetStreamOneInitDistrib(double x, double v)
{
    //double ue = 1.0 + 0.02 * cos( 2.0 * M_PI * x / L);
    double f1 = 0.5 * sqrt( m_e / (2 * M_PI * T_e) ) * exp(-m_e * pow(v - u1, 2) / 2 / T_e);
    //return ue / L * f1;
    return f1;
}
double Input::GetStreamTwoInitDistrib(double x, double v)
{

    double f2 = 0.5 * sqrt( m_e / (2 * M_PI * T_e) ) * exp(-m_e * pow(v - u2, 2) / 2 / T_e);
    return f2;
}
double Input::GetBackgroundIonDensity(double x)
{
    return 1.0;
}
void Input::Initialize()
{
    Particles electrons1(N_e, q_e, m_e, "e1");
    Particles electrons2(N_e, q_e, m_e, "e2");

    //设置随机数
    struct timeb time_seed;
    ftime(&time_seed);
    gsl_rng_default_seed = (time_seed.time * 1000 + time_seed.millitm);
    gsl_rng *r;
    r = gsl_rng_alloc(gsl_rng_default);

    //initialize electrons1
    for(int i = 0; i < electrons1.num; i++)
    {
        double temp_vx = gsl_rng_uniform(r) * vx_width + vx_min;
        double temp_x = gsl_rng_uniform(r) * L + x_min;

        while(gsl_rng_uniform(r) * max_probability_density > GetStreamOneInitDistrib(temp_x, temp_vx))
        {
            temp_x = gsl_rng_uniform(r) * L + x_min;
            temp_vx = gsl_rng_uniform(r) * vx_width + vx_min;
        }
        electrons1.rv[i].x = temp_x;
        electrons1.rv[i].vx = temp_vx;
    }
    //initialize electrons2
    for(int i = 0; i < electrons2.num; i++)
    {
        double temp_vx = gsl_rng_uniform(r) * vx_width + vx_min;
        double temp_x = gsl_rng_uniform(r) * L + x_min;

        while(gsl_rng_uniform(r) * max_probability_density > GetStreamTwoInitDistrib(temp_x, temp_vx))
        {
            temp_x = gsl_rng_uniform(r) * L + x_min;
            temp_vx = gsl_rng_uniform(r) * vx_width + vx_min;
        }
        electrons2.rv[i].x = temp_x;
        electrons2.rv[i].vx = temp_vx;
    }
    gsl_rng_free(r);

    species.push_back(electrons1);
    species.push_back(electrons2);
}
