#include<cmath>
#include<string>
#include<iostream>
#include<fstream>
#include<iomanip>
#include "input.h"
#include "particles.h"
#include<gsl/gsl_rng.h>
#include<sys/timeb.h>
using namespace std;

void Input::PrintSpecialInformation()
{
    cout << "This is the code simulated the formation of kappa plasma!" << endl;
    cout << "Uae = " << uae << " Uai = " << uai << endl;
}
double Input::GetElecInitDistrib(double x, double v)
{
    double ue = 1.0 + uae * cos(k * x);
    double f = sqrt( m_e / (2 * M_PI * T_e) ) * exp(-m_e * pow(v - u, 2) / 2 / T_e);
    return ue / L * f;
}
double Input::GetBackgroundIonDensity(double x)
{
    double ui = 1.0 + uai * cos(k * x);
    return ui;
}
void Input::Initialize()
{
    Particles electrons(N_e, q_e, m_e, "electrons");

    if(if_continued)
    {
        ifstream lastdata_x("./lastdatax");
        ifstream lastdata_v("./lastdatav");
        for (int i = 0; i < electrons.num; i++)
        {
            lastdata_x >> electrons.rv[i].x;
            lastdata_v >> electrons.rv[i].vx;
        }
    }
    else
    {
        //设置随机数
        struct timeb time_seed;
        ftime(&time_seed);
        gsl_rng_default_seed = (time_seed.time * 1000 + time_seed.millitm);
        gsl_rng *r;
        r = gsl_rng_alloc(gsl_rng_default);

        //initialize electrons1
        for(int i = 0; i < electrons.num; i++)
        {
            double temp_vx = gsl_rng_uniform(r) * vx_width + vx_min;
            double temp_x = gsl_rng_uniform(r) * L + x_min;

            while(gsl_rng_uniform(r) * max_probability_density > GetElecInitDistrib(temp_x, temp_vx))
            {
                temp_x = gsl_rng_uniform(r) * L + x_min;
                temp_vx = gsl_rng_uniform(r) * vx_width + vx_min;
            }
            electrons.rv[i].x = temp_x;
            electrons.rv[i].vx = temp_vx;
        }
    }
    species.push_back(electrons);
}
