//测试，调试，以及输出文件
#ifndef _diagnose_h
#define _diagnose_h
#include<fstream>
#include<string>
#include<vector>
#include"particles.h"

using namespace std;

void OutputData(string filename, vector<double> a);
vector<double> GetParticlesVX(Particles testp);
vector<double> GetParticlesVY(Particles testp);
vector<double> GetParticlesX(Particles testp);
vector<double> GetParticlesY(Particles testp);

#endif
