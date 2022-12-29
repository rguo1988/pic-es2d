#ifndef _partition2d_h
#define _partition2d_h

#include <cmath>

class PartitionToGrids
{
  public:
    int x_idx[4];
    int y_idx[4];
    double contrib[4];
    PartitionToGrids(double _dx, double _dy, double x, double y);
};
PartitionToGrids::PartitionToGrids(double dx, double dy, double x, double y)
{
    int j = floor(x / dx);
    int k = floor(y / dy);
    double X_j = j * dx;
    double Y_k = k * dy;
    double X_j1 = (j + 1) * dx;
    double Y_k1 = (k + 1) * dy;
    //0:(j,k);1:(j+1,k);2:(j,k+1);3:(j+1,k+1)
    x_idx[0] = j;
    x_idx[1] = j + 1;
    x_idx[2] = j;
    x_idx[3] = j + 1;
    y_idx[0] = k;
    y_idx[1] = k;
    y_idx[2] = k + 1;
    y_idx[3] = k + 1;
    contrib[0] = (X_j1 - x) * (Y_k1 - y) / dx / dy;
    contrib[1] = (x - X_j) * (Y_k1 - y) / dx / dy;
    contrib[2] = (X_j1 - x) * (y - Y_k) / dx / dy;
    contrib[3] = (x - X_j) * (y - Y_k) / dx / dy;
}
#endif
