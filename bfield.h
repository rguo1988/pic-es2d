//
//
#ifndef _bfield_h
#define _bfield_h

class ConstMagneticField//
{
  public:
    double Bx;
    double By;
    double Bz;
    ConstMagneticField(double bx = 0, double by = 0, double bz = 1);
};
#endif
