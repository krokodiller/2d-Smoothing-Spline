#pragma once
#include <vector>
class Hermit
{
public:
    static double her1d(int n, double x, double xi, double step);
    static double her2d(int n, double x, double xi, double stepx, double y, double yi, double stepy);
    static std::vector<double> densityMatrix1d(double lambda, double h);
    static std::vector<double> massMatrix1d(double gamma, double h);
    
    static int mu(int i);
    static int nu(int i);

};

