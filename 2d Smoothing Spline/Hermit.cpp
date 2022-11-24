#include "Hermit.h"
#include <iostream>


int Hermit::mu(int i)
{
    return 2 * (((i - 1) / 4) % 2) + ((i - 1) % 2) + 1;
}

int Hermit::nu(int i)
{
    return 2 * ((i - 1) / 8) + ((i - 1) / 2) % 2 + 1;
}

double Hermit::her1d(int n, double x, double xi, double step)
{
    double k = (x - xi) / step;
    switch (n)
    {
    case 0:
        return 1 - 3 * k * k + 2 * k * k * k;
    case 1:
        return step * (k - 2 * k * k + k * k * k);
    case 2:
        return 3 * k * k - 2 * k * k * k;
    case 3:
        return step * (-(k * k) + k * k * k);
    default:
        std::cerr << "incorrect xi n(1d).\n";
        return 0;
    }
}

double Hermit::her2d(int n, double x, double xi, double stepx, double y, double yi, double stepy)
{
    n = n + 1;
    if (n >= 1 && n <= 16)
    {
        return her1d(mu(n) - 1, x, xi, stepx) * her1d(nu(n) - 1, y, yi, stepy);
    }
    else
    {
        std::cerr << "incorrect xi n(2d).\n";
        return 0;
    }
}

std::vector<double> Hermit::densityMatrix1d(double lambda, double h)
{
    std::vector<double> matrix(16);
    matrix[0] = 36;
    matrix[1] = 3 * h;
    matrix[2] = -36;
    matrix[3] = 3 * h;
    matrix[4] = 3 * h;
    matrix[5] = 4 * h * h;
    matrix[6] = -3 * h;
    matrix[7] = -(h * h);
    matrix[8] = -36;
    matrix[9] = -3 * h;
    matrix[10] = 36;
    matrix[11] = -3 * h;
    matrix[12] = 3 * h;
    matrix[13] = -(h * h);
    matrix[14] = -3 * h;
    matrix[15] = 4 * h * h;

    for (size_t i = 0; i < matrix.size(); i++)
    {
        matrix[i] = (matrix[i] * lambda) / (30 * h);
    }
    return matrix;
}

std::vector<double> Hermit::massMatrix1d(double gamma, double h)
{
    std::vector<double> matrix(16);
    matrix[0] = 156;
    matrix[1] = 22 * h;
    matrix[2] = 54;
    matrix[3] = -13 * h;
    matrix[4] = 22 * h;
    matrix[5] = 4 * h * h;
    matrix[6] = 13 * h;
    matrix[7] = -3 * (h * h);
    matrix[8] = 54;
    matrix[9] = 13 * h;
    matrix[10] = 156;
    matrix[11] = -22 * h;
    matrix[12] = -13 * h;
    matrix[13] = -3 * (h * h);
    matrix[14] = -22 * h;
    matrix[15] = 4 * h * h;

    for (size_t i = 0; i < matrix.size(); i++)
    {
        matrix[i] = (matrix[i] * gamma * h) / (420);
    }
    return matrix;
}