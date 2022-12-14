#include "Spline.h"
#include "Solver.h"
#include <iostream>
#include <fstream>
#include <iomanip>
constexpr double size = 32;
int main()
{
    constexpr double eps = 1e-20;
    constexpr double alpha = 0.0000001;
    //std::vector<double> x = { 0.25, 0.75, 1.25, 1.75, 0.25, 0.75, 1.25, 1.75, 0.25, 0.75, 1.25, 1.75, 0.25, 0.75, 1.25, 1.75 };
    //std::vector<double> y = { 0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75, 1.25, 1.25, 1.25, 1.25, 1.75, 1.75, 1.75, 1.75 };
    //std::vector<double> f = { 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2 };
    //std::vector<double> w = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    std::vector<double> x, y, f, w;

    double xs, ys, fs, ws;
    std::ifstream s("femsolution.txt");
    while (s >> xs >> ys >> fs >> ws)
    {
        x.push_back(xs);
        y.push_back(ys);
        f.push_back(fs);
        w.push_back(ws);
    }
    s.close();

    std::vector<double> dervs, dx, dy;
    std::ifstream d("femderivatives.txt");
    while (d >> xs >> ys >> fs)
    {
        dx.push_back(xs);
        dy.push_back(ys);
        dervs.push_back(fs);
    }


    Spline spline(x, y, f, w);
    spline.setGridx(0, 4, 4 / size);
    spline.setGridy(0, 4, 4 / size);
    spline.setAlpha(alpha);
    spline.assembleGrid();
    
    spline.solveProblem();

    spline.printPointsToFile("points.txt");
    spline.printSolutionToFileAsMatrix("solution.txt", 0.01);

    //for (size_t i = 0; i < f.size(); i++)
    //{
    //    std::cout.precision(10);
    //    std::cout << std::left << std::setw(16) << x[i] << " " << std::setw(16) << y[i] << " "
    //        << std::setw(16) << f[i] << " "
    //        << std::setw(16) << spline.getSplineSolution(x[i], y[i]) << "\n";
    //}
    //std::cout << "\n";
    //for (size_t i = 0; i < dervs.size(); i++)
    //{
    //    std::cout.precision(10);
    //    std::cout << std::left << std::setw(16) << dx[i] << " " << std::setw(16) << dy[i] << " "
    //        << std::setw(16) << dervs[i] << " "
    //        << std::setw(16) << spline.getSplineDerivative(dx[i], dy[i]) << "\n";
    //}

    std::ofstream q("pointsdervs.txt");
    std::ofstream wq("pointssols.txt");

    for (size_t i = 0; i < dervs.size(); i++)
    {
        wq.precision(10);
        wq << spline.getSplineSolution(dx[i], dy[i]) << "\n";
    }
    std::cout << "\n";
    for (size_t i = 0; i < dervs.size(); i++)
    {
        q.precision(10);
        q << spline.getSplineDerivative(dx[i], dy[i]) << "\n";
    }
    
}