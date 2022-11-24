#include "Spline.h"
#include "Solver.h"
#include <iostream>
int main()
{
    constexpr double eps = 0.00001;
    constexpr int alpha = 1;
    std::vector<double> x = { 0.25, 0.75, 1.25, 1.75, 0.25, 0.75, 1.25, 1.75, 0.25, 0.75, 1.25, 1.75, 0.25, 0.75, 1.25, 1.75 };
    std::vector<double> y = { 0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75, 1.25, 1.25, 1.25, 1.25, 1.75, 1.75, 1.75, 1.75 };
    std::vector<double> f = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    std::vector<double> w = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    //std::vector<double> x = { 0.5, 0.5, 1.5, 1.5};
    //std::vector<double> y = { 0.5, 1.5, 0.5, 1.5 };
    //std::vector<double> f = { 2, 1, 1, 2 };
    //std::vector<double> w = { 1, 1, 1, 1 };
    Spline spline(x, y, f, w);
    spline.setGridx({ 0, 1, 2 });
    spline.setGridy(0, 2, 2);
    spline.setAlpha(alpha);
    spline.assembleGrid();
    
    spline.solveProblem();

    spline.printPointsToFile("points.txt");
    spline.printSolutionToFileAsMatrix("solution.txt", 0.01);
    

    for (size_t i = 0; i < x.size(); i++)
    {
        std::cout << x[i] << " " << y[i] << " " << spline.getSplineSolution(x[i], y[i]) << "\n";
    }
}