#pragma once
#include <vector>
#include <string>
#include "FiniteElement.h"
#include "Solver.h"

using Node = Point;
class Spline
{
public:
    Spline();
    Spline(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& f, const std::vector<double>& w);
    void initialize(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& f, const std::vector<double>& w);
    void setGridx(std::vector<double> grid);
    void setGridx(double x0, double x1, double step);
    void setGridy(std::vector<double> grid);
    void setGridy(double y0, double y1, double step);

    void setAlpha(double alpha);

    void assembleGrid();

    std::vector<double> getAssembledMatrix();

    std::vector<double> getFirstMatrix();
    std::vector<double> getSecondMatrix();
    std::vector<double> getBVector();

    std::vector<double> getQVector();

    void solveProblem();

    double getSplineSolution(double x, double y);
    
    //format: x y value
    void printSolutionToFile(std::string name, double precision);

    //format: leftBottom.x leftBottom.y rightTop.x rightTop.y
    //v00 ... v0n
    //::: vmm :::
    //vn0 ... vnn
    void printSolutionToFileAsMatrix(std::string name, double precision);

    //format: x y value
    void printPointsToFile(std::string name);

private:
    void calculateMainMatrix();
    void calculateFirstMatrix();
    void calculateSecondMatrix();

    void calculateBVector();
    void calculateQVector();


    std::vector<double> getLocalFirstMatrix(int currentFE);
    std::vector<double> getLocalB(int currentFE);

    std::vector<double> getLocalSecondMatrix(int currentFE);

    void initializePoints();
    std::vector<int> getGlobalPositionMatrix(int currentFE);
    std::vector<int> getGlobalPositionVector(int currentFE);

    std::vector<FiniteElement> fes;
    std::vector<Node> nodes;
    std::vector<Point> unusedPoints;
    std::vector<bool> isPointUsed;

    std::vector<double> x, y, f, w;
    std::vector<double> gridx, gridy;

    std::vector<double> mainMatrix, firstMatrix, secondMatrix, b, q;

    double alpha = 0;

    int numberOfFiniteElementsx = 0, numberOfFiniteElementsy = 0;
    int globalMatrixRank = 0;

    bool isMainMatrixCalculated = 0, isFirstMatrixCalculated = 0,
        isSecondMatrixCalculated = 0, isBVectorCalculated = 0, isQVectorCalculated = 0;

};

