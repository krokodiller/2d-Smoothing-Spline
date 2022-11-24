#include "Spline.h"
#include "Hermit.h"
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <ranges>
#include <fstream>

constexpr double eps = 0.00001;

Spline::Spline()
{

}

Spline::Spline(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& f, const std::vector<double>& w)
{
    if (x.size() == y.size() && y.size() == f.size() && f.size() == w.size())
    {
        this->x = x;
        this->y = y;
        this->f = f;
        this->w = w;
        initializePoints();
    }
    else
    {
        throw std::invalid_argument("vectors are different sizes");
    }
}

void Spline::initialize(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& f, const std::vector<double>& w)
{
    if (x.size() == y.size() && y.size() == f.size() && f.size() == w.size())
    {
        this->x = x;
        this->y = y;
        this->f = f;
        this->w = w;
        initializePoints();
    }
    else
    {
        throw std::invalid_argument("vectors are different sizes");
    }
}

void Spline::setGridx(std::vector<double> grid)
{
    if (grid.empty())
    {
        throw std::invalid_argument("grid is empty");
    }
    if (!std::is_sorted(grid.begin(), grid.end()))
    {
        std::cerr << "grid is not sorted. Sorting...\n";
        std::sort(grid.begin(), grid.end());
    }
    if (std::adjacent_find(grid.begin(), grid.end()) != grid.end())
    {
        std::cerr << "grid values are not uniqe. Removing duplicates...\n";
        auto last = std::unique(grid.begin(), grid.end());
        grid.erase(last, grid.end());
    }
    this->gridx = grid;
    numberOfFiniteElementsx = grid.size() - 1;
}

void Spline::setGridy(std::vector<double> grid)
{
    if (grid.empty())
    {
        throw std::invalid_argument("grid is empty");
    }
    if (!std::is_sorted(grid.begin(), grid.end()))
    {
        std::cerr << "grid is not sorted. Sorting...\n";
        std::sort(grid.begin(), grid.end());
    }
    if (std::adjacent_find(grid.begin(), grid.end()) != grid.end())
    {
        std::cerr << "grid values are not uniqe. Removing duplicates...\n";
        auto last = std::unique(grid.begin(), grid.end());
        grid.erase(last, grid.end());
    }
    this->gridy = grid;
    numberOfFiniteElementsy = grid.size() - 1;
}

void Spline::setGridx(double x0, double x1, double step)
{
    gridx.clear();
    int i = 0;
    for (double x = x0; x <= x1 + eps; x = x0 + step * i)
    {
        gridx.push_back(x);
        i++;
    }
    numberOfFiniteElementsx = gridx.size() - 1;
}

void Spline::setGridy(double y0, double y1, double step)
{
    gridy.clear();
    int i = 0;
    for (double y = y0; y <= y1 + eps; y = y0 + step * i)
    {
        gridy.push_back(y);
        i++;
    }
    numberOfFiniteElementsy = gridy.size() - 1;
}

void Spline::setAlpha(double alpha)
{
    this->alpha = alpha;
}

void Spline::assembleGrid()
{
    for (size_t i = 0; i < gridy.size(); i++)
    {
        for (size_t j = 0; j < gridx.size(); j++)
        {
            nodes.push_back(Node(gridx[j], gridy[i]));
        }
    }

    for (size_t i = 0; i < numberOfFiniteElementsy; i++)
    {
        for (size_t j = 0; j < numberOfFiniteElementsx; j++)
        {
            fes.push_back(FiniteElement(nodes[i * gridx.size() + j], nodes[i * gridx.size() + j + 1], 
                nodes[(i + 1) * gridx.size() + j], nodes[(i + 1) * gridx.size() + j + 1]));
        }
    }

    globalMatrixRank = nodes.size() * 4;
    mainMatrix.resize(globalMatrixRank * globalMatrixRank);
    firstMatrix.resize(globalMatrixRank * globalMatrixRank);
    secondMatrix.resize(globalMatrixRank * globalMatrixRank);
    b.resize(globalMatrixRank);
    q.resize(globalMatrixRank);
}

std::vector<double> Spline::getAssembledMatrix()
{
    if (isMainMatrixCalculated)
    {
        return mainMatrix;
    }
    else
    {
        calculateMainMatrix();
        return mainMatrix;
    }
}

void Spline::calculateMainMatrix()
{
    if (isMainMatrixCalculated)
    {
        return;
    }
    if (!isFirstMatrixCalculated)
    {
        calculateFirstMatrix();
    }
    if (!isSecondMatrixCalculated)
    {
        calculateSecondMatrix();
    }
    //mainMatrix = firstMatrix;
    std::transform(firstMatrix.begin(), firstMatrix.end(), secondMatrix.begin(), mainMatrix.begin(), std::plus<double>());
    isMainMatrixCalculated = true;
}

std::vector<double> Spline::getFirstMatrix()
{
    if (isFirstMatrixCalculated)
    {
        return firstMatrix;
    }
    else
    {
        calculateFirstMatrix();
        return firstMatrix;
    }
}

void Spline::calculateFirstMatrix()
{
    if (isFirstMatrixCalculated)
    {
        return;
    }
    for (size_t i = 0; i < fes.size(); i++)
    {
        std::vector<int> localMatrixMapping = getGlobalPositionMatrix(i);
        std::vector<double> localMatrix = getLocalFirstMatrix(i);
        for (size_t j = 0; j < localMatrix.size(); j++)
        {
            firstMatrix[localMatrixMapping[j]] += localMatrix[j];
        }

        for (size_t q = 0; q < globalMatrixRank; q++)
        {
            for (size_t z = 0; z < globalMatrixRank; z++)
            {
                std::cout.precision(4);
                std::cout.width(6);
                std::cout << firstMatrix[q * globalMatrixRank + z] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    isFirstMatrixCalculated = true;
}

std::vector<double> Spline::getSecondMatrix()
{
    if (isSecondMatrixCalculated)
    {
        return secondMatrix;
    }
    else
    {
        calculateSecondMatrix();
        return secondMatrix;
    }
}

void Spline::calculateSecondMatrix()
{
    if (isSecondMatrixCalculated)
    {
        return;
    }
    for (size_t i = 0; i < fes.size(); i++)
    {
        std::vector<int> localMatrixMapping = getGlobalPositionMatrix(i);
        std::vector<double> localMatrix = getLocalSecondMatrix(i);
        for (size_t j = 0; j < localMatrix.size(); j++)
        {
            secondMatrix[localMatrixMapping[j]] += localMatrix[j];
        }
    }
    isSecondMatrixCalculated = true;
}

std::vector<double> Spline::getBVector()
{
    if (isBVectorCalculated)
    {
        return b;
    }
    else
    {
        calculateBVector();
        return b;
    }
}

void Spline::calculateBVector()
{
    if (isBVectorCalculated)
    {
        return;
    }
    for (size_t i = 0; i < fes.size(); i++)
    {
        std::vector<int> localVectorMapping = getGlobalPositionVector(i);
        std::vector<double> localVector = getLocalB(i);
        for (size_t j = 0; j < localVector.size(); j++)
        {
            b[localVectorMapping[j]] += localVector[j];
        }

        for (size_t q = 0; q < globalMatrixRank; q++)
        {
            std::cout.precision(4);
            std::cout.width(6);
            std::cout << b[q] << "\n";
        }
        std::cout << "\n";

    }
    isBVectorCalculated = true;
}

std::vector<double> Spline::getQVector()
{
    if (isQVectorCalculated)
    {
        return q;
    }
    else
    {
        calculateQVector();
        return q;
    }
}

void Spline::calculateQVector()
{
    if (isQVectorCalculated)
    {
        return;
    }
    if (!isMainMatrixCalculated)
    {
        calculateMainMatrix();
    }
    if (!isBVectorCalculated)
    {
        calculateBVector();
    }

    q = Solver::gauss(mainMatrix, b);
}

void Spline::solveProblem()
{
    calculateQVector();
}

double Spline::getSplineSolution(double x, double y)
{
    int size = 16;
    auto whichFE = [&]()
    {
        for (int i = 0; i < fes.size(); i++)
        {
            if (x >= fes[i].leftBottom().x - eps && x <= fes[i].rightBottom().x + eps &&
                y >= fes[i].leftBottom().y - eps && y <= fes[i].leftTop().y + eps)
            {
                return i;
            }
        }
        return -1;
    };

    int currentFE = whichFE();

    if (currentFE == -1)
    {
        std::cerr << "trying to get spline solution from outside of area. Returning...\n";
        return 0;
    }
    double value = 0;
    double stepx = fes[currentFE].rightBottom().x - fes[currentFE].leftBottom().x,
        stepy = fes[currentFE].leftTop().y - fes[currentFE].leftBottom().y,
        xi = fes[currentFE].rightBottom().x,
        yi = fes[currentFE].leftTop().y;
    std::vector<int> indices = getGlobalPositionVector(currentFE);
    for (size_t i = 0; i < size; i++)
    {
        value += Hermit::her2d(i, x, xi, stepx, y, yi, stepy) * q[indices[i]];
    }
    return value;
}

void Spline::printSolutionToFile(std::string name, double precision)
{
    std::ofstream f(name);

    double leftBorder = fes.front().leftBottom().x,
        rightBorder = fes.back().rightBottom().x,
        bottomBorder = fes.front().leftBottom().y,
        topBorder = fes.back().leftTop().y;

    for (double i = leftBorder; i <= rightBorder + eps; i += precision)
    {
        for (double j = bottomBorder; j <= topBorder + eps; j += precision)
        {
            f << i << " " << j << " " << getSplineSolution(i, j) << "\n";
        }
    }
    f.close();
}

void Spline::printSolutionToFileAsMatrix(std::string name, double precision)
{
    std::ofstream f(name);

    double leftBorder = fes.front().leftBottom().x,
        rightBorder = fes.back().rightBottom().x,
        bottomBorder = fes.front().leftBottom().y,
        topBorder = fes.back().leftTop().y;

    f << leftBorder << " " << bottomBorder << " " << rightBorder << " " << topBorder << "\n";

    for (double i = leftBorder; i <= rightBorder + eps; i += precision)
    {
        for (double j = bottomBorder; j <= topBorder + eps; j += precision)
        {
            f << getSplineSolution(i, j) << " ";
        }
        f << "\n";
    }
    f.close();
}

void Spline::printPointsToFile(std::string name)
{
    std::ofstream f(name);
    for (size_t i = 0; i < this->f.size(); i++)
    {
        f << x[i] << " " << y[i] << " " << this->f[i] << "\n";
    }
    f.close();
}

std::vector<double> Spline::getLocalFirstMatrix(int currentFE)
{
    size_t size = 16;
    std::vector<double> matrix(size * size);
    std::vector<int> indices;
    double stepx = fes[currentFE].rightBottom().x - fes[currentFE].leftBottom().x,
        stepy = fes[currentFE].leftTop().y - fes[currentFE].leftBottom().y,
        xi = fes[currentFE].rightBottom().x,
        yi = fes[currentFE].leftTop().y;
    for (size_t i = 0; i < x.size(); i++)
    {
        bool betweenX = x[i] >= fes[currentFE].leftBottom().x - eps && x[i] <= fes[currentFE].rightBottom().x + eps;
        bool betweenY = y[i] >= fes[currentFE].leftBottom().y - eps && y[i] <= fes[currentFE].leftTop().y + eps;
        if (betweenX && betweenY && isPointUsed[i] == false)
        {
            indices.push_back(i);
            isPointUsed[i] = true;
        }
    }

    /*auto contains = unusedPoints | std::views::filter([this, currentFE](Point p) {
        bool betweenX = p.x >= fes[currentFE].leftBottom().x - eps && p.x <= fes[currentFE].rightBottom().x + eps;
        bool betweenY = p.y >= fes[currentFE].leftBottom().y - eps && p.y <= fes[currentFE].leftTop().y + eps;
        return betweenX && betweenY;
        });

    std::vector<Point> points(std::begin(contains), std::end(contains));

    auto removed = std::remove_if(unusedPoints.begin(), unusedPoints.end(), [points](Point x) {
        for (size_t i = 0; i < points.size(); i++)
        {
            if (x == points[i])
            {
                return true;
            }
        }
        return false;
        });

    unusedPoints.erase(removed, unusedPoints.end());*/

    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            for (size_t k = 0; k < indices.size(); k++)
            {
                matrix[i * size + j] += w[indices[k]] * Hermit::her2d(i, x[indices[k]], xi, stepx, y[indices[k]], yi, stepy) *
                    Hermit::her2d(j, x[indices[k]], xi, stepx, y[indices[k]], yi, stepy);
            }
        }
    }

    if (currentFE == fes.size() - 1)
    {
        isPointUsed.assign(isPointUsed.size(), false);
    }

    return matrix;
}

std::vector<double> Spline::getLocalSecondMatrix(int currentFE)
{
    int size = 16;
    int miniSize = 4;
    std::vector<double> matrix(size * size);

    std::vector<double> Gx = Hermit::densityMatrix1d(1, fes[currentFE].rightBottom().x - fes[currentFE].leftBottom().x);
    std::vector<double> Gy = Hermit::densityMatrix1d(1, fes[currentFE].leftTop().y - fes[currentFE].leftBottom().y);

    std::vector<double> Mx = Hermit::massMatrix1d(1, fes[currentFE].rightBottom().x - fes[currentFE].leftBottom().x);
    std::vector<double> My = Hermit::massMatrix1d(1, fes[currentFE].leftTop().y - fes[currentFE].leftBottom().y);

    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            matrix[i * size + j] = alpha * (Gx[(Hermit::mu(i + 1) - 1) * miniSize + (Hermit::mu(j + 1) - 1)] *
                My[(Hermit::nu(i + 1) - 1) * miniSize + (Hermit::nu(j + 1) - 1)] +
                Mx[(Hermit::mu(i + 1) - 1) * miniSize + (Hermit::mu(j + 1) - 1)] *
                Gy[(Hermit::nu(i + 1) - 1) * miniSize + (Hermit::nu(j + 1) - 1)]);
        }
    }

    return matrix;
}

std::vector<double> Spline::getLocalB(int currentFE)
{
    size_t size = 16;
    std::vector<double> vector(size);
    std::vector<int> indices;
    double stepx = fes[currentFE].rightBottom().x - fes[currentFE].leftBottom().x,
        stepy = fes[currentFE].leftTop().y - fes[currentFE].leftBottom().y,
        xi = fes[currentFE].rightBottom().x,
        yi = fes[currentFE].leftTop().y;
    for (size_t i = 0; i < x.size(); i++)
    {
        bool betweenX = x[i] >= fes[currentFE].leftBottom().x - eps && x[i] <= fes[currentFE].rightBottom().x + eps;
        bool betweenY = y[i] >= fes[currentFE].leftBottom().y - eps && y[i] <= fes[currentFE].leftTop().y + eps;
        if (betweenX && betweenY && isPointUsed[i] == false)
        {
            indices.push_back(i);
            isPointUsed[i] = true;
        }
    }

    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < indices.size(); j++)
        {
            vector[i] += w[indices[j]] * Hermit::her2d(i, x[indices[j]], xi, stepx, y[indices[j]], yi, stepy) * f[indices[j]];
        }
    }

    if (currentFE == fes.size() - 1)
    {
        isPointUsed.assign(isPointUsed.size(), false);
    }

    return vector;
}

void Spline::initializePoints()
{
    isPointUsed.resize(x.size(), false);
    for (size_t i = 0; i < x.size(); i++)
    {
        unusedPoints.push_back(Point(x[i], y[i]));
    }
}

std::vector<int> Spline::getGlobalPositionMatrix(int currentFE)
{
    int size = 16;
    std::vector<int> mappingNodes;
    int firstNode = (currentFE + currentFE / numberOfFiniteElementsx) * 4;
    for (size_t i = 0; i < 4; i++)
    {
        mappingNodes.push_back(firstNode + i);
    }
    int secondNode = (currentFE + currentFE / numberOfFiniteElementsx + 1) * 4;
    for (size_t i = 0; i < 4; i++)
    {
        mappingNodes.push_back(secondNode + i);
    }
    int thirdNode = (currentFE + numberOfFiniteElementsx + currentFE / numberOfFiniteElementsx + 1) * 4;
    for (size_t i = 0; i < 4; i++)
    {
        mappingNodes.push_back(thirdNode + i);
    }
    int forthNode = (currentFE + numberOfFiniteElementsx + currentFE / numberOfFiniteElementsx + 1 + 1) * 4;
    for (size_t i = 0; i < 4; i++)
    {
        mappingNodes.push_back(forthNode + i);
    }
    
    std::vector<int> returnIndices(size * size);

    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            returnIndices[i * size + j] = mappingNodes[i] * globalMatrixRank + mappingNodes[j];
        }
    }

    return returnIndices;

}

std::vector<int> Spline::getGlobalPositionVector(int currentFE)
{
    int size = 16;
    std::vector<int> mappingNodes;
    int firstNode = (currentFE + currentFE / numberOfFiniteElementsx) * 4;
    for (size_t i = 0; i < 4; i++)
    {
        mappingNodes.push_back(firstNode + i);
    }
    int secondNode = (currentFE + currentFE / numberOfFiniteElementsx + 1) * 4;
    for (size_t i = 0; i < 4; i++)
    {
        mappingNodes.push_back(secondNode + i);
    }
    int thirdNode = (currentFE + numberOfFiniteElementsx + currentFE / numberOfFiniteElementsx + 1) * 4;
    for (size_t i = 0; i < 4; i++)
    {
        mappingNodes.push_back(thirdNode + i);
    }
    int forthNode = (currentFE + numberOfFiniteElementsx + currentFE / numberOfFiniteElementsx + 1 + 1) * 4;
    for (size_t i = 0; i < 4; i++)
    {
        mappingNodes.push_back(forthNode + i);
    }

    return mappingNodes;
}