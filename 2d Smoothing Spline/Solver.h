#pragma once
#include <vector>

class Solver
{
public:
    // works with only square matrices
    static std::vector<double> gauss(std::vector<double> matrix, std::vector<double> vector);
};

