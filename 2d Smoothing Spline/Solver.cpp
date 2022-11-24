#include "Solver.h"
#include <iostream>

constexpr double eps = 0.00001;

std::vector<double> Solver::gauss(std::vector<double> matrix, std::vector<double> vector)
{
    if (vector.size() * vector.size() != matrix.size())
    {
        std::cerr << "vector and matrix have different sizes. Returning...\n";
        return vector;
    }

    std::vector<double> matrixCopy = matrix;
    std::vector<double> vectorCopy = vector;

    int size = vector.size();

    for (size_t k = 0; k < size; k++)
    {
        for (size_t i = 0; i < size; i++)
        {
            if (fabs(matrix[k * size + k]) < eps)
            {
                std::cerr << "diagonal element is zero at " << k << ". Returning...\n";
                return vector;
            }
            matrixCopy[k * size + i] = matrixCopy[k * size + i] / matrix[k * size + k];
        }
        vectorCopy[k] = vectorCopy[k] / matrix[k * size + k];

        for (size_t i = k + 1; i < size; i++)
        {
            double K = matrixCopy[i * size + k] / matrixCopy[k * size + k];
            for (size_t j = 0; j < size; j++)
            {
                matrixCopy[i * size + j] = matrixCopy[i * size + j] - matrixCopy[k * size + j] * K;
            }
            vectorCopy[i] = vectorCopy[i] - vectorCopy[k] * K;
        }
        for (size_t i = 0; i < size; i++)
        {
            for (size_t j = 0; j < size; j++)
            {
                matrix[i * size + j] = matrixCopy[i * size + j];
            }
            vector[i] = vectorCopy[i];
        }
    }

    for (int k = size - 1; k >= 0; k--)
    {
        for (int i = size - 1; i >= 0; i--)
        {
            matrixCopy[k * size + i] = matrixCopy[k * size + i] / matrix[k * size + k];
        }
        vectorCopy[k] = vectorCopy[k] / matrix[k * size + k];

        for (int i = k - 1; i >= 0; i--)
        {
            double K = matrixCopy[i * size + k] / matrixCopy[k * size + k];
            for (int j = size - 1; j >= 0; j--)
            {
                matrixCopy[i * size + j] = matrixCopy[i * size + j] - matrixCopy[k * size + j] * K;
            }
            vectorCopy[i] = vectorCopy[i] - vectorCopy[k] * K;
        }
    }

    return vectorCopy;
}