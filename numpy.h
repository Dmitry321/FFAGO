#ifndef NUMPY_H
#define NUMPY_H

#include "eigen/Eigen/Dense"

using namespace Eigen;


class numpy
{
public:
    numpy();

    Eigen::MatrixXf reshape(
        Eigen::MatrixXf &x,
        uint32_t r,
        uint32_t c
        );

    Eigen::MatrixXd reshape(
        Eigen::MatrixXd &x,
        uint32_t r,
        uint32_t c
        );


};

#endif // NUMPY_H
