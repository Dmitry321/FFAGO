#include "numpy.h"

numpy::numpy()
{

}

MatrixXf numpy::reshape(MatrixXf &x, uint32_t r, uint32_t c)
{
    //Map<MatrixXf> rx(x.data(), c, r);

     MatrixXf rx = Map<MatrixXf>(x.data(), c, r);
    rx.transposeInPlace();
    return rx;
}

MatrixXd numpy::reshape(MatrixXd &x, uint32_t r, uint32_t c)
{
    MatrixXd rx = Map<MatrixXd>(x.data(), c, r);

    rx.transposeInPlace();

    return rx;
}


