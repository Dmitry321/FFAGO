#ifndef MATRIXMULTY_H
#define MATRIXMULTY_H
#include <functional>
#include <QVector>
//#include <QtConcurrent>
#include <QThreadPool>
#include <QReadWriteLock>

class MatrixMulty
{
public:
    MatrixMulty();

    template <typename T, typename U>
    auto innerProduct(QVector<T> const &vector1, QVector<U> const &vector2) -> decltype(T() * U())
    {
        Q_ASSERT_X(vector1.size() == vector2.size(), "template <typename T, typename U> auto innerProduct(QVector<T> const &vector1, QVector<U> const &vector2) -> decltype(T() * U())", qPrintable("incompatible sizes of vectors, size of vector1 = " + QString::number(vector1.size()) + ", size of vector2 = "+ QString::number(vector1.size())));

        double sum(0.0);

        for (int i = 0; i < vector1.size(); ++i)
            sum += vector1.at(i) * vector2.at(i);

        return sum;
    }

//    template <typename T, typename U>
//    void _multiplicationAuxiliaryFunction(QVector<int> const &indexesToProcess, QVector<QVector<T> > const &matrix1, QVector<QVector<U> > const &matrix2, QVector<QVector<decltype(T() * U())> > &resultMatrix)
//    {
//        for (int i = 0; i < indexesToProcess.size(); ++i) {

//            int currentIndex = indexesToProcess.at(i);

//            QVector<U> currentColumnIndexOfMatrix2 = MatrixOperations::columnVector(currentIndex, matrix2);

//            QVector<decltype(T() * U())> &currentResultVector(resultMatrix[currentIndex]);

//            for (int j = 0; j < matrix1.size(); ++j)
//            currentResultVector[j] = innerProduct(matrix1.at(j), currentColumnIndexOfMatrix2);
//        }
//    }

    template <typename T, typename U>
    void innerSumm(QVector<T> const &vector1, QVector<U> const &vector2, QVector<decltype(T() + U())> &resultMatrix) //-> decltype(T() + U())
    {
        Q_ASSERT_X(vector1.size() == vector2.size(), "template <typename T, typename U> auto innerProduct(QVector<T> const &vector1, QVector<U> const &vector2) -> decltype(T() + U())",
                   qPrintable("incompatible sizes of vectors, size of vector1 = " + QString::number(vector1.size()) + ", size of vector2 = "+ QString::number(vector1.size())));

        resultMatrix = vector1;

        for (int i = 0; i < vector1.size(); ++i)
            resultMatrix[i]+=  vector2.at(i);


    }

    template <typename T, typename U>
    void constMulti(QVector<T> const &vector1, U const &vector2, QVector<decltype(T() * U())> &resultMatrix) //-> decltype(T() + U())
    {
        resultMatrix = vector1;

        for (T &val: resultMatrix) {
             val = vector2*val;
        }
//        for (int i = 0; i < vector1.size(); ++i)
//            resultMatrix[i]*=  vector2;
    }


};

#endif // MATRIXMULTY_H
