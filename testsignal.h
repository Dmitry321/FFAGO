#ifndef TESTSIGNAL_H
#define TESTSIGNAL_H

#include <QVector>


class testSignal
{
public:
    testSignal();
    void genSinus(QVector<float> &sinVect, int num, float period, float tay);
};

#endif // TESTSIGNAL_H
