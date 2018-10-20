#include "testsignal.h"
#include <QtMath>

testSignal::testSignal()
{

}

//=========================================================================
// генерация тестового импульсного сигнала
//=========================================================================
void testSignal::genSinus(QVector<float>  &y, int num, float period, float tay)
{

    //float sampl = 1000.0/titul.tay; // перевод в секунды в -1 степени

    //period=0.3;
    float sigAmp=10e-45;
    float freq = 2.f*M_PI/period;
    //num = 100000;
    //tay = 0.0001;

    float x;
    //QVector<float> x(num,1);
    y.fill(1,num);
    int phase = 0;
    double ti = 0.;

    for(int j=0; j<num; j++)
    {
        x=j*tay + ti + phase;

       //y = sigAmp*exp(100.0*sin(x*freq-Pi/4.0)) + level; // + noiseAmp*sin(x*freq2);//noiseAmp*qrand();

        y[j] = sigAmp*qExp(100.0f*qSin(x*freq-M_PI/4.f));
    }

}

