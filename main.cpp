#include <QCoreApplication>
#include "ffasettings.h"
#include "ffa.h"
#include "matrixmulty.h"
#include <QDebug>
#include "ffa_tools.h"
#include "testsignal.h"
#include "mydisplay.h"
#include "myqvector.h"

//#include "eigen/Eigen/Dense"
#include <iostream>

//using Eigen::MatrixXf;

//using namespace Eigen;
using namespace std;

// for redirect message
void myMessageOutput(QtMsgType type, const QMessageLogContext &context, const QString &msg)
{
    QByteArray localMsg = msg.toLocal8Bit();
    switch (type) {
    case QtDebugMsg:
        fprintf(stderr, "%s\n", localMsg.constData());
        //fprintf(stderr, "Debug: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtInfoMsg:
        fprintf(stderr, "Info: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtWarningMsg:
        fprintf(stderr, "Warning: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtCriticalMsg:
        fprintf(stderr, "Critical: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtFatalMsg:
        fprintf(stderr, "Fatal: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    }
}

int main(int argc, char *argv[])
{
    qInstallMessageHandler(myMessageOutput);
    QCoreApplication a(argc, argv);

  //  FFASettings ffaset;

    //FFA testForLoop;
   // mydisplay dd;

    QVector<float> test1({1., 3., 3., 0., 5, 3., 4., 2., 3., 1., 2., 4., 1., 0., 5, 3., 4., 2., 3.});


    MyQVector<float> myvect1(5,4);

    myvect1.arange(1,12,1);

    myvect1.show();

    myvect1.reshape(3,4);

    myvect1.show();

    myvect1.transpon();

    myvect1.show();

    MyQVector<float> mxv = myvect1.minAxis(1);

    mxv.show();

    myvect1.transpon();


    myvect1.show();

    myvect1.masked_where(5);

    myvect1.show_masked();

    myvect1(1,2) = 23;

    QVector<float> testRand({0.82307711, 0.34292748, 0.51316598, 0.47400193, 0.53924523,0.45135059, 0.2796427 , 0.28789632, 0.32794904, 0.82076353});

    myvect1 = testRand;

    myvect1.show();

    mxv = myvect1.medianAxis();

    mxv.show();

//    mxv.arange(-10);

//    mxv.show();

//    mxv = mxv*mxv;

//    mxv.show();



//    mxv.show_masked();

//    qDebug() << myvect1.getColumn(0);

//    qDebug() << myvect1.getColumn(1);

//    qDebug() << myvect1.getColumn(2);

//    qDebug() << myvect1.getColumn(-1);

//    qDebug() << myvect1.getColumn(10);

//    qDebug() << myvect1.getRow(0);

//    qDebug() << myvect1.getRow(1);

//    qDebug() << myvect1.getRow(-1);

    return 0;




//    QVector<float> test2(19,111);
//    ffa_tools ft;

//    QVector< QVector<float> > matrrix1;
//    QVector< QVector<float> > matrrix2;

//    qDebug() << "test1";
//    dd.dispVector(test1);
//    qDebug() << "test2";
//    dd.dispVector(test2);

//    ft.reshape(test1, matrrix1, 5, 10.f);
//    ft.reshape(test2, matrrix2, 5, 11.f);

//    qDebug() << "matrrix1";
//    dd.dispMatr(matrrix1);
//    qDebug() << "matrrix2";
//    dd.dispMatr(matrrix2);



////    int M = matrrix2.at(2).size();
////    int N = matrrix1.size();
////    int nCol = matrrix1.first().size();

////    int maxAvailSize = qMin(M,nCol);
////    qDebug() << M << nCol << maxAvailSize;

////    ft.replaceSubMatr(matrrix1, matrrix2.at(2),0,4,-1);
////    dd.dispMatr(matrrix1);

////    qDebug() << "standardDeviationInRow(matrrix1)" << ft.standardDeviationInRow(matrrix1);

////    qDebug() << "standardDeviation" << ft.standardDeviation(matrrix1.at(0));

////    QVector<bool> mask =  ft.masked_where(test1, 3.f);
////    dd.dispVector(mask);

////    QVector< QVector<bool> > maskmtr = ft.masked_where(matrrix1, 3.f);
////    dd.dispMatr(maskmtr);

////    dd.dispMatr(ft.MatrTranspon(matrrix1));

////    QVector<float> devidend = ft.VectorDevision(test2,test1,mask);
////    qDebug() << "ft.VectorDevision(test2,test1)";
////    dd.dispVector(devidend);

//    QVector<float>  snr1 = ft.simple_SNR_A(matrrix1,0.5f,0.1f);
//    dd.dispVector(snr1);

//    qDebug() << "std::log(8)" << std::log2f(8);

    MyQVector<float> testqv(4,4);

//    testqv.fill(4);
//    testqv.show();
//    testqv.fill(4,3,5);
//    testqv.show();
//    testqv.fill(4,8);
//    testqv.show();
//    testqv.fill(10,8,1);
//    testqv.show();

    //testqv = test1;
   // testqv.show();


    //qDebug() << "testqv.getVector()" << testqv.getVector();


    MyQVector<float> testqv1(0,4);
//    testqv1.fill(6);
//    testqv1.show();


//    testqv1.arrange(11);

//    testqv = testqv1;

    testqv1.fill(4,3,5);

    testqv = testqv1;

    //testqv.show();

    testqv[1]=5.f;

   // testqv.show();

    //qDebug() << " testqv[1] " << testqv[1];

    testqv(2,2) = 32.f;
    testqv(1,3) = 13.f;
    testqv(0,0) = 0.0f;
    testqv(0,1) = 0.1f;
    testqv(0,2) = 0.2f;
    testqv(0,3) = 0.3f;
    testqv(1,0) = 1.0f;
    testqv(1,1) = 1.1f;
    testqv(1,0) = 1.0f;

   // testqv1.show();

    testqv.show();

    MyQVector<float> sum;

    //sum.show();

  //  testqv+=testqv1;
    qDebug() << "Sum =";

   // testqv.show();
    testqv1.fill(10,3,5);
    sum = testqv+testqv1;
    qDebug() << "Sum2 =";
    sum.show();

    MyQVector<float> diff;

    testqv1.fill(2,3,4);

    testqv.fill(1,3,4);

    diff = testqv1 * testqv;

    qDebug() << "diff =";
    diff.show();

    testqv1(2,2) = 12.f;

    MyQVector<float> mult = testqv1*3.f - diff*2.f + 3;

    mult.show();

    QVector<float> tesgetvect = mult.getVector(1);
    qDebug() << "mult.getVector(): " << tesgetvect.size() << "  " << tesgetvect;

    qDebug() << "result of reshape " << mult.reshape(4,3);

    mult.show();

    //testqv1.arrange(12);

    testqv1.reshape(4,3);

   // testqv1 = ft.reshape(mult, 5,7);

    testqv1.show();

    //mult.show();

    testqv1.transpon();

    testqv1.show();

//    MyQVector<float> testqv2(4,0);
//    testqv2.fill(5);
//    testqv2.show();

//    testqv2.arrange(11);
//    testqv2.show();



//    qDebug() << rowMax.look_for_nan(test2);
//    qDebug() << rowMax.look_for_nan(test12);


//    test12[2]=11;
//    test12[7]=3;
//    //dd.dispVector(test2);
//    dd.dispVector(test12);

//    QVector<float> test3 = rowMax.add2Vect(test2,test12);
//    qDebug() << "rowMax.add2Vect(test2,test12)";
//    dd.dispVector(test3);

//    test12 = rowMax.multVectOnValue(test3, 2.f);

//    dd.dispVector(test12);


//    QVector< QVector<float> > matrrix2;

 //   rowMax.reshape(test2, matrrix2, 5, 10.f);


 //   dd.dispMatr(matrrix2);
   // qDebug() << rowMax.MatrMaxInRow(matrrix2);

    //qDebug() << rowMax.MatrMinInRow(matrrix2);

    //qDebug() << "median in row" << rowMax.MatrMedianInRow(matrrix2);

    //dd.dispMatr(rowMax.MatrTranspon(matrrix2));

    //qDebug() << rowMax.MatrMaxInRow(rowMax.MatrTranspon(matrrix2));

    //dd.dispVector( rowMax.MatrToVect(rowMax.MatrTranspon(matrrix2)));
//    QVector<float> test4({11., 12., 13., 14., 15, 16., 17., 18.});
//    dd.dispVector(test4);
//    rowMax.replaceSubVector(test4, test2, 2, 4);
//    dd.dispVector(test4);

    return 0;

    //QVector<float> test1({2., 4., 1., 3., 3., 2., 4., 3., 4., 2., 3., 1., 2., 4., 1., 3., 3., 2., 4., 3., 4., 2., 3., 1., 2., 4., 1.});
    //ffa_tools minmax;

    //qDebug() << "minmax.VecMax(test1) = " << minmax.VecMax(test1);
    //qDebug() << "minmax.VecMin(test1) = " << minmax.VecMin(test1);

    //return 0;

//    QVector<float> test; //({1., 3., 3., 2., 4., 3., 4., 2., 3., 1., 2., 4., 1., 5.});

//    testSignal tsig;

//    float period = 0.3f;
//    float dt = 0.05f;
//    int num = 40;

//    float testT = num*dt;

//    tsig.genSinus(test, num, period, dt);

    //float testT = 0.7;

    //testForLoop.startFFA(test, testT);

//    MatrixMulty testM;

//    qDebug() << testM.innerProduct(test,test2);

//    QVector<float> sm;

//    testM.innerSumm(test,test2,sm);

//    testM.constMulti(test,3,sm);


//    qDebug() << sm;


    //qDebug() << ffaset.p_ranges;

    //qDebug() << ffaset.dt_list;

    //qDebug() << ffaset.win_detrend;

    //    ffa_tools ft;




    //    cout << "factors test " << ft.factors(1000) << endl;

    //    numpy np;

    //    MatrixXf mnp = MatrixXf::Random(1,14);

    //    mnp << 1., 3., 3., 2., 4., 3., 4., 2., 3., 1., 2., 4., 1., 5.;

    //    cout  << "Numpy test reshape initial matrix mnp" << std::endl;
    //    cout  << mnp << std::endl;

    //    MatrixXf mresape = np.reshape(mnp, mnp.size()/2,2);

    //    //mresape.transposeInPlace();

    //    cout  << "Numpy after reshape : "  << endl;
    //    cout  << mresape << endl;

    //    MatrixXf  msum = mresape.rowwise().sum();
    //    cout  << msum << endl;

    //    MatrixXf m(2,2);
    //    m(0,0) = 3;
    //    m(1,0) = 2.5;
    //    m(0,1) = -1;
    //    m(1,1) = m(1,0) + m(0,1);




        //std::cout  << m << std::endl;

    //    VectorXf v(3);
    //    v << 1, 2, 3;

    //    QVector<float> y(4);
    //    memcpy(y.data(), m.data(), sizeof(float) * 4);

    //    QVector<float> z(3);

    //    memcpy(z.data(), v.data(), sizeof(float) * 3);

    //    //qDebug() << z;
    //    srand((unsigned int) time(0));
    //    MatrixXd m1 = MatrixXd::Random(3,3);
    //    m1 = (m1 + MatrixXd::Constant(3,3,1.2)) * 50;
    //    //std::cout << "m1 =" << std::endl << m1 << std::endl;
    //    VectorXd v1(3);
    //    v1 << 1, 2, 3;
    //    std::cout << "m1 * v1 = " << std::endl << m1 * v1 << std::endl;


    //    MatrixXf ca = MatrixXf::Random(3,2);
    //    cout << "Here is the matrix ca\n" << ca << endl;
    //    cout << "Here is the matrix ca^T\n" << ca.transpose() << endl;
    //    cout << "Here is the conjugate of ca\n" << ca.conjugate() << endl;
    //    cout << "Here is the matrix ca^*\n" << ca.adjoint() << endl;

    return a.exec();
}
