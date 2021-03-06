#include <QCoreApplication>
//#include "ffa.h"
#include <QDebug>
#include "testsignal.h"
#include "myqvector.h"
#include "ffastages.h"


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

    QVector<float> test1({1., 3., 3., 0., 5, 3., 4., 2., 3., 1., 2., 4., 1., 0., 5, 3., 4., 2., 3.});

    MyQVector<float> myvect1(5,4);

//    myvect1.arange(1,162,2);

//    myvect1.show();

//    myvect1.reshape(-1, 5,true,23.1f);

//    myvect1.show();

//    myvect1.transpon();

//    myvect1.show();

    MyQVector<float> mxv = myvect1.maxAxis(1);

    //mxv.show();

    ffaStages fs;

    testSignal ts;

    ts.genSinus(test1,20,0.3f,0.05f);

    myvect1 = test1;

    myvect1.show();

    mxv = myvect1.maxAxis();

    mxv.show();

    mxv = myvect1.meanAxis();

    mxv.show();

    float sigma_total = myvect1.standardDeviation();


    fs.ffaCodeStage1(myvect1, 0.05, 1, sigma_total, 0.2, 0.4, 40);


    return 0;

    myvect1.show();

    myvect1.masked_where(5);

    myvect1.show_masked();

    myvect1(1,2) = 23;

    QVector<float> testRand({0.82307711f, 0.34292748f, 0.51316598f, 0.47400193f, 0.53924523f,0.45135059f, 0.2796427f , 0.28789632f, 0.32794904f, 0.82076353f});

    myvect1 = testRand;

    myvect1.show();

    mxv = myvect1.medianAxis();

    mxv.show();

    mxv = myvect1.meanAxis(1);
    mxv.show();

   // myvect1.normalize();

   // myvect1.show();





    //    MatrixXf ca = MatrixXf::Random(3,2);
    //    cout << "Here is the matrix ca\n" << ca << endl;
    //    cout << "Here is the matrix ca^T\n" << ca.transpose() << endl;
    //    cout << "Here is the conjugate of ca\n" << ca.conjugate() << endl;
    //    cout << "Here is the matrix ca^*\n" << ca.adjoint() << endl;

    return a.exec();
}
