#ifndef FFA_H
#define FFA_H

#include <QObject>
#include <QVector>
#include <QDebug>
#include <cmath>
#include "myqvector.h"
#include "mymath.h"
//#include <QtMath>

class FFA : public QObject
{
    Q_OBJECT
public:
    explicit FFA(QObject *parent = nullptr);

//    template<typename T>
//    MyQVector<T> XWrap2(const MyQVector<T> &x, int P0, T fill_value=0, bool pow2=false)
//    {
//        int ncad = x.size(); // # Number of cadences
//        Q_ASSERT_X(ncad >= P0 , "void XWrap2(const QVector<T> &x, QVector< QVector<T> > &xwrap, int P0, T fill_value=0, bool pow2=false)",
//                   qPrintable(" Lenght of vector = " + QString::number(ncad) + " less then period in point P0 = " + QString::number(P0) + "... "));
//        int nrow = ncad/P0;

//        MyQVector<T> xwrap = x;



//        xwrap.reshape(-1,P0);


//        nrow = xwrap.size(0);
//        //qDebug() << "ncad = " << ncad << " nrow = " << nrow << " P0 = " << P0;
//        if(pow2)
//        {
//            MyMath pw;
//            int k = static_cast<int>(std::ceil(std::log2(nrow)));
//            int nrow2 = static_cast<int>(pw.pow2(static_cast<uint64_t>(k)));
//          //  qDebug() << "k = " << k << " nrow2 = " << nrow2;
//            QVector< QVector<float> > mfill;
//            mfill.fill(QVector<float>().fill(fill_value,P0),nrow2-nrow);
//            xwrap.append(mfill);
//        }
//    }

signals:



private:

};

#endif // FFA_H
