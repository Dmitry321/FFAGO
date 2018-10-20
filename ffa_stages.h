#ifndef FFA_STAGES_H
#define FFA_STAGES_H
#include <QVector>
#include "ffa_tools.h"
#include <QDebug>

class ffa_stages
{
public:
    ffa_stages();
    QVector<QVector<float> > FFA(const QVector<QVector<float> > &XW);

    //def XWrap2(x,P0,fill_value=0,pow2=False):
    //    """
    //    Extend and wrap array.

    //    Fold array every y indecies.  There will typically be a hanging
    //    part of the array.  This is padded out.
    //    Parameters
    //    ----------
    //    x     : input vector
    //    P0    : Base period, units of elements
    //    pow2  : If true, pad out nRows so that it's the next power of 2.
    //    Return
    //    ------
    //    xwrap : Wrapped array.
    //    """


    template<typename T>
    void XWrap2(const QVector<T> &x, QVector< QVector<T> > &xwrap, int P0, T fill_value=0, bool pow2=false)
    {
        ffa_tools ffat;
        int ncad = x.size(); // # Number of cadences
        Q_ASSERT_X(ncad >= P0 , "void XWrap2(const QVector<T> &x, QVector< QVector<T> > &xwrap, int P0, T fill_value=0, bool pow2=false)",
                   qPrintable(" Lenght of vector = " + QString::number(ncad) + " less then period in point P0 = " + QString::number(P0) + "... "));
        int nrow = ncad/P0;

        ffat.reshape(x,xwrap,P0,fill_value);
        nrow = xwrap.size();
        //qDebug() << "ncad = " << ncad << " nrow = " << nrow << " P0 = " << P0;
        if(pow2)
        {
            int k = qCeil(std::log2(nrow));
            int nrow2 = static_cast<int>(qPow(2,k));
          //  qDebug() << "k = " << k << " nrow2 = " << nrow2;
            QVector< QVector<float> > mfill;
            mfill.fill(QVector<float>().fill(fill_value,P0),nrow2-nrow);
            xwrap.append(mfill);
        }
    }

    inline constexpr std::uint64_t pow2(std::uint64_t i) const;

private:
    QVector<QVector<float> > FFAShiftAdd(QVector<QVector<float> > &XW0, const int &stage);
    void FFAGroupShiftAdd(const QVector<QVector<float> > &group0, QVector<QVector<float> > &group, int nRowGroup, int nColGroup);
};

#endif // FFA_STAGES_H
