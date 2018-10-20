#include "ffa_tools.h"


ffa_tools::ffa_tools()
{

}

QVector<int> ffa_tools::factors(int x)
{
    QVector<int> facts;
    for(int i = 1; i < x+1; i++ )
    {
        if (x % i == 0)
        {
            facts.append(i);
        }
    }
    return facts;
}


int ffa_tools::select_factor(QVector<float> &ts, const int &m, const int &mm)
{
//    """ This will delete the last element of ts until it has a factor within [m,mm]
//        There is a maximum number of element that you delete in order to get a factor.
//        select_factor(ts, m, mm):
//        Inputs:
//            ts :  array (time series usually)
//            m  : minimum factor
//            mm : maximum factor
//        Returns:
//            ts : array (possibly shorter)
//            x[0]: first factor that is within the range [minimum,maximum]
//        """
    QVector<int> a = factors(ts.length());   // a=np.array(factors(len(ts)))
    QVector<int> x;
    vector_range_select(a,m,mm,x);           // x = a[(a >=m) & (a <=mm)]
    qDebug() << "a=" << a;
    qDebug() << "x=" << x;
    int counts = 0;

    while(x.length()==0 && counts <50 && !ts.isEmpty())
    {
        ts.removeLast();
        //qDebug() << "ts= " << ts;
        a = factors(ts.length());
        vector_range_select(a,m,mm,x);
        counts+=1;
        if (counts >= 50)
        {
//            "Having a hard time finding a downsampling factor to match the desired time 			resolution. "
//                    print "Try changing : "+'\n'+'1) the desired time resolution '+'\n'+\
//                    ' 2) by how much you are welling to vary from that time resolution'+'\n'+\
//            ' 3) How much sample you are willing to delete (default = max 40) '
            return 0;
        }
    }
    if(x.length())
        return x.first();
    else
        return 0;
}

//template <typename T>
void ffa_tools::reshape(const QVector<float> &x, QVector<QVector<float> > &res, int nrow, float fv)
{
    reshapeT(x,res,nrow,fv);
//    if(pow2)
//    {
//        int P0 = res.first().size();
//        int k = qCeil(std::log2(nrow));
//        int nrow2 = static_cast<int>(qPow(2,k));


//        QVector< QVector<float> > mfill;
//        mfill.fill(QVector<float>().fill(fv,P0),nrow2-nrow);
//        qDebug() << nrow2 << mfill;
//        res.append(mfill);

//    }
}

void ffa_tools::reshape(const QVector<double> &x, QVector<QVector<double> > &res, int nrow, double fv)
{
    reshapeT(x,res,nrow,fv);
}

void ffa_tools::reshape(const QVector<int> &x, QVector<QVector<int> > &res, int nrow, int fv)
{
    reshapeT(x,res,nrow,fv);
}
