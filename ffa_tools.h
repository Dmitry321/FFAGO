#ifndef FFA_TOOLS_H
#define FFA_TOOLS_H

#include <QVector>
#include <QtMath>
#include <algorithm>
#include <QDebug>


class ffa_tools
{
public:
    explicit ffa_tools();

    int select_factor(QVector<float> &ts, const int &m, const int &mm);
    QVector<int> factors(int x);

    // Matrix creation with size nRow X mCol and fill int witch fill_value
    template<typename T>
    QVector< QVector<T> > initMatrix(const int &nRow, const int &mCol, T fill_value)
    {
        QVector< QVector<T> > mtr;
        mtr.fill(QVector<T>().fill(fill_value,mCol),nRow);
        return  mtr;
    }

//    template<typename T>
//    MyQVector<T> reshape(const MyQVector<T> &vec, int row, int col)
//    {
//        MyQVector<T> res(vec);
//        res.reshape(row,col);
//        return res;
//    }

//    QVector< QVector<bool> > initMatrixBool(const int &nRow, const int &mCol, bool fill_value = false)
//    {
//        QVector< QVector<bool> > mtr;
//        mtr.fill(QVector<bool>().fill(fill_value,mCol),nRow);
//        return  mtr;
//    }



// http://amin-ahmadi.com/2016/01/19/how-to-find-min-and-max-in-a-vector-qt-stl/
    template<typename T>
    T VecMin(const QVector<T> &vec)
    {
        T min = *std::min_element(vec.constBegin(), vec.constEnd());
        return min;
    }

    template<typename T>
    T VecMax(const QVector<T> &vec)
    {
        T max = *std::max_element(vec.constBegin(), vec.constEnd());
        return max;
    }

    template<typename T>
    QVector<bool> masked_where(const QVector<T> &vec, const T condition, bool maskVal = true)
    {
        int N = vec.size();
        QVector<bool> maskvec(N,!maskVal);
        int k=0;
        for (const T &v : qAsConst(vec))
        {
            if(qFuzzyCompare(v,condition))
            {
                maskvec[k] = maskVal;
            }
            k++;
        }
        return  maskvec;
    }

    template<typename T>
    QVector< QVector<bool> > masked_where(QVector< QVector<T> > &matr, const T condition, bool maskVal = true)
    {
        int nRow = matr.size();
        int nCol = matr.first().size();
        QVector< QVector<bool> >  maskmatr;

        maskmatr = initMatrix(nRow,nCol,false);
        int k=0;
        for (const QVector<T> &v : qAsConst(matr))
        {
            maskmatr[k] = masked_where(v, condition, maskVal);
            k++;
        }
        return maskmatr;
    }


    template<typename T>
    QVector<int> MatrMaxIndexInRow(const QVector<QVector<T> > &matr)
    {
        int M = matr.size();
        int k=0;
        QVector<int> resInd(M,0);
        T maxV;
        for (const QVector<T> &v : qAsConst(matr))
        {
            maxV = VecMax(v);
            resInd[k] = v.indexOf(maxV);
            k++;
        }
        return resInd;
    }

    template<typename T>
    QVector<int> MatrMinIndexInRow(const QVector<QVector<T> > &matr)
    {
        int M = matr.size();
        int k=0;
        QVector<int> resInd(M,0);
        T minV;
        for (const QVector<T> &v : qAsConst(matr))
        {
            minV = VecMin(v);
            resInd[k] = v.indexOf(minV);
            k++;
        }
        return resInd;
    }

    template<typename T>
    QVector<T> MatrMaxInRow(const QVector<QVector<T> > &matr)
    {
        int M = matr.size();
        int k=0;
        QVector<T> resMax(M,0);
        for (const QVector<T> &v : qAsConst(matr))
        {
            resMax[k] = VecMax(v);
            k++;
        }

// Deprecated
//        Q_FOREACH(const QVector<T> &v, matr)
//        {
//            resMax[k] = VecMax(v);
//            k++;
//            //resMax.append(VecMax(v));
//        }
        return resMax;
    }

    template<typename T>
    QVector<T> MatrMinInRow(const QVector<QVector<T> > &matr)
    {
        int M = matr.size();
        int k=0;
        QVector<T> resMin(M,0);
        for (const QVector<T> &v : qAsConst(matr))
        {
            resMin[k] = VecMin(v);
            k++;
        }
        return resMin;
    }

    template<typename T>
    T median(const QVector<T> &vector)
    {
        T med=static_cast<T>(0);
        int vsize = vector.size();
        if(vsize <=0)
            return med;
        if(vsize == 1)
            return vector.at(0);

        QVector<T> vsort = vector;
        //qSort(vsort);

        std::sort(vsort.begin(), vsort.end());
        if(vsize % 2 == 0)
        {
            med = (vsort.at(vsize/2 - 1) + vsort.at(vsize/2))/2.0;
        }
        else
        {
            med = vsort.at(vsize/2);
        }
        return med;
    }

    template<typename T>
    QVector<T> MatrMedianInRow(const QVector<QVector<T> > &matr)
    {
        int M = matr.size();


        int k=0;
        QVector<T> resMed(M,0);
        for (const QVector<T> &v : qAsConst(matr))
        {
            resMed[k] = median(v);
            k++;
        }


//        Q_FOREACH(const QVector<T> &v, matr)
//        {
//            resMed[k] = median(v);
//            k++;
//        }
        return resMed;
    }

    template<typename T>
    QVector<T> add2Vect(const QVector<T> &v1, const QVector<T> &v2, int znak=1, QVector<bool> mask = QVector<bool>(1,false))
    {
        int M = v1.size();
        int S = mask.size();
        //int N = v2.size();
        Q_ASSERT_X(M == v2.size(), "QVector<T> add2Vect(const QVector<T> &v1, const QVector<T> &v2)",
                   qPrintable(" Lenght of vector1 must be the same as cector2 "));
        int k=0;
        T cznak = static_cast<T>(znak);
        if(S!=M)
        {
            mask.fill(false, M);
        }
        QVector<T> res(M,0);
        for (const T &val2 : qAsConst(v2))
        {
            if(!mask.at(k))
            {
                res[k] = v1.at(k) + cznak * val2;
            }
            k++;
        }
        return res;
    }

    template<typename T>
    QVector<T> addValToVec(const QVector<T> &vec, const T &constValue, int vecZnak = 1)
    {
        int n = vec.size();
        QVector<T> res(n);
        int k=0;
        for (const T &val2 : qAsConst(vec))
        {
            res[k] = static_cast<T>(vecZnak) * val2 + constValue;
            k++;
        }
        return res;
    }

    template<typename T>
    QVector<T> multVectOnValue(const QVector<T> &vec, const T &constValue)
    {
        int n = vec.size();
        QVector<T> res(n);
        int k=0;
        for (const T &val2 : qAsConst(vec))
        {
            res[k] = val2 * constValue;
            k++;
        }
        return res;
    }

    template<typename T>
    QVector<T> devVectOnValue(const QVector<T> &vec, const T &constValue)
    {
        Q_ASSERT_X(!qFuzzyCompare(constValue, static_cast<T>(0)), "QVector<T> devVectOnValue(const QVector<T> &vec, const T &constValue)",
                   qPrintable(" Devision by zero "));

        int n = vec.size();
        QVector<T> res(n);
        int k=0;
        if(qFuzzyCompare(constValue, static_cast<T>(0)))
        {
            return res;
        }
        for (const T &val2 : qAsConst(vec))
        {
            res[k] = val2 / constValue;
            k++;
        }
        return res;
    }

    template<typename T>
    QVector<T> VectorDevision(const QVector<T> &numerator, const QVector<T> &denumerator, QVector<bool> mask = QVector<bool>(1,false))
    {
        int N = numerator.size();
        int M = denumerator.size();
        int S = mask.size();

        if(N!=M)
        {
            return QVector<T>(1,0);
        }
        QVector<T> devres(N,0);
        if(S!=N)
        {
            mask.fill(false, N);
        }
        int k = 0;
        for (const T &v : qAsConst(denumerator))
        {
            if(qFuzzyCompare(0.0,v))
                qWarning("Devision by zero in vector index: %d", k);
            if(!mask.at(k) && !qFuzzyCompare(0.0,v))
            {
                devres[k] = numerator.at(k)/v;
            }
            k++;
        }
        return devres;
    }

    template<typename T>
    QVector<T> VectorMultiplication(const QVector<T> &vec1, const QVector<T> &vec2, QVector<bool> mask = QVector<bool>(1,false))
    {
        int N = vec1.size();
        int M = vec2.size();
        int S = mask.size();

        if(N!=M)
        {
            return QVector<T>(1,0);
        }

        QVector<T> multres(N,0);

        if(S!=N)
        {
            mask.fill(false, N);
        }

        int k = 0;
        for (const T &v1 : qAsConst(vec1))
        {
            if(!mask.at(k) )
            {
                multres[k] = vec2.at(k)*v1;
            }
            k++;
        }
        return multres;
    }


    template<typename T>
    QVector<QVector<T> > MatrTranspon(const QVector<QVector<T> > &matr)
    {
        int nRow = matr.size();
        int nCol = matr.first().size();
        QVector<QVector<T> > trans;
        trans = initMatrix(nCol,nRow,static_cast<T>(0));

       // trans.fill(QVector<T>().fill(0,nRow),nCol);

        int k = 0;
        int j = 0;
        for (const QVector<T> &v : qAsConst(matr))
        //Q_FOREACH(const QVector<T> &v, matr)
        {
            j=0;
            for (const T &val : qAsConst(v))
            //Q_FOREACH(const T &val, v)
            {
                trans[j][k] = val;
                j++;
            }
            k++;
        }
        return trans;
    }

    template<typename T>
    QVector<T>  MatrToVect(const QVector<QVector<T> > &matr)
    {
        int nRow = matr.size();
        int nCol = matr.first().size();
        QVector<T>  trans(nRow*nCol,0);

        int k = 0;
        for (const QVector<T> &v : qAsConst(matr))
        //Q_FOREACH(const QVector<T> &v, matr)
        {
            for (const T &val : qAsConst(v))
            //Q_FOREACH(const T &val, v)
            {
                trans[k] = val;
                k++;
            }
        }
        return trans;
    }

    template <typename T>
    void downsample(QVector<T> &vector, const int &factor/*, QVector<T> &res*/)
    {
        if(factor>0)
        {
             QVector<T> res;
            //qDebug() << "downsample " << vector.length() % factor;

            Q_ASSERT_X(vector.length() % factor == 0 , "QVector<float> FFA::downsample(QVector<float> &vector, float &factor)",
                       qPrintable(" Lenght of vector " + QString::number(vector.size()) + " is not divisible by factor " + QString::number(factor)));

            int downLength = vector.size()/factor;

            res.fill(0.0, downLength);

            for(int i=0; i<downLength; i++)
            {
                for(int j=0; j<factor; j++)
                {
                    res[i] += vector.at(j+i*factor);
                }
            }
            vector = res;
        }
    }
    void reshape(const QVector<int> &x, QVector<QVector<int> > &res, int nrow, int fv);
    void reshape(const QVector<double> &x, QVector<QVector<double> > &res, int nrow, double fv);
    void reshape(const QVector<float> &x, QVector<QVector<float> > &res, int nrow, float fv);
    template<typename T>
    void normalize(QVector<T> &vector)
    {
        //QVector<T> vsort = vector;
       // qSort(vsort); deprecated

        //std::sort(vsort.begin(), vsort.end());

        T vMax = VecMax(vector);
                //vsort.last();
        //for (const T &val : qAsConst(vector))
        for(int i=0; i<vector.size(); i++)
        {
            vector[i] =  vector.at(i)/vMax;
        }
    }



    template<typename T>
    void vector_range_select(const QVector<T> &src, const T &minr, const T &maxr, QVector<T> &res)
    {
        res.resize(0);
        for (const T &v : qAsConst(src))
        //Q_FOREACH(const T &v, src)
        {
            if(v>=minr && v<=maxr)
                res.append(v);
        }
    }

    template<typename T>
    T standardDeviation(const QVector<T> &vector, T xavg, int n)
    {
        T sum = 0;
        for (const T &xi : qAsConst(vector))
        //Q_FOREACH(const T &xi, vector)
        {
            sum += (xi - xavg) * (xi - xavg);
        }
        //sum /= static_cast<T>(n);
        //sum = qSqrt(sum);
        return  qSqrt(sum / static_cast<T>(n));
    }

    template <typename T>
    T meanVec(const QVector<T>& vector)
    {
        int n = vector.size();
        T avg = std::accumulate(vector.begin(), vector.end(), static_cast<T>(0)) / static_cast<T>(n);
        return avg;
    }

    template <typename T>
    QVector<T> MatrMeanInRow(const QVector< QVector<T> >& matr)
    {
        int M = matr.size();
        int k=0;
        QVector<T> resMean(M,0);
        for (const QVector<T> &v : qAsConst(matr))
        {
            resMean[k] = meanVec(v);
            k++;
        }
        return resMean;
    }

    template <typename T>
    QPair<T, T> standardDeviation(const QVector<T>& vector)
    {
        int n = vector.size();
        T avg = meanVec(vector); //std::accumulate(vector.begin(), vector.end(), static_cast<T>(0)) / static_cast<T>(n);
        return QPair<T, T>(standardDeviation(vector, avg, n), avg);
    }

    template <typename T>
    QVector<T> standardDeviationInRow(const QVector< QVector<T> >& matr)
    {
        int nRow = matr.size();
        //int nCol = matr.first().size();
        QPair<T, T> stdmean;
        int k = 0;
        QVector<T> res(nRow);
        for (const QVector<T> &v : qAsConst(matr))
        {
            stdmean = standardDeviation(v);
            res[k] = stdmean.first;
            k++;
        }
        return res;
    }


    template <typename T>
    void maxMinInWindow(const QVector<T>& vector, const int &windSize, QVector<T> &minV, QVector<T> &maxV)
    {
        int n = vector.size();
        minV.resize(n);
        maxV.resize(n);
        Q_ASSERT_X(windSize <= n, "void maxMinInWindow(const QVector<T>& vector, const int &windSize)",
                   qPrintable(" Window size bigger then time series data size "));

        QVector<T> tmp(windSize,0);
        int ktmp = 0;
        int nextWin = 0;
        for(int i=0; i<n; i++)
        {
            tmp[ktmp] = vector.at(i);
            ktmp++;
            if(ktmp >= windSize)
            {
                ktmp=0;
                T maxv = VecMax(tmp);
                T minV = VecMin(tmp);
                for(int j = 0; j< windSize; j++)
                {
                    int cur = j+nextWin*windSize;
                    if(cur >= n)
                        break;
                    minV[cur] = minV;
                    maxV[cur] = maxV;
                }
                nextWin++;
                tmp.fill(0,windSize);
            }
        }

    }

    template <typename T>
    void np_arange(const int &M, QVector<T> &vec)
    {
        vec.resize(M);
        for(int i = 0; i < M; i++)
        {
            vec[i] = static_cast<T>(i);
        }
    }


//#==================	Signal-to-noise functions	==================
//if Config_ffa.metric == 'A':
//	def simple_SNR(folds, sigma, added_profs):
//    		""" Return a very simple signal-to-noise for a profile.
//        	Works for narrow duty-cycle since  the S/N=Max-Med/std
//		For each M profiles, returns a value of SNR (i.e, output is a list of lenght M)
//    		"""
    template <typename T>
    QVector<T> simple_SNR_A(const QVector<QVector<T> > &folds, float sigma, float added_profs)
    {
        int M = folds.size();
        //int P0 = folds.first().size();
        T prof_std = static_cast<T>(1.0)/(static_cast<T>(qSqrt(static_cast<qreal>(M - added_profs))*static_cast<qreal>(sigma)));
        QVector<T> snr;
        snr = add2Vect(MatrMaxInRow(folds), MatrMedianInRow(folds), -1);

        snr = multVectOnValue(snr, prof_std);
        Q_ASSERT_X(!look_for_nan(snr), "QVector<T> simple_SNR(const QVector<QVector<T> > &folds, float sigma, float added_profs)",
                   qPrintable(" Signal to noise = nan "));

        return snr;
    }

    template <typename T>
    void simple_SNR_off_pulse(const QVector<QVector<T> > &folds, QVector<QVector<T> > &off_pulse)
    {
        int M = folds.size();
        int P0 = folds.first().size();

        int width = int(0.1*P0);
        if( width < 1)
            width = 1;
        QVector<int> vwidth(M, width);

        //qDebug() << "QVector<T> vwidth(M, width)" << vwidth.size();
        QVector<int> testIndex = MatrMaxIndexInRow(folds);
        //qDebug() << "MatrMaxIndexInRow(folds)" << testIndex;
        //qDebug() << "MatrMaxIndexInRow(folds) size " << testIndex.size();
        QVector<int> id_low = add2Vect(MatrMaxIndexInRow(folds),vwidth,-1);                   //.argmax(axis=1)-width
        QVector<int> id_high = add2Vect(MatrMaxIndexInRow(folds),vwidth);
        //QVector< QVector<T> > off_pulse;
        off_pulse = initMatrix(M, P0, static_cast<T>(0));
        //qDebug() << "id_low" << id_low << "id_high" << id_high;
        for (int i=0; i< M; i++)
        {
            if( (id_low.at(i) >=0) && (id_high.at(i)<=P0) )
            {
                replaceSubMatr(off_pulse, folds.at(i), i, 0, id_low.at(i));
                replaceSubMatr(off_pulse, folds.at(i), i, id_high.at(i), -1);
                //off_pulse[i][0:id_low[i]] = folds[i][0:id_low[i]]
                //off_pulse[i][id_high[i]:-1] = folds[i][id_high[i]:-1]
            }
            else if ((id_low.at(i) <0) && (id_high.at(i)<=P0))
            {
                replaceSubMatr(off_pulse, folds.at(i), i, id_high.at(i), id_low.at(i));
                //off_pulse[i][id_high[i]:id_low[i]] = folds[i][id_high[i]:id_low[i]]
            }
            else if ( (id_low.at(i) >=0) && (id_high.at(i)>P0) )
            {
                int hi = id_high.at(i)-P0;
                replaceSubMatr(off_pulse, folds.at(i), i, hi, id_low.at(i));
                //off_pulse[i][hi:id_low[i]] = folds[i][hi:id_low[i]]
            }
        }
    }

    /*
      """ Return a very simple signal-to-noise for a profile.
            Works for narrow duty-cycle since  the S/N is max_value/std
        For each M profiles, returns a value of SNR (i.e, output is a list of lenght M)
        """
    */
    template <typename T>
    QVector<T> simple_SNR_B(const QVector<QVector<T> > &folds, float sigma, float added_profs)
    {
        int M = folds.size();
        QVector< QVector<T> > off_pulse;
        simple_SNR_off_pulse(folds, off_pulse);

        //qDebug() << "            simple_SNR_B" << M << " " << off_pulse.size();

        QVector<T> onesVec(M,static_cast<T>(1));
        QVector<T> std_off_pulse = standardDeviationInRow(off_pulse);
        QVector<bool> masked = masked_where(std_off_pulse, static_cast<T>(0));  //                    masked =  np.ma.masked_where(off_pulse == 0, off_pulse)
        QVector<T> prof_std = VectorDevision(onesVec,std_off_pulse,masked);  //        prof_std = np.ones(M)/(off_pulse.std(axis=1))
        QVector<T> snr;

        snr = add2Vect(MatrMaxInRow(folds),MatrMedianInRow(off_pulse),-1,masked); //                    snr = (folds.max(axis=1)-np.ma.median(masked,axis=1))*prof_std
        snr = VectorMultiplication(snr, prof_std,masked);

        //look_for_nan(snr);
        Q_ASSERT_X(!look_for_nan(snr), "QVector<T> simple_SNR(const QVector<QVector<T> > &folds, float sigma, float added_profs)",
                   qPrintable(" Signal to noise = nan "));
        return snr;

    }

//    def simple_SNR(folds, sigma_total, added_profs):
//                """ Return a very simple signal-to-noise for a profile.
//            For each M profiles, returns a value of SNR (i.e, output is a list of lenght M)
//            When calculating the median and the standard deviation of a profile, it exludes
//            a 20% window around the peak. S/N = Max-Med/std
//                """
    template <typename T>
    QVector<T> simple_SNR_C(const QVector<QVector<T> > &folds, float sigma_total, float added_profs)
    {
        int M = folds.size();
        //int P0 = folds.first().size();
        QVector< QVector<T> > off_pulse;
        simple_SNR_off_pulse(folds, off_pulse);

        QVector<T> onesVec(M,static_cast<T>(1));

        //QVector<T> std_off_pulse = standardDeviationInRow(off_pulse);
        //QVector<bool> masked = masked_where(std_off_pulse, static_cast<T>(0));

        T denumerator = sigma_total * static_cast<T>(qSqrt(0.8*static_cast<double>(M - added_profs)));
        denumerator = static_cast<T>(1)/denumerator;
        QVector<T> prof_std = multVectOnValue(onesVec,denumerator);  //  prof_std = np.ones(M)/(sigma_total*np.sqrt((0.8*(M-added_profs))))
        QVector<T> snr;
        snr = add2Vect(MatrMaxInRow(folds),MatrMedianInRow(off_pulse),-1);
        snr = VectorMultiplication(snr, prof_std);  //  snr = (folds.max(axis=1)-np.ma.median(masked,axis=1))*prof_std
        Q_ASSERT_X(!look_for_nan(snr), "QVector<T> simple_SNR(const QVector<QVector<T> > &folds, float sigma, float added_profs)",
                   qPrintable(" Signal to noise = nan "));
        return snr;
}
    // Replace sub vector in vector for sub colunm interval
    // from i_start to i_end, i_enf not included. If i_end < 0 it mean - for end of vector.
    template <typename T>
    void replaceSubVector(QVector<T> &lvect, const QVector<T> &rvect, const int &i_start,  const int &i_end)
    {
        int M = rvect.size();
        int N = lvect.size();
        int maxAvailSize = qMin(M,N);
        int ci_end = i_end;
        if(i_end < 0)
        {
            ci_end = maxAvailSize;
        }
        if(i_start < ci_end && i_start>=0 && ci_end<=maxAvailSize)
        {
            for(int i = i_start; i < ci_end; i++)
            {
                lvect[i] = rvect.at(i);
            }
        }
    }

    // Replace sub vector in 2D matrix on row - "row" and sub colunm interval
    // from i_start to i_end, i_enf not included. If i_end < 0 it mean - for end of row.
    template <typename T>
    void replaceSubMatr(QVector< QVector<T> > &lmatr, const QVector<T> &rvect, const int & row, const int &i_start, const int &i_end)
    {
        int M = rvect.size();
        int N = lmatr.size();
        int nCol = lmatr.first().size();
        int ci_end = i_end;
        int crow = row;

        int maxAvailSize = qMin(M,nCol);
        if(i_end < 0)
        {
            ci_end = nCol + i_end + 1;
        }
        // if row = -1 then use last row
        if(row < 0)
        {
            crow = N + row;
        }
        if(i_start==ci_end)
            return;

        Q_ASSERT_X(crow >=0 && crow < N && i_start>=0 && i_start < ci_end  && ci_end <= maxAvailSize, "void replaceSubMatr(QVector< QVector<T> > &lmatr, const QVector<T> &rvect, const int & row, const int &i_start, const int &i_end)",
                   qPrintable(" Some error in index sub matrix. LMatrix size " + QString::number(N) + "x" + QString::number(nCol) +
                              ", vector replase size " + QString::number(M) + " row for replays " + QString::number(crow) + " i_start " +
                              QString::number(i_start) +  " i_end " + QString::number(ci_end) ));

        if(crow >=0 && crow < N && i_start>=0 && i_start < ci_end  && ci_end <= maxAvailSize)
        {
            for(int i = i_start; i < ci_end; i++)
            {
                lmatr[crow][i] = rvect.at(i);
            }
        }
    }


    template <typename T>
    bool look_for_nan(const QVector<T> &vec)
    {
        bool flag = false;
        for (const T &xi : qAsConst(vec))
        {
            if( std::isnan(xi) )
            {
                flag = true;
            }
        }
        return flag;
    }


//    		M, P0 = folds.shape
//    		prof_std = 1.0/(np.sqrt(M-added_profs)*sigma)
//    		snr = (folds.max(axis=1)-np.median(folds, axis=1))*prof_std
//    		look_for_nan(snr)
//return snr

private:
    template<typename T>
    void reshapeT(const QVector<T> &x, QVector<QVector<T> > &res, int nrow, T fill_val)
    {
        res.resize(0);
        QVector<T> xx = x;
        int initSize = x.size();
        int colSize = 0;
        int add_position = initSize;

        if(initSize % nrow !=0 )
        {

            while(add_position % nrow !=0 )
            {
                add_position++;
                xx.append(fill_val);
            }
            //qDebug() << "          add_position " << add_position;
        }
        //qDebug() << xx;
        colSize = add_position/nrow;
        for(int i = 0; i < colSize; i++)
        {

            res.append(xx.mid(0+nrow*i,nrow));
        }
    }
};

#endif // FFA_TOOLS_H
