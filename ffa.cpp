#include "ffa.h"

#include <QtMath>
#include "ffa_stages.h"
#include "mydisplay.h"




FFA::FFA(QObject *parent) : QObject(parent)
{

}

void FFA::startFFAtest(QVector<float> &ts, float &T)
{

    ffa_tools ft;
    QVector< QVector<float> > tss;

    QVector<float> test({2., 11.,  0.,  1.,  3.,  3.,  2.,  4.,  3.,  4.,  2.,  3.,  1.,  2.,  4.,  1.,  5.,  6.});

    //qDebug() << test;

    //qDebug() << "after sort " << ts;

    //qDebug() << "reshape test ";

    ffa_stages FFAstage;
    float fillvl = 0.0;
    //FFAstage.XWrap2(test,tss,7,fillvl);

    mydisplay dd;

    //dd.dispMatr(tss);

    //int start = 1;
    //int stop = 2;
    //int len  = stop - start;
    //qDebug() << "Select range ";
    //dd.dispMatr(tss.mid(start, len));


    FFAstage.XWrap2(test,tss,7,fillvl,true);

    //qDebug() << "XWrap2";

    //dd.dispMatr(tss);

    FFAstage.FFA(tss);

    //qDebug() << tss;

    //ft.normalize(ts);


    //qDebug() << "Median= " << ft.median(ts);

    //qDebug() << "Normaliz = " << ts;

//    QPair<float, float> stdavg = standardDeviation(ts);
//    qDebug() << stdavg << "stdavg.first = " << stdavg.first;

//    downsample(ts, 2);
//    qDebug() << "after downsample";
//    qDebug() << ts;
}


// T    : float, the lenght of the observation (T in sec)
void FFA::startFFA(QVector<float> &ts, float &T)
{
    ffa_tools ft;
    ft.normalize(ts);
    // dt - Width of each time series bin (sec)
    float dt = T/ts.size();
    //int N = int(T/dt);

    int minimum_dwn = 0;
    int maximum_dwn = 0;
    QPair<float, float> stdavg = ft.standardDeviation(ts);
    float sigma_total = stdavg.first;

//# count_lim: used in stage 2 and 3; how many consecutive downsamplings
    int count_lim = 2;
    int dwn_ideal = 0;
//# Going through subsequent sub-ranges of periods (set in config_ffa.cfg)
//    # Each range of periods has it own initial sampling interval
    for(int num=0; num <  ffs.p_ranges.size(); num++)
    {
        if( num > 0)
        {
            dwn_ideal = int(ffs.dt_list.at(num)/dt);
            if( (dwn_ideal == 1)  || (dwn_ideal == 0) )
                dwn_ideal = 2;

            minimum_dwn = dwn_ideal - int(dwn_ideal*0.05);
            maximum_dwn = dwn_ideal + int(dwn_ideal*0.15);

            int dwn = ft.select_factor(ts,minimum_dwn,maximum_dwn);

            // // test only
            qDebug() << "dwn= " << dwn << "minimum_dwn " << minimum_dwn << "maximum_dwn " << maximum_dwn << "dt = " << dt << "T = " << T;

            ft.downsample(ts, dwn);

            sigma_total *= qSqrt(dwn);
            dt = T/static_cast<float>(ts.length());
        }

        qDebug() << "  Folding, period range of " << ffs.p_ranges.at(num) << "...";

        //all_SNs_x1, all_Ps_x1, dts_x1 =  fs.ffa_code_stage1(ts, dt,T, sigma_total,p_ranges[num][0],\
        //p_ranges[num][1], count_lim,name)

        ffa_code_stage1(ts, dt,T, sigma_total,ffs.p_ranges.at(num).at(0),ffs.p_ranges.at(num).at(1), count_lim);



    }



}

//def ffa_code_stage1(data ,dt , T , sigma_total,p_min, p_max, count_lim, name):
//	"""
//	ffa_code_stage1 (data , dt , T, period , N , p_min , p_max , count_lim , name):
//		- data		:  Time series
//		- dt		:  Sampling interval (s)
//		- T		:  Total observative time (s)
//		- p_min		:  Minimum period in the subset of trial periods (s)
//		- p_max		:  Maximum period in the subset of trial periods (s)
//		- count_lim	:  int 1 or 2, used in stage2 and stage 3
//				   if count_lim =1, goes to 4*dt and 9*dt
//				   if count_lim =2, goes to 8*dt and 27*dt
//		- name		:  Name of the beam (without the extension)
//	Returns
//	"""
//	# --------------------	   FFA Stage 1	----------------------

void FFA::ffa_code_stage1(QVector<float> &data, const float &dt, const float &TT, const float &sigma_total, const float &p_min , const float &p_max, const int &count_lim)
{
    ffa_tools ft;
    ffa_stages fs;
    mydisplay dd;
    float fill_value = ft.median(data);
    int N = int(TT/dt);
    int P0_start = qFloor(static_cast<qreal>(p_min)/static_cast<qreal>(dt));
    int P0_end = qCeil(static_cast<qreal>(p_max)/static_cast<qreal>(dt));
    //SNs1 = []
    // all_Ps1 = []
    for(int p0 = P0_start; p0 < P0_end; p0++)
    {
        float  M_real = static_cast<float>(N)/p0;
        int stepen = qFloor(static_cast<double>(std::log2f(M_real))) + 1;
        //qDebug() << std::log(M_real);
        float added_profs = static_cast<float>(qPow(2,stepen)) - M_real; //            added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
        if(p0==0 || p0 ==1)
        {
            qDebug() << "         It tried to fold with period = 0 or 1 bin";
            continue;
        }
        //float p0_sec = p0*dt;
        QVector< QVector<float> > xwrap;
        fs.XWrap2(data,xwrap, p0, fill_value, true); //            xwrap = FFA.XWrap2(data,p0, fill_value=fill_value, pow2=True)
        //dd.dispVector(data);
        //dd.dispMatr(xwrap);
        QVector< QVector<float> > folds = fs.FFA(xwrap); //            folds = FFA.FFA(xwrap)
        dd.dispMatr(folds); //            M = folds.shape[0]
        int M = folds.size();

        QVector<float> SN = ft.simple_SNR_A(folds, sigma_total, added_profs); //            SN = f.simple_SNR(folds, sigma_total ,added_profs)
        QVector<float> arr_vec;
        ft.np_arange(M, arr_vec);

        QVector<float> P = ft.addValToVec(ft.devVectOnValue(arr_vec,static_cast<float>(M-1)), static_cast<float>(p0)); // P = p0 + (np.arange(M, dtype=np.float) / (M-1))
        QVector<float> Psec = ft.multVectOnValue(P,dt);


    }



//            if p0==0 or p0 ==1:
//                print '                  It tried to fold with period = 0 or 1 bin'
//                continue
//            p0_sec = p0*dt





//            Psec=P*dt
//            SNs1.append(SN)
//            all_Ps1.extend(Psec)
//        SNs1 = np.concatenate(SNs1)
//        dts = [dt]*len(SNs1)
//        return np.array([SNs1]),np.array([all_Ps1]), dts

//    #______________________________________________________________________


}

//def XWrap2(x,P0,fill_value=0,pow2=False):
//    """
//    Extend and wrap array.

//    Fold array every y indecies.  There will typically be a hanging
//    part of the array.  This is padded out.
//    Parameters
//    ----------
//    x     : input
//    P0    : Base period, units of elements
//    pow2  : If true, pad out nRows so that it's the next power of 2.
//    Return
//    ------
//    xwrap : Wrapped array.
//    """
//template <typename T>
//void FFA::XWrap2(QVector<T> &x, const int &P0, int fill_value, bool pow2)
//{
//    int ncad = x.size();
//    int nrow = qFloor(ncad/P0) + 1;
//    int nExtend = nrow * P0 - ncad; // # Pad out remainder of array with 0s.
//    QVector<T> pad;
//    pad.fill(fill_value);

//    x.append(pad);


//}




