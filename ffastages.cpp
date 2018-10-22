#include "ffastages.h"
#include <QtMath>
#include "mymath.h"
#include "ffa.h"

ffaStages::ffaStages(QObject *parent) : QObject(parent)
{

}

/*
"""
Make histogram
check for why sigma isnt sig*sqrt{dwn}
"""
def ffa_code_stage1(data ,dt , T , sigma_total,p_min, p_max, count_lim, name):
"""
    ffa_code_stage1 (data , dt , T, period , N , p_min , p_max , count_lim , name):
        - data		:  Time series
        - dt		:  Sampling interval (s)
        - T		:  Total observative time (s)
        - p_min		:  Minimum period in the subset of trial periods (s)
        - p_max		:  Maximum period in the subset of trial periods (s)
        - count_lim	:  int 1 or 2, used in stage2 and stage 3
                   if count_lim =1, goes to 4*dt and 9*dt
                   if count_lim =2, goes to 8*dt and 27*dt
        - name		:  Name of the beam (without the extension)
    Returns
    """
*/
void ffaStages::ffaCodeStage1(MyQVector<float> &data, const float &dt, const float &T, const float &sigma_total, const float &p_min , const float &p_max, const int &count_lim)
{
// # --------------------	   FFA Stage 1	----------------------
    auto fill_value = data.medianV();                              // fill_value = np.median(data)
    float N = T/dt;                                                // N = T/dt


    int P0_start = qFloor(static_cast<qreal>(p_min)/static_cast<qreal>(dt)); //P0_start, P0_end = np.floor(p_min/dt), np.ceil(p_max/dt)
    int P0_end = qCeil(static_cast<qreal>(p_max)/static_cast<qreal>(dt));
    //SNs1 = []
    // all_Ps1 = []
    for(auto p0 = P0_start; p0 <= P0_end; p0++)  // P0s = np.arange(P0_start,P0_end,1) for p0 in P0s
    {
        float  M_real = N/p0;                    //  M_real = float(float(N)/p0)

        int stepen = qFloor(static_cast<double>(std::log2f(M_real))) + 1;
        float added_profs = static_cast<float>(pow2(static_cast<uint64_t>(stepen))) - M_real; //added_profs = 2**(int(math.floor(math.log(M_real,2)) + 1)) - M_real
        if(p0==0 || p0 ==1)   // if p0==0 or p0 ==1:
        {
            qDebug() << "         It tried to fold with period = 0 or 1 bin";
            continue;
        }
        float p0_sec = p0*dt;        //p0_sec = p0*dt
        MyQVector<float> xwrap =  data;
        xwrap.reshape(-1,p0,true,fill_value); //xwrap = FFA.XWrap2(data,p0, fill_value=fill_value, pow2=True)
        qDebug() << "void ffaStages::ffaCodeStage1 " << p0 << "added_profs = " << added_profs << "p0_sec = " << p0_sec;
        xwrap.show();
    }


}
/*
    for p0 in P0s:



        folds = FFA.FFA(xwrap)
        M = folds.shape[0]
        SN = f.simple_SNR(folds, sigma_total ,added_profs)
        P = p0 + (np.arange(M, dtype=np.float) / (M-1))
        Psec=P*dt
        SNs1.append(SN)
        all_Ps1.extend(Psec)
    SNs1 = np.concatenate(SNs1)
    dts = [dt]*len(SNs1)
    return np.array([SNs1]),np.array([all_Ps1]), dts

#______________________________________________________________________
*/
