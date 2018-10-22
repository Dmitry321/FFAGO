#include "ffago.h"

FFAGO::FFAGO(QObject *parent) : QObject(parent)
{
    
}

FFAGO::~FFAGO()
{
    
}
/* print "  Folding, period range of ", p_ranges[num], " ..."
        all_SNs_x1, all_Ps_x1, dts_x1 =  fs.ffa_code_stage1(ts, dt,T, sigma_total,p_ranges[num][0],\
p_ranges[num][1], count_lim,name)
*/
void FFAGO::ffaCodeStage1(MyQVector<float> &vector, const float &dt, const float &TT, const float &sigma_total, const float &p_min, const float &p_max, const int &count_lim)
{
    
}
