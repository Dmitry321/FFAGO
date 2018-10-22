#ifndef FFAGO_H
#define FFAGO_H

#include <QObject>
#include "myqvector.h"
#include "ffasettings.h"

class FFAGO : public QObject
{
    Q_OBJECT
public:
    explicit FFAGO(QObject *parent = nullptr);
    ~FFAGO();

signals:

public slots:

private:
    MyQVector<float> all_SNs_x1, all_SNs1;
    MyQVector<float> all_Ps_x1,all_Ps1;
    MyQVector<float> all_SNs_x2,all_SNs_2;
    MyQVector<float> all_Ps_x2,all_Ps_2;
    MyQVector<float> all_SNs_x3,all_SNs_3;
    MyQVector<float> all_Ps_x3,all_Ps_3;
    MyQVector<float> SNs1;
    MyQVector<float> SNs2_phase1, SNs2_phase2;
    MyQVector<float> SNs4_phase1, SNs4_phase2 ,SNs4_phase3, SNs4_phase4;
    MyQVector<float> SNs8_phase1 ,SNs8_phase2 ,SNs8_phase3, SNs8_phase4;
    MyQVector<float> SNs8_phase5 ,SNs8_phase6 ,SNs8_phase7, SNs8_phase8;
    MyQVector<float> SNs3_phase1 ,SNs3_phase2, SNs3_phase3;
    MyQVector<float> SNs9_phase1 ,SNs9_phase2 ,SNs9_phase3 ,SNs9_phase4;
    MyQVector<float> SNs9_phase5 ,SNs9_phase6 ,SNs9_phase7 ,SNs9_phase8 ,SNs9_phase9;
    MyQVector<float> SNs27_phase1, SNs27_phase2, SNs27_phase3, SNs27_phase4;
    MyQVector<float> SNs27_phase5, SNs27_phase6, SNs27_phase7, SNs27_phase8;
    MyQVector<float> SNs27_phase9, SNs27_phase10, SNs27_phase11, SNs27_phase12;
    MyQVector<float> SNs27_phase13, SNs27_phase14, SNs27_phase15, SNs27_phase16;
    MyQVector<float> SNs27_phase17, SNs27_phase18, SNs27_phase19, SNs27_phase20;
    MyQVector<float> SNs27_phase21, SNs27_phase22, SNs27_phase23, SNs27_phase24;
    MyQVector<float> SNs27_phase25, SNs27_phase26,SNs27_phase27;
    MyQVector<float> Ps1, Ps2, Ps4, Ps8, Ps3, Ps9, Ps27;
    MyQVector<float> dts_x1,dts_1;
    MyQVector<float> dts_x2,dts_2;
    MyQVector<float> dts_x3,dts_3;
    MyQVector<float> dt1s,dt2s,dt4s,dt8s,dt3s,dt9s,dt27s;

    FFASettings ffs;

    void ffaCodeStage1(MyQVector<float> &vector, const float &dt, const float &TT, const float &sigma_total, const float &p_min, const float &p_max, const int &count_lim);
};

#endif // FFAGO_H
