#ifndef FFA_H
#define FFA_H

#include <QObject>
#include <QVector>
#include <QDebug>
#include "ffasettings.h"

#include "ffa_tools.h"



class FFA : public QObject
{
    Q_OBJECT
public:
    explicit FFA(QObject *parent = nullptr);

    void startFFAtest(QVector<float> &ts, float &T);
signals:

public slots:
    void startFFA(QVector<float> &ts, float &T);

private:
    QVector<float> all_SNs_x1, all_SNs1;
    QVector<float> all_Ps_x1,all_Ps1;
    QVector<float> all_SNs_x2,all_SNs_2;
    QVector<float> all_Ps_x2,all_Ps_2;
    QVector<float> all_SNs_x3,all_SNs_3;
    QVector<float> all_Ps_x3,all_Ps_3;
    QVector<float> SNs1;
    QVector<float> SNs2_phase1, SNs2_phase2;
    QVector<float> SNs4_phase1, SNs4_phase2 ,SNs4_phase3, SNs4_phase4;
    QVector<float> SNs8_phase1 ,SNs8_phase2 ,SNs8_phase3, SNs8_phase4;
    QVector<float> SNs8_phase5 ,SNs8_phase6 ,SNs8_phase7, SNs8_phase8;
    QVector<float> SNs3_phase1 ,SNs3_phase2, SNs3_phase3;
    QVector<float> SNs9_phase1 ,SNs9_phase2 ,SNs9_phase3 ,SNs9_phase4;
    QVector<float> SNs9_phase5 ,SNs9_phase6 ,SNs9_phase7 ,SNs9_phase8 ,SNs9_phase9;
    QVector<float> SNs27_phase1, SNs27_phase2, SNs27_phase3, SNs27_phase4;
    QVector<float> SNs27_phase5, SNs27_phase6, SNs27_phase7, SNs27_phase8;
    QVector<float> SNs27_phase9, SNs27_phase10, SNs27_phase11, SNs27_phase12;
    QVector<float> SNs27_phase13, SNs27_phase14, SNs27_phase15, SNs27_phase16;
    QVector<float> SNs27_phase17, SNs27_phase18, SNs27_phase19, SNs27_phase20;
    QVector<float> SNs27_phase21, SNs27_phase22, SNs27_phase23, SNs27_phase24;
    QVector<float> SNs27_phase25, SNs27_phase26,SNs27_phase27;
    QVector<float> Ps1, Ps2, Ps4, Ps8, Ps3, Ps9, Ps27;
    QVector<float> dts_x1,dts_1;
    QVector<float> dts_x2,dts_2;
    QVector<float> dts_x3,dts_3;
    QVector<float> dt1s,dt2s,dt4s,dt8s,dt3s,dt9s,dt27s;

    FFASettings ffs;

    void ffa_code_stage1(QVector<float> &vector, const float &dt, const float &TT, const float &sigma_total, const float &p_min, const float &p_max, const int &count_lim);


};

#endif // FFA_H
