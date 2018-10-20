#ifndef FFASETTINGS_H
#define FFASETTINGS_H

#include <QObject>
#include <QVector>

class FFASettings : public QObject
{
    Q_OBJECT

public:
    explicit FFASettings(QObject *parent = nullptr);

    QVector< QVector<float> >  p_ranges;
    QVector<float> dt_list;
    float win_detrend; // "Time window for detrending (in sec)" twice the largest trial period
    float SN_tresh; // "Signal-to-noise treshold for picking candidates"
    float mindc; 	 //   = 0.5     	# min_dc = 0.5 , 1 or 1.5 (this is the min duty-cycle in %)
//    "Minimum duty-cycle to look for. Default is 0.5%. "\
//                "Options are: 0.5, 1., 1.5. It multiplies the list "\
//                "of minimum sampling intervals that has to be tested in each "\
//    "subranges of periods by 2 (mindc=1) or 3 (mindc =1.5)"
    int numdms;      //= 2         # relevant for sifting with ffa_final.py

    char metric; //= 'A' # 'A', 'B', or 'C', see Parent et al. 2018 for details on the metrics.
private:
    void reallocateMassiv(int &row, int &col, QVector<QVector<float> > &marray);

signals:

public slots:
};

#endif // FFASETTINGS_H
