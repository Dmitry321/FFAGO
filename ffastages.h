#ifndef FFASTAGES_H
#define FFASTAGES_H

#include <QObject>
#include "myqvector.h"


class ffaStages : public QObject
{
    Q_OBJECT
public:
    explicit ffaStages(QObject *parent = nullptr);

    void ffaCodeStage1(MyQVector<float> &data, const float &dt, const float &T, const float &sigma_total, const float &p_min, const float &p_max, const int &count_lim);
signals:

public slots:


};

#endif // FFASTAGES_H
