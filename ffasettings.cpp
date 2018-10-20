#include "ffasettings.h"


FFASettings::FFASettings(QObject *parent) : QObject(parent),
    SN_tresh(6.f),
    mindc(0.5f),
    numdms(2),
    metric('A')
{
    //reallocateMassiv(6, 2, p_ranges);

    p_ranges.push_back(QVector<float>({0.1f, 0.5f}) );
    //p_ranges.push_back(QVector<float>({0.5f, 1.f}) );
    //p_ranges.push_back(QVector<float>({1.f, 2.f}) );
    //p_ranges.push_back(QVector<float>({2.f, 5.f}) );
    //p_ranges.push_back(QVector<float>({5.f, 10.f}) );
    //p_ranges.push_back(QVector<float>({10.f, 15.f}) );
    //p_ranges.push_back(QVector<float>({15.f, 30.f}) );

    dt_list = QVector<float>({0.002f, 0.005f, 0.01f, 0.02f, 0.05f, 0.075f});

    win_detrend = p_ranges.last().last() * 2.f; //# twice the largest trial period
}

void FFASettings::reallocateMassiv(int &row, int &col, QVector< QVector<float> > &marray)
{
    QVector<float> fillArray;
    // Готовим обнуленный массив для всех колонок
    fillArray.fill(0.0, col);

    //двумерный массив для хранения всех строк
    marray.fill(fillArray, row);
}
