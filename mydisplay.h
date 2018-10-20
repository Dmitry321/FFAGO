#ifndef MYDISPLAY_H
#define MYDISPLAY_H
#include <QDebug>


class mydisplay
{
public:
    mydisplay();
    template<typename T>
    void dispMatr(const QVector<QVector<T> > &x)
    {
        qDebug() << "Matrix " << x.size() << "x" << x.first().size();
        QString strOut = "[";
        Q_FOREACH(const QVector<T> &xi, x)
        {
            strOut.append("[");
            Q_FOREACH(const T &val, xi)
            {
                strOut.append(QString::number(val));
                strOut.append(" ");
            }
            strOut.append("]");
            qDebug() << strOut;
            strOut =" ";
            //sum += (xi - xavg) * (xi - xavg);
        }
    }
    template<typename T>
    void dispVector(const QVector<T> &x)
    {
        qDebug() << "Vector " << x.size();
        QString strOut = "[";
        Q_FOREACH(const T &val, x)
        {
            strOut.append(QString::number(val));
            strOut.append(" ");
        }
        strOut.append("]");
        qDebug() << strOut;
    }
};

#endif // DISPLAY_H
