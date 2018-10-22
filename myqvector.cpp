#include "myqvector.h"


//template<class A>
//bool MyQVector<A>::reshape(int row, int col)
//{
//    int allSize = nnRow*nnCol;
//    if(row <0 && col > 1)
//    {
//        row = allSize%col
//                ? allSize/col
//                : allSize/col + 1;
//    }
//    int newsize = row*col;
//    if(newsize > allSize)
//    {
//        QVector<A> tmp(newsize);
//        int k = 0;
//        for (const auto &v : qAsConst(matrix))
//        {
//            tmp[k]=v;
//            k++;
//        }
//        tmp.swap(matrix);
//        nnRow = qAbs(row);
//        nnCol = qAbs(col);
//        if(nnCol > 1)
//            flag_isMatrix=true;
//        else
//            flag_isMatrix=false;
//        return true;
//    }
//    if(newsize==allSize && nnCol!=col)
//    {
//        nnRow = qAbs(row);
//        nnCol = qAbs(col);
//        if(nnCol > 1)
//            flag_isMatrix=true;
//        else
//            flag_isMatrix=false;
//        return true;
//    }
//    return false;
//}

//template<class A>
//QVector<A> MyQVector<A>::getRow(int row)
//{
//    QVector<A> res(nnCol);
//    if(nnRow>1 && nnCol>1)
//    {
//        if(row<0 || row>=nnRow)
//            row=nnRow-1;
//        res = matrix.mid(row * nnCol, nnCol);
//    }
//    return res;
//}

