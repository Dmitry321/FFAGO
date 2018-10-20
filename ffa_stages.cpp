#include "ffa_stages.h"
#include <QDebug>
#include <cmath>
#include <QtMath>
#include "mydisplay.h"

ffa_stages::ffa_stages()
{

}

//def FFA(XW):
//    """
//    Fast Folding Algorithm
//    Consider an evenly-spaced timeseries of length N.  We can fold it
//    on P0, by creating a new array XW, shape = (P0,M), M = N/P0.
//    There are M ways to fold XW, yielding the following periods
//    P = P0 + i / M - 1
//    Where i ranges from 0 to M-1.  Summing all these values requires P
//    * M**2 = N**2 / P0 sums.  If M is a power of 2, the FFA eliminates
//    some of the redundant summing and requires only N log2 (N/P0)
//    sums.
//    Algorithm
//    ---------
//    The columns of XW are shifted and summed in a pairwise fashion.
//    - `FFAButterfly` : for a group of size nGroup, `FFAButterfly`
//      computes the amount pairwise combinations of rows and the amount
//      the second row is shifted.
//    - `FFAShiftAdd` : Adds the rows in the manner specified by
//      `FFAButterfly`
        
//    Parameters
//    ----------
//    XW : Wrapped array folded on P0.  shape(XW) = (P0,M) and M must be
//         a power of 2
    
//    Returns
//    -------
//    XWFS : XW Folded and Summed. 
//    References
//    ----------
//    [1] Staelin (1969)
//    [2] Kondratiev (2009)
//    """
QVector< QVector<float> > ffa_stages::FFA(const QVector<QVector<float> > &XW)
{
    mydisplay dd;
    QVector<float> res;
    int nRow = XW.size();
    int P0 = XW.first().size();
    qDebug() << "nRow = " << nRow << "P0 = " << P0;
    double nStageD = std::log2(static_cast<double>(nRow));
    Q_ASSERT_X(qAbs(nStageD - qRound(nStageD)) <=0.1e-15, "QVector<float> ffa_stages::FFA(const QVector<QVector<float> > &XW)",
               qPrintable(" nRow must be power of 2 but nStage = " + QString::number(nStageD)));
    //qDebug() << "nStage " << nStage;
    int nStage = static_cast<int>(nStageD);
    QVector< QVector<float> > XWFS = XW;

    for(int stage = 1; stage <= nStage; stage++)
    {
        //XWFS = FFAShiftAdd(XWFS.astype(np.float32),stage)
        XWFS = FFAShiftAdd(XWFS,stage);
        //qDebug() << "        stage " << stage;
        //dd.dispMatr(XWFS);
    }


    return XWFS;
}

// fast int power of 2
uint64_t ffa_stages::pow2(uint64_t i) const
{
    return std::uint64_t(1) << i;
}


//def FFAShiftAdd(cnp.ndarray[cnp.float32_t, ndim=2] XW0,
//                int stage):
//    """
//    FFA Shift and Add
//    Shuffle pairwise add the rows of the FFA data array corresponding
//    to stage

//    Parameters
//    ----------
//    XW0   : array
//    stage : The stage in the FFA.  An integer ranging from 1 to K
//            where 2**K = M

//    Returns
//    -------
//    XW    : Shifted and added array
//    Test Cases
//    ----------
//    >>> tfind.FFAShiftAdd(eye(4),1)
//    >>> array([[ 1.,  1.,  0.,  0.],
//               [ 2.,  0.,  0.,  0.],
//               [ 0.,  0.,  1.,  1.],
//               [ 0.,  0.,  2.,  0.]])
//    """
QVector< QVector<float> > ffa_stages::FFAShiftAdd(QVector< QVector<float> > &XW0, const int &stage)
{

    int nRow,nCol,nGroup,nRowGroup,iGroup,start,stop,lenz;
    QVector< QVector<float> > XW, XWtmp;
    nRow = XW0.size();
    nCol = XW0.first().size();
    nRowGroup = static_cast<int>(pow2(static_cast<uint64_t>(stage)));

    nGroup = nRow/nRowGroup;
    qDebug() << "nRowGroup test pow2 = " << nRowGroup << "nGroup = nRow/nRowGroup = " << nGroup;
    XW.fill(QVector<float>().fill(0.f, nCol), nRow);
    //qDebug() << "XW zeros " << XW;
    for(int iGroup = 0; iGroup < nGroup; iGroup++)
    {
        start = iGroup*nRowGroup;
        stop  = (iGroup+1)*nRowGroup;
        lenz = stop - start;
        //dd.dispMatr(XW0);
        qDebug()  << "start = " << start << "stop = " << stop;
        //dd.dispMatr(XW0.mid(start,1));
        XWtmp = XW.mid(start,lenz);

        FFAGroupShiftAdd(XW0.mid(start, lenz),XWtmp,nRowGroup,nCol);
        //dd.dispMatr(XW);
        //dd.dispMatr(XWtmp);
        for(int i=0; i < lenz; i++)
        {
            XW[start+i] = XWtmp.at(i);
        }

        //break;

        //FFAGroupShiftAdd(XW0[start:stop],XW[start:stop],nRowGroup,nCol);
    }
    return XW;
}

//def FFAGroupShiftAdd(cnp.ndarray[cnp.float32_t, ndim=2,mode='c'] group0,
//                     cnp.ndarray[cnp.float32_t, ndim=2,mode='c'] group,
//                     int nRowGroup,
//                     int nColGroup):

//    """
//    FFA Shift and Add
//    Add the rows of `group` to each other.

//    Parameters
//    ----------
//    group0 : Initial group before shuffling and adding.
//             shape(group0) = (M,P0) where M is a power of 2.
//    """

void ffa_stages::FFAGroupShiftAdd(const QVector< QVector<float> > &group0, QVector< QVector<float> > &group, int nRowGroup, int nColGroup)
{
    int iA,iB,Bs,jB;
    int nRowGroupOn2 = nRowGroup / 2; // # Half the group size

    //    # Grow group by the maximum shift value
    //    # Loop over rows in group
    for(int i = 0; i < nRowGroup; i++)
    {
        iA = i/2; //                # Row in the group that A is draw from
        iB = iA + nRowGroupOn2;   //# Row in group that B is drawn from
        Bs = (i + 1) / 2;
        //        # Loop over the columns in the group
        for(int j = 0; j < nColGroup; j++)
        {
            jB = (j + Bs + nColGroup) % nColGroup;
            //qDebug() << "i = " << i << "j=" << j;
            group[i][j] = group0.at(iA).at(j) + group0.at(iB).at(jB);
        }
    }
}
//    cdef int iRow,iCol,iA,iB,Bs,i,j,jB
//    cdef int nRowGroupOn2 = nRowGroup / 2 # Half the group size

//    # Grow group by the maximum shift value
//    # Loop over rows in group
//    for i in range(nRowGroup):
//        iA = i/2                 # Row in the group that A is draw from
//        iB = iA + nRowGroupOn2   # Row in group that B is drawn from
//        Bs = (i + 1) / 2
//        # Loop over the columns in the group
//        for j in range(nColGroup):
//            jB = (j + Bs + nColGroup) % nColGroup
//group[i,j] = group0[iA,j] + group0[iB,jB]
