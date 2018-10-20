#ifndef MYQVECTOR_H
#define MYQVECTOR_H

#include <QObject>
#include <QVector>
#include <QDebug>

template<class A>
class MyQVector
{

public:
    explicit MyQVector(int nRow=1, int nCol=1, A fill_val = static_cast<A>(0))
    {
        nnRow = 1;
        nnCol = 1;
        flag_isMatrix=false;

        if(nRow>0)
        {
            nnRow = nRow;
        }
        if(nCol>0)
        {
            nnCol = nCol;
            if(nnCol > 1)
                flag_isMatrix=true;
        }
        matrix.fill(fill_val,nnRow*nnCol);


//        if(flag_isMatrix)
//            matrix.resize(nnRow*nnCol);
//        else
//            matrix.resize(nnRow);


    }

//    MyQVector<A> reshape(const MyQVector<A> &vec, int row, int col)
//    {
//        MyQVector<A> res(vec);
//        res.reshape(row,col);
//        return res;
//    }

    //===================================================================================================
    // reshape MyQVector matrix with new row and column
    //===================================================================================================
    bool reshape(int row, int col)
    {
        int newsize = row*col;
        if(newsize > nnRow*nnCol)
        {
            QVector<A> tmp(newsize);
            int k = 0;
            for (const A &v : qAsConst(matrix))
            {
                tmp[k]=v;
                k++;
            }
            tmp.swap(matrix);
            nnRow = qAbs(row);
            nnCol = qAbs(col);
            if(nnCol > 1)
                flag_isMatrix=true;
            else
                flag_isMatrix=false;
            return true;
        }
        if(newsize==nnRow*nnCol && nnCol!=col)
        {
            nnRow = qAbs(row);
            nnCol = qAbs(col);
            if(nnCol > 1)
                flag_isMatrix=true;
            else
                flag_isMatrix=false;
            return true;
        }
        return false;
    }

    //===================================================================================================
    // get column vector from MyQvector matrix by index
    //===================================================================================================
    QVector<A> getColumn(int col=0)
    {
        QVector<A> res(nnRow);
        if(nnRow>1 && nnCol>1)
        {
            if(col<0 || col>=nnCol)
                col=nnCol-1;
            for(int row=0; row< nnRow; row++)
            {
                res[row] = matrix.at(col + row * nnCol);
            }
        }
        return res;
    }

    //===================================================================================================
    // get row vector from MyQvector matrix by index
    //===================================================================================================
    QVector<A> getRow(int row=0)
    {
        QVector<A> res(nnCol);
        if(nnRow>1 && nnCol>1)
        {
            if(row<0 || row>=nnRow)
                row=nnRow-1;
            res = matrix.mid(row * nnCol, nnCol);
        }
        return res;
    }


    //===================================================================================================
    // get sub vector from MyQvector
    //===================================================================================================
    QVector<A> getVector(int i_start = 0, int i_end=-1)
    {

        int newsize = matrix.size();

        if(i_start==0 && i_end<0)
        {
            return matrix;
        }
        else if (i_start>=0 && i_end > i_start && i_end <=newsize)
        {
            newsize=i_end-i_start;

            return matrix.mid(i_start,newsize);
        }
        else if (i_start>=0 && i_end < 0 && newsize+i_end >0 )
        {
            newsize= newsize + i_end - i_start+1;
            return matrix.mid(i_start,newsize);
        }
        return matrix;
    }

    bool isMatrix()
    {
        return flag_isMatrix;
    }

    //===================================================================================================
    // fill matrix with values
    //===================================================================================================
    void fill(const A &fill_val, int nRow = 1, int nCol=1)
    {   
        if(nRow<=0)
            nRow=1;
        if(nCol<=0)
            nCol=1;

        nnRow=nRow;
        nnCol=nCol;

        if(nnCol>1)
           flag_isMatrix=true;
        else
           flag_isMatrix=false;

        matrix.fill(fill_val,nnRow*nnCol);
    }

    //===================================================================================================
    // maximum for vector
    //===================================================================================================
    A maxVec(const QVector<A> &vect) const
    {
        return *std::max_element(vect.constBegin(), vect.constEnd());
    }

    //===================================================================================================
    // maximum for Matrix row or column; axis = 0 or 1
    //===================================================================================================
    MyQVector<A> maxAxis(int axis = 0)
    {
        if(!flag_isMatrix || nnRow==1 || nnCol == 1)
        {
            MyQVector<A> mxtmp;
            mxtmp.fill(maxVec(matrix),1,1);
            return mxtmp;
        }
        else
        {
            QVector<A> tmp;
            MyQVector<A> resMAx;

            switch (axis) {
            case 0:
                resMAx.fill(0,nnCol);

                for(int col=0; col < nnCol; col++)
                {
                        tmp = getColumn(col);
                        resMAx[col] = maxVec(tmp);
                }
                break;
            case 1:
                resMAx.fill(0,nnRow);
                for(int row=0; row < nnRow; row++)
                {
                    tmp = matrix.mid(row*nnCol,nnCol);
                    resMAx[row] = maxVec(tmp);
                }
                break;
            default:
                break;
            }
            return resMAx;
        }
    }

    //===================================================================================================
    // minimum for vector
    //===================================================================================================
    A minVec(const QVector<A> &vect) const
    {
        return *std::min_element(vect.constBegin(), vect.constEnd());
    }

    //===================================================================================================
    // minimum for Matrix row or column; axis = 0 or 1
    //===================================================================================================
    MyQVector<A> minAxis(int axis = 0)
    {
        if(!flag_isMatrix || nnRow==1 || nnCol == 1)
        {
            MyQVector<A> mntmp;
            mntmp.fill(minVec(matrix),1,1);
            return mntmp;
        }
        else
        {
            QVector<A> tmp;
            MyQVector<A> resMin;

            switch (axis) {
            case 0:
                resMin.fill(0,nnCol);
                for(int col=0; col < nnCol; col++)
                {
                        tmp = getColumn(col);
                        resMin[col] = minVec(tmp);
                }
                break;
            case 1:
                resMin.fill(0,nnRow);
                for(int row=0; row < nnRow; row++)
                {
                    tmp = matrix.mid(row*nnCol,nnCol);
                    resMin[row] = minVec(tmp);
                }
                break;
            default:
                break;
            }
            return resMin;
        }
    }
//===================================================================================================
// http://www.cplusplus.com/forum/beginner/232940/    median for odd and even time series
//===================================================================================================
    A median( const QVector<A> &v ) const
    {
        QVector<A> seq = v;

        const auto n = seq.size() ;

        std::nth_element( seq.begin(), seq.begin() + n/2, seq.end() );
        const A m1 = seq[n/2] ;

        if( n%2 == 1 ) return m1 ; // if n is odd

        // even n
        std::nth_element( seq.begin(), seq.begin() + (n-1)/2, seq.end() );
        const A d2 = static_cast<A>(2);
        return ( m1 + seq[ (n-1)/2 ] ) / d2 ;
    }

 //===================================================================================================
 // median for Matrix row or column; axis = 0 or 1
 //===================================================================================================
    MyQVector<A>  medianAxis(int axis = 0)
    {
        if(!flag_isMatrix || nnRow==1 || nnCol == 1)
        {
            MyQVector<A> medtmp;
            medtmp.fill(median(matrix),1,1);
            return medtmp;
        }
        else
        {
            QVector<A> tmp;
            MyQVector<A> medtmp;

            switch (axis) {
            case 0:  // column
                medtmp.fill(0,nnCol);
                for(int col=0; col < nnCol; col++)
                {
                        tmp = getColumn(col);
                        medtmp[col] = median(tmp);
                }
                break;
            case 1: // row
                medtmp.fill(0,nnRow);
                for(int row=0; row < nnRow; row++)
                {
                    tmp = matrix.mid(row*nnCol,nnCol);
                    medtmp[row] = median(tmp);
                }
                break;
            default:
                break;
            }
            return medtmp;
        }
    }

    //===================================================================================================
    // transpon Matrix
    //===================================================================================================
    void transpon()
    {
        int tmprow = nnRow;
        if(nnRow==1 || nnCol==1)
        {
            nnRow = nnCol;
            nnCol = tmprow;

            if(nnCol<=1)
                flag_isMatrix = false;
            else
                flag_isMatrix = true;
        }
        else if(nnRow > 1 && nnCol > 1)
        {
            int k=0;
            QVector<A> tmp(nnRow*nnCol);
            for(int col=0; col< nnCol; col++)
            {
                for(int row=0; row< nnRow; row++)
                {

                    tmp[k] = matrix.at(col + row * nnCol);
                    k++;
                }
            }
            nnRow = nnCol;
            nnCol = tmprow;
            qDebug() << tmp;
            tmp.swap(matrix);
            flag_isMatrix = true;
        }
    }

    //===================================================================================================
    // generate range vector from 0 to num integer value
    //===================================================================================================
    void arange(const int &num)
    {
        if(flag_isMatrix)
        {
            flag_isMatrix=false;
            nnCol=1;
        }
        int cnum = num;
        int sign = 1;
        if(num < 0)
        {
            cnum = (-1)*num;
            sign = -1;
        }
        matrix.resize(cnum);
        nnRow = cnum;
        for(int i=0; i<cnum; i++)
        {
            matrix[i] = static_cast<A>(i*sign);
        }
    }

    //===================================================================================================
    // Overload generate range vector from n_start to n_end integer value with step delta
    //===================================================================================================
    void arange(const int &n_start, const int &n_end, int delta = 1)
    {

        int num = n_end - n_start + delta;

        if( num > 0 && delta > 0)
        {
            num = num/delta;

            if(flag_isMatrix)
            {
                flag_isMatrix=false;
                nnCol=1;
            }
            matrix.resize(num);
            nnRow = num;

            for(int i=0; i<num; i++)
            {
                matrix[i] = static_cast<A>(n_start) + static_cast<A>(delta)*static_cast<A>(i);
            }
        }
    }

    //===================================================================================================
    // size of matrix row or column
    //===================================================================================================
    inline int size(const int &naxis) const
    {
        switch (naxis) {
        case 0:
            return nnRow;
            break;
        case 1:
            return nnCol;
            break;
        default:
            return 0;
            break;
        }
    }

    //===================================================================================================
    // size of matrix as a vector size
    //===================================================================================================
    int msize() const
    {
        return nnCol*nnRow;
    }

    //===================================================================================================
    // size of matrix as a QPair<row,col>
    //===================================================================================================
    QPair<int, int> sizeQ()
    {
        return QPair<int, int>(nnRow, nnCol);
    }
    //===================================================================================================
    // analog of pythom numpy masked_where
    //===================================================================================================
    void masked_where(const A condition, bool maskVal = true)
    {
        matrix_mask.fill(!maskVal,nnCol*nnRow);
        int k=0;
        for (const A &v : qAsConst(matrix))
        {
            if(qFuzzyCompare(v,condition))
            {
                matrix_mask[k] = maskVal;
            }
            k++;
        }
    }

    //===================================================================================================
    // Overload operator =  For assignment QVector<A> to MyQVector object
    //===================================================================================================
    inline void operator = (const QVector<A> &Vect ) {
        if(Vect.size() > 0)
        {
             matrix = Vect;
             nnCol=1;
             nnRow = Vect.size();
             flag_isMatrix=false;
         }
        else{
            qDebug() << "inline void operator = (const QVector<A> &Vect ) " << "Error Size of rvector <=0 ";
        }
    }

    //===================================================================================================
    // Overload operator =  For assignment one MyQVector object to another
    //===================================================================================================
    inline void operator = (MyQVector<A> &Vect )
    {
             nnRow=Vect.size(0);
             nnCol = Vect.size(1);
             matrix = Vect.getVector();
             flag_isMatrix=Vect.isMatrix();
    }

    //===================================================================================================
    // Overload operator =
    //===================================================================================================
    MyQVector<A> &operator=(const MyQVector<A>  &v)
    {
         QVector<A> tmp(v.matrix);
         tmp.swap(this->matrix);
         this->flag_isMatrix = v.flag_isMatrix;
         this->nnCol = v.nnCol;
         this->nnRow = v.nnRow;
         return *this;
    }

    //===================================================================================================
    // Overload operator []  For assignment to matrix element
    //===================================================================================================
    A& operator[](int n) {
        int m_size = matrix.size();
        if( n >= m_size )
            n = m_size - 1;
        //qDebug() << "A& operator[](int n)  m_size = matrix.size() = " << m_size << " n= " << n;
        return (matrix.begin()[n]);
    }

    //===================================================================================================
    // Overload operator []  For get matrix element value
    //===================================================================================================
    const A& operator[](int n) const
    {
        int m_size = matrix.size();
        if( n >= m_size )
            n = m_size - 1;
        return matrix.at(n);

    }

    //===================================================================================================
    // Overload at() For get matrix element value as vector index
    //===================================================================================================
    A at(int n) const
    {
        int m_size = matrix.size();
        if( n >= m_size )
            n = m_size - 1;
        return matrix.at(n);
    }

    //===================================================================================================
    // Overload at() For get matrix element value as matrix nxm indexes
    //===================================================================================================
    A at(int row, int col) const
    {
        int m_size = matrix.size();
        int n_index = nnCol * row + col;
        if( n_index >= m_size )
            n_index = m_size - 1;
        return matrix.at(n_index);
    }

    //===================================================================================================
    // Overload operator ()  For assignment matrix element value by matrix indexes
    //===================================================================================================
    A& operator()(int row, int col)
    {
        int n_index = nnCol * row + col;
        int m_size = matrix.size();
        if(n_index >= m_size)
            n_index = m_size - 1;
        return (matrix.begin()[n_index]);
    }

    //===================================================================================================
    // Overload operator ()  For get matrix element value by matrix indexes
    //===================================================================================================
    const A& operator()(int row, int col) const
    {
        int n_index = nnCol * row + col;
        int m_size = matrix.size();
        if(n_index >= m_size)
            n_index = m_size - 1;
        return matrix.at(n_index);
    }

    //===================================================================================================
    // Overload operator + For adds two MyVectors object as  2d arrays
    //===================================================================================================
    inline MyQVector<A> operator+(const MyQVector<A> &l) const
    { MyQVector<A> n = *this; n += l; return n; }

    //===================================================================================================
    // Overload operator -
    //===================================================================================================
    inline MyQVector<A> operator-(const MyQVector<A> &l) const
    { MyQVector<A> n = *this; n -= l; return n; }

    //===================================================================================================
    // Overload operator +  to add value to each element of matrix
    //===================================================================================================
    inline MyQVector<A> operator+(const A &l) const
    { MyQVector<A> n = *this; n += l; return n; }

    //===================================================================================================
    // Overload operator -  to del value from each element of matrix
    //===================================================================================================
    inline MyQVector<A> operator-(const A &l) const
    { MyQVector<A> n = *this; n -= l; return n; }

    //===================================================================================================
    // Overload operator *  to multiple two matrixes
    //===================================================================================================
    inline MyQVector<A> operator*(const MyQVector<A> &l) const
    { MyQVector<A> n = *this; n *= l; return n; }

    //===================================================================================================
    // Overload operator *  to multiple matrix by value
    //===================================================================================================
    inline MyQVector<A> operator*(const A &l) const
    { MyQVector<A> n = *this; n *= l; return n; }

    //===================================================================================================
    // Overload operator +=
    //===================================================================================================
    inline MyQVector<A> &operator+=(const MyQVector<A> &v)
    {
        int nRow1 =  this->size(0);
        int nRow2 = v.size(0);
        int nCol1 = this->size(1);
        int nCol2 = v.size(1);
        if(nRow1 != nRow2 || nCol1 != nCol2)
        {
            qDebug() << "FATAL vector::operator+=(const vector &) size mismatch: " <<
                        nRow1 << " " << " != " << nRow2 << "\n";
            exit(1);
        }
        for(int i = 0; i < nRow1*nCol1; i++) this->operator[](i) += v.at(i);

        return *this;
    }

    //===================================================================================================
    // Overload operator +=
    //===================================================================================================
    inline MyQVector<A> &operator+=(const A &v)
    {
        int nRow1 =  this->size(0);

        int nCol1 = this->size(1);

        if(nRow1 <=0 || nCol1 <=0 )
        {
            qDebug() << "FATAL vector::operator+=(const vector &) size mismatch: " <<
                        nRow1 << " x " << nCol1 << "\n";
            exit(1);
        }
        for(int i = 0; i < nRow1*nCol1; i++) this->operator[](i) += v;

        return *this;
    }

    //===================================================================================================
    // Overload operator -=
    //===================================================================================================
    inline MyQVector<A> &operator-=(const A &v)
    {
        int nRow1 =  this->size(0);

        int nCol1 = this->size(1);

        if(nRow1 <=0 || nCol1 <=0 )
        {
            qDebug() << "FATAL vector::operator+=(const vector &) size mismatch: " <<
                        nRow1 << " x " << nCol1 << "\n";
            exit(1);
        }
        for(int i = 0; i < nRow1*nCol1; i++) this->operator[](i) -= v;

        return *this;
    }

    //===================================================================================================
    // Overload operator -=
    //===================================================================================================
    inline MyQVector<A> &operator-=(const MyQVector<A> &v)
    {
        int nRow1 =  this->size(0);
        int nRow2 = v.size(0);
        int nCol1 = this->size(1);
        int nCol2 = v.size(1);
        if(nRow1 != nRow2 || nCol1 != nCol2)
        {
            qDebug() << "FATAL vector::operator-=(const vector &) size mismatch: " <<
                        nRow1 << " " << " != " << nRow2 << "\n";
            exit(1);
        }
        for(int i = 0; i < nRow1*nCol1; i++) this->operator[](i) -= v.at(i);

        return *this;
    }

    //===================================================================================================
    // Overload operator *=
    //===================================================================================================
    inline MyQVector<A> &operator*=(const MyQVector<A> &v)
    {
        int nRow1 =  this->size(0);
        int nRow2 = v.size(0);
        int nCol1 = this->size(1);
        int nCol2 = v.size(1);
        if(nRow1 != nRow2 || nCol1 != nCol2)
        {
            qDebug() << "FATAL vector::operator*=(const vector &) size mismatch: " <<
                        nRow1 << " " << " != " << nRow2 << "\n";
            exit(1);
        }
        for(int i = 0; i < nRow1*nCol1; i++) this->operator[](i) *= v.at(i);

        return *this;
    }

    //===================================================================================================
    // Overload operator *=
    //===================================================================================================
    inline MyQVector<A> &operator*=(const A &v)
    {
        int nRow1 =  this->size(0);

        int nCol1 = this->size(1);

        if(nRow1 <=0 || nCol1 <= 0)
        {
            qDebug() << "FATAL vector::operator*=(const vector &) size mismatch: " <<
                        nRow1 << " x "  << nCol1 << "\n";
            exit(1);
        }
        for(int i = 0; i < nRow1*nCol1; i++) this->operator[](i) *= v;

        return *this;
    }


    //===================================================================================================
    // For debug to show matrix
    //===================================================================================================
    void show()
    {
        QString strOut = "[";
        if(flag_isMatrix || nnCol>1)
        {
            qDebug() << "Matrix size of " << nnRow << "x" << nnCol;

            int k=0;
            for (const A &val : qAsConst(matrix))
            //Q_FOREACH(const A &val, matrix)
            {
                if(k==0)
                    strOut.append("[");
                else if(k>=nnCol)
                {
                    strOut.append("]");
                    qDebug() << strOut;
                    strOut = "[";
                    k=0;
                }
                strOut.append(QString::number(val));
                strOut.append(" ");
                k++;

            }
            strOut.append("]]");
            qDebug() << strOut;

        }
        else
        {
            qDebug() << "Vector size of " <<  nnRow;

            for (const A &val : qAsConst(matrix))
            {
                strOut.append(QString::number(val));
                strOut.append(" ");
            }
            strOut.append("]");
            qDebug() << strOut;

        }


    }

    //===================================================================================================
    // For debug to show matrix mask
    //===================================================================================================
    void show_masked()
    {
        QString strOut = "[";
        if(flag_isMatrix || nnCol>1)
        {
            qDebug() << "Mask Matrix size of " << nnRow << "x" << nnCol;

            int k=0;
            for (const bool &val : qAsConst(matrix_mask))
            {
                if(k==0)
                    strOut.append("[");
                else if(k>=nnCol)
                {
                    strOut.append("]");
                    qDebug() << strOut;
                    strOut = "[";
                    k=0;
                }
                if(val)
                    strOut.append("True");
                else
                    strOut.append("False");
                strOut.append(" ");
                k++;
            }
            strOut.append("]]");
            qDebug() << strOut;

        }
        else
        {
            qDebug() << "Mask Vector size of " <<  nnRow;

            for (const bool &val : qAsConst(matrix_mask))
            {
                if(val)
                    strOut.append("True");
                else
                    strOut.append("False");
                strOut.append(" ");
            }
            strOut.append("]");
            qDebug() << strOut;
        }
    }

public slots:

private:
    bool flag_isMatrix;
    int nnRow;
    int nnCol;
    QVector<A> matrix;
    QVector<bool> matrix_mask;


};

#endif // MYQVECTOR_H
