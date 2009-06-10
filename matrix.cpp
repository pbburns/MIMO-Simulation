//
//  Operator/Method                          Description
//  ---------------                          -----------
//   operator ()   :   This function operator can be used as a
//                     two-dimensional subscript operator to get/set
//                     individual matrix elements.
//
//   operator !    :   This operator has been used to calculate inversion
//                     of matrix.
//
//   operator ~    :   This operator has been used to return transpose of
//                     a matrix.
//

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>
#define _NO_THROW               throw ()
#define _THROW_MATRIX_ERROR     throw (matrix_error)
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))
using namespace std;
namespace math {
class matrix_error : public logic_error
{
    public:
	matrix_error (const string& what_arg) : logic_error( what_arg) {}
};
#define REPORT_ERROR(ErrormMsg)  throw matrix_error( ErrormMsg);
#define MAT_TEMPLATE  template <class T>
#define matrixT  matrix<T>

MAT_TEMPLATE
class matrix
{
public:
   // Constructors
   matrix (const matrixT& m);
   matrix (size_t row = 6, size_t col = 6);

   // Destructor
   ~matrix ();

   // Assignment operators
   matrixT& operator = (const matrixT& m) _NO_THROW;

   // Value extraction method
   size_t RowNo () const { return _m->Row; }
   size_t ColNo () const { return _m->Col; }

   // Subscript operator
   T& operator () (size_t row, size_t col) _THROW_MATRIX_ERROR;
   T  operator () (size_t row, size_t col) const _THROW_MATRIX_ERROR;

   // Unary operators
   matrixT operator + () _NO_THROW { return *this; }
   matrixT operator - () _NO_THROW;

   // Combined assignment - calculation operators
   matrixT& operator += (const matrixT& m) _THROW_MATRIX_ERROR;
   matrixT& operator -= (const matrixT& m) _THROW_MATRIX_ERROR;
   matrixT& operator *= (const matrixT& m) _THROW_MATRIX_ERROR;
   matrixT& operator *= (const T& c) _NO_THROW;
   matrixT& operator /= (const T& c) _NO_THROW;
   matrixT& operator ^= (const size_t& pow) _THROW_MATRIX_ERROR;

   // Miscellaneous -methods
   void Unit (const size_t& row) _NO_THROW;
   void Unit () _NO_THROW;
   matrixT Inv () _THROW_MATRIX_ERROR;

private:
    struct base_mat
    {
	T **Val;
	size_t Row, Col, RowSiz, ColSiz;
	int Refcnt;

	base_mat (size_t row, size_t col, T** v)
	{
	    Row = row; RowSiz = row;
	    Col = col; ColSiz = col;
	    Refcnt = 1;

	    Val = new T* [row];
	    size_t rowlen = col * sizeof(T);

	    for (size_t i=0; i < row; i++)
	    {
			Val[i] = new T [col];
			if (v) memcpy( Val[i], v[i], rowlen);
	    }
	}
	~base_mat ()
	{
	    for (size_t i=0; i < RowSiz; i++)
			delete [] Val[i];
	    delete [] Val;
	}
    };
    base_mat *_m;

    void clone ();
    int pivot (size_t row);
};

// constructor
MAT_TEMPLATE inline
matrixT::matrix (size_t row, size_t col)
{
  _m = new base_mat( row, col, 0);
}

// copy constructor
MAT_TEMPLATE inline
matrixT::matrix (const matrixT& m)
{
    _m = m._m;
    _m->Refcnt++;
}

// Internal copy constructor
MAT_TEMPLATE inline void
matrixT::clone ()
{
    _m->Refcnt--;
    _m = new base_mat( _m->Row, _m->Col, _m->Val);
}

// destructor
MAT_TEMPLATE inline
matrixT::~matrix ()
{
   if (--_m->Refcnt == 0) delete _m;
}

// assignment operator
MAT_TEMPLATE inline matrixT&
matrixT::operator = (const matrixT& m) _NO_THROW
{
    m._m->Refcnt++;
    if (--_m->Refcnt == 0) delete _m;
    _m = m._m;
    return *this;
}


// subscript operator to get/set individual elements
MAT_TEMPLATE inline T&
matrixT::operator () (size_t row, size_t col) _THROW_MATRIX_ERROR
{
   if (row >= _m->Row || col >= _m->Col)
      REPORT_ERROR( "matrixT::operator(): Index out of range!");
   if (_m->Refcnt > 1) clone();
   return _m->Val[row][col];
}

// subscript operator to get/set individual elements
MAT_TEMPLATE inline T
matrixT::operator () (size_t row, size_t col) const _THROW_MATRIX_ERROR
{
   if (row >= _m->Row || col >= _m->Col)
      REPORT_ERROR( "matrixT::operator(): Index out of range!");
   return _m->Val[row][col];
}

// logical equal-to operator
MAT_TEMPLATE inline bool
operator == (const matrixT& m1, const matrixT& m2) _NO_THROW
{
   if (m1.RowNo() != m2.RowNo() || m1.ColNo() != m2.ColNo())
      return false;

   for (size_t i=0; i < m1.RowNo(); i++)
      for (size_t j=0; j < m1.ColNo(); j++)
	      if (m1(i,j) != m2(i,j))
	         return false;

   return true;
}

// logical no-equal-to operator
MAT_TEMPLATE inline bool
operator != (const matrixT& m1, const matrixT& m2) _NO_THROW
{
    return (m1 == m2) ? false : true;
}

// combined addition and assignment operator
MAT_TEMPLATE inline matrixT&
matrixT::operator += (const matrixT& m) _THROW_MATRIX_ERROR
{
   if (_m->Row != m._m->Row || _m->Col != m._m->Col)
      REPORT_ERROR( "matrixT::operator+= : Inconsistent matrix sizes in addition!");
   if (_m->Refcnt > 1) clone();
   for (size_t i=0; i < m._m->Row; i++)
      for (size_t j=0; j < m._m->Col; j++)
	 _m->Val[i][j] += m._m->Val[i][j];
   return *this;
}

// combined subtraction and assignment operator
MAT_TEMPLATE inline matrixT&
matrixT::operator -= (const matrixT& m) _THROW_MATRIX_ERROR
{
   if (_m->Row != m._m->Row || _m->Col != m._m->Col)
      REPORT_ERROR( "matrixT::operator-= : Inconsistent matrix sizes in subtraction!");
   if (_m->Refcnt > 1) clone();
   for (size_t i=0; i < m._m->Row; i++)
      for (size_t j=0; j < m._m->Col; j++)
	 _m->Val[i][j] -= m._m->Val[i][j];
   return *this;
}

// combined scalar multiplication and assignment operator
MAT_TEMPLATE inline matrixT&
matrixT::operator *= (const T& c) _NO_THROW
{
    if (_m->Refcnt > 1) clone();
    for (size_t i=0; i < _m->Row; i++)
	for (size_t j=0; j < _m->Col; j++)
	    _m->Val[i][j] *= c;
    return *this;
}

// combined matrix multiplication and assignment operator
MAT_TEMPLATE inline matrixT&
matrixT::operator *= (const matrixT& m) _THROW_MATRIX_ERROR
{
   if (_m->Col != m._m->Row)
      REPORT_ERROR( "matrixT::operator*= : Inconsistent matrix sizes in multiplication!");

   matrixT temp(_m->Row,m._m->Col);

   for (size_t i=0; i < _m->Row; i++)
      for (size_t j=0; j < m._m->Col; j++)
      {
         temp._m->Val[i][j] = T(0);
         for (size_t k=0; k < _m->Col; k++)
            temp._m->Val[i][j] += _m->Val[i][k] * m._m->Val[k][j];
      }
   *this = temp;

   return *this;
}

// combined scalar division and assignment operator
MAT_TEMPLATE inline matrixT&
matrixT::operator /= (const T& c) _NO_THROW
{
    if (_m->Refcnt > 1) clone();
    for (size_t i=0; i < _m->Row; i++)
	for (size_t j=0; j < _m->Col; j++)
	    _m->Val[i][j] /= c;

    return *this;
}

// combined power and assignment operator
MAT_TEMPLATE inline matrixT&
matrixT::operator ^= (const size_t& pow) _THROW_MATRIX_ERROR
{
	matrixT temp(*this);

	for (size_t i=2; i <= pow; i++)
      *this = *this * temp;

	return *this;
}

// unary negation operator
MAT_TEMPLATE inline matrixT
matrixT::operator - () _NO_THROW
{
   matrixT temp(_m->Row,_m->Col);

   for (size_t i=0; i < _m->Row; i++)
      for (size_t j=0; j < _m->Col; j++)
	 temp._m->Val[i][j] = - _m->Val[i][j];

   return temp;
}

// binary addition operator
MAT_TEMPLATE inline matrixT
operator + (const matrixT& m1, const matrixT& m2) _THROW_MATRIX_ERROR
{
   matrixT temp = m1;
   temp += m2;
   return temp;
}

// binary subtraction operator
MAT_TEMPLATE inline matrixT
operator - (const matrixT& m1, const matrixT& m2) _THROW_MATRIX_ERROR
{
   matrixT temp = m1;
   temp -= m2;
   return temp;
}

// binary scalar multiplication operator
MAT_TEMPLATE inline matrixT
operator * (const matrixT& m, const T& no) _NO_THROW
{
   matrixT temp = m;
   temp *= no;
   return temp;
}


// binary scalar multiplication operator
MAT_TEMPLATE inline matrixT
operator * (const T& no, const matrixT& m) _NO_THROW
{
   return (m * no);
}

// binary matrix multiplication operator
MAT_TEMPLATE inline matrixT
operator * (const matrixT& m1, const matrixT& m2) _THROW_MATRIX_ERROR
{
   matrixT temp = m1;
   temp *= m2;
   return temp;
}

// binary scalar division operator
MAT_TEMPLATE inline matrixT
operator / (const matrixT& m, const T& no) _NO_THROW
{
    return (m * (T(1) / no));
}


// binary scalar division operator
MAT_TEMPLATE inline matrixT
operator / (const T& no, const matrixT& m) _THROW_MATRIX_ERROR
{
    return (!m * no);
}

// binary matrix division operator
MAT_TEMPLATE inline matrixT
operator / (const matrixT& m1, const matrixT& m2) _THROW_MATRIX_ERROR
{
    return (m1 * !m2);
}

// binary power operator
MAT_TEMPLATE inline matrixT
operator ^ (const matrixT& m, const size_t& pow) _THROW_MATRIX_ERROR
{
   matrixT temp = m;
   temp ^= pow;
   return temp;
}

// unary transpose operator
MAT_TEMPLATE inline matrixT
operator ~ (const matrixT& m) _NO_THROW
{
   matrixT temp(m.ColNo(),m.RowNo());

   for (size_t i=0; i < m.RowNo(); i++)
      for (size_t j=0; j < m.ColNo(); j++)
      {
         T x = m(i,j);
	      temp(j,i) = x;
      }
   return temp;
}

// unary inversion operator
MAT_TEMPLATE inline matrixT
operator ! (const matrixT m) _THROW_MATRIX_ERROR
{
   matrixT temp = m;
   return temp.Inv();
}

// inversion function
MAT_TEMPLATE inline matrixT
matrixT::Inv () _THROW_MATRIX_ERROR
{
   size_t i,j,k;
   T a1,a2,*rowptr;

   if (_m->Row != _m->Col)
      REPORT_ERROR( "matrixT::operator!: Inversion of a non-square matrix");

   matrixT temp(_m->Row,_m->Col);
   if (_m->Refcnt > 1) clone();
   temp.Unit();
   for (k=0; k < _m->Row; k++)
   {
      int indx = pivot(k);
      if (indx == -1) 
	  {
		  temp(0,0) = -1010;
		  return temp;
	  }
      if (indx != 0)
      {
	      rowptr = temp._m->Val[k];
	      temp._m->Val[k] = temp._m->Val[indx];
	      temp._m->Val[indx] = rowptr;
      }
      a1 = _m->Val[k][k];
      for (j=0; j < _m->Row; j++)
      {
	      _m->Val[k][j] /= a1;
	      temp._m->Val[k][j] /= a1;
      }
      for (i=0; i < _m->Row; i++)
	   if (i != k)
	   {
	      a2 = _m->Val[i][k];
	      for (j=0; j < _m->Row; j++)
	      {
	         _m->Val[i][j] -= a2 * _m->Val[k][j];
	         temp._m->Val[i][j] -= a2 * temp._m->Val[k][j];
	      }
	   }
   }
   return temp;
}

// set this matrix to unity
MAT_TEMPLATE inline void
matrixT::Unit (const size_t& row) _NO_THROW
{
    if (row != _m->Row || row != _m->Col)
	realloc( row, row);
	
    if (_m->Refcnt > 1) 
	clone();

    for (size_t i=0; i < _m->Row; i++)
	for (size_t j=0; j < _m->Col; j++)
	    _m->Val[i][j] = i == j ? T(1) : T(0);
    return;
}

// set this matrix to unity
MAT_TEMPLATE inline void
matrixT::Unit () _NO_THROW
{
    if (_m->Refcnt > 1) clone();   
    size_t row = min(_m->Row,_m->Col);
    _m->Row = _m->Col = row;

    for (size_t i=0; i < _m->Row; i++)
	for (size_t j=0; j < _m->Col; j++)
	    _m->Val[i][j] = i == j ? T(1) : T(0);
    return;
}

// private partial pivoting method
MAT_TEMPLATE inline int
matrixT::pivot (size_t row)
{
  int k = int(row);
  double amax,temp;

  amax = -1;
  for (size_t i=row; i < _m->Row; i++)
    if ( (temp = abs( _m->Val[i][row])) > amax && temp != 0.0)
     {
       amax = temp;
       k = i;
     }
  if (_m->Val[k][row] == T(0))
     return -1;
  if (k != int(row))
  {
     T* rowptr = _m->Val[k];
     _m->Val[k] = _m->Val[row];
     _m->Val[row] = rowptr;
     return k;
  }
  return 0;
}

}


