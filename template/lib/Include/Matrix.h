// Matrix.h: interface for the CMatrix class.
// Part of KenTools library.
// Written by Kwan-Yee Kenneth Wong 2000.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#include <afxtempl.h>

template <class Type> 
class CMatrix : public CObject
{

public:
	// constructors
	CMatrix();
	CMatrix(int r, int c);
	CMatrix(int r, int c, Type v);
	CMatrix(const CMatrix<Type>& M);
	CMatrix(int Len);	// as vector

	// destructors
	virtual ~CMatrix();

	// initializations
		// if the size of the current allocated memory is the same as requested,
		// SetSize will only change the values of m_bRows and m_nCols and all the
		// values in the original matrix are retained

		// set the size of the matrix
	virtual BOOL SetSize(int Rows, int Cols);
		// set the size of the matrix and initialize every element to Val
	virtual BOOL SetSize(int Rows, int Cols, Type Val);
		// set the size of the vector
	virtual BOOL SetSize(int Len);
		// set every element to zero
	virtual void SetZeros();

	// retrievers
		// return number of rows
	inline int Rows() const {return m_nRows;}
		// return number of columns
	inline int Cols() const {return m_nCols;}
		// return the length of the vector
	inline int Len() const {return m_nRows*m_nCols;}
		// check if the matrix is valid
	inline BOOL isValid() const {return m_nRows!=0;}
		// return a pointer to a row
	inline const Type* Row(int r) const {return (m_pData + r*m_nCols);}
	inline Type* Row(int r) {return (m_pData + r*m_nCols);}
		// return a pointer to the internal buffer
	inline const Type* Buffer() const {return m_pData;}
	inline Type* Buffer() {return m_pData;}
		// retrieve the sub-matrix
	virtual BOOL SubMat(int r, int rr, int c, int cc, CMatrix<Type>& M) const;
	virtual CMatrix<Type> SubMat(int r, int rr, int c, int cc) const;
		// return a row as a separate vector
	virtual BOOL GetRow(int r, CMatrix<Type>& M) const;
	virtual CMatrix<Type> GetRow(int r) const;
		// return a column as a separate vector
	virtual BOOL GetCol(int c, CMatrix<Type>& M) const;
	virtual CMatrix<Type> GetCol(int c) const;
		// return the vector norm
	virtual double GetVecNorm() const;
		// return the homogeneous vector
	virtual CMatrix<Type> Homogeneous(Type lastdim = 1);

	// manipulators
		// check if the index is okay
	virtual BOOL isInside(int r, int c) const;
		// return the trace (sum of the diagonal elements)
	virtual Type Trace() const;
		// perform QR decomposition
	virtual BOOL QR(CMatrix<Type>& Q, CMatrix<Type>& R) const;
		// perform modified QR decomposition such that R is lower-triangular
	virtual BOOL QR2(CMatrix<Type>& Q, CMatrix<Type>& R) const;
		// perform Cholesky decomposition where L is lower-triangular and M = L.L^T 
		// only the upper triangular part of M will be considered and used in the decomposition
	virtual BOOL Chol(CMatrix<Type>& L) const;
		// perform singular value decomposition
	virtual BOOL SVD(CMatrix<Type>& U, CMatrix<Type>& D, CMatrix<Type>& V) const;
		// perform singular value decomposition with descending singular values in D
	virtual BOOL SVD2(CMatrix<Type>& U, CMatrix<Type>& D, CMatrix<Type>& V) const;
		// return the determinant of the matrix
	virtual double Det() const;
		// return the transpose
	virtual BOOL Transpose(CMatrix<Type>& M) const;
	virtual CMatrix<Type> Transpose() const;
		// return the inverse
	virtual BOOL Inverse(CMatrix<Type>& M) const;
	virtual CMatrix<Type> Inverse() const;
		// solve for X where M*X = B
	virtual BOOL Solve(const CMatrix<Type>& B, CMatrix<Type>& X) const;
	virtual CMatrix<Type> Solve(const CMatrix<Type>& B) const;
		// form a diagonal matrix assuming the values from matrix M row-wise,
		// failure will clean up this matrix
	virtual CMatrix<Type>& Diag(const CMatrix<Type>& M);
		// form a diagonal matrix assuming the values from array M with length l,
		// failure will clean up this matrix
	virtual CMatrix<Type>& Diag(const Type M[], int l);
		// form a diagonal matrix assuming the same value V,
		// the size of the matrix must be set beforehand
	virtual CMatrix<Type>& Diag(Type V);
		// form a vector from the diagonal entries of the matrix
	virtual CMatrix<Type> Diag() const;
		// set values of the submatrix to M,
		// failure will keep this matrix untouched
	virtual BOOL Set(int r, int c, const CMatrix<Type>& M);
		// set values of the submatrix to M(r1:r2,c1:c2)
		// failure will keep this matrix untouched
	virtual BOOL Set(int r, int c, const CMatrix<Type>& M, int r1, int r2, int c1, int c2);
		// set values of the submatrix to T,
		// failure will keep this matrix untouched
	virtual BOOL Set(int r, int rr, int c, int cc, const Type& T);
		// set the row of M,
		// failure will keep this matrix untouched
	virtual BOOL SetRow(int r, const CMatrix<Type>& Row);
		// set the column of M
		// failure will keep this matrix untouched
	virtual BOOL SetCol(int c, const CMatrix<Type>& Col);
		// swap rows r and rr
		// failure will keep this matrix untouched
	virtual BOOL SwapRows(int r, int rr);
		// swap columns c and cc
		// failure will keep this matrix untouched
	virtual BOOL SwapCols(int c, int cc);
		// swap values i and j
		// failure will keep this matrix untouched
	virtual BOOL Swap(int i, int j);
		// repeat matrix by rxc
	virtual BOOL RepMat(int r, int c, CMatrix<Type>& M) const;
	virtual CMatrix<Type> RepMat(int r, int c) const;
		// for vector cross product: M.Tx()*V = M x V, with both M and V are 3x1 vectors
	virtual BOOL Tx(CMatrix<Type>& M) const;
	virtual CMatrix<Type> Tx() const;
		// convolution and polynomial multiplication, treat all matrices as vectors
	virtual BOOL Conv(const CMatrix<Type>& M, CMatrix<Type>& Result) const;
	virtual CMatrix<Type> Conv(const CMatrix<Type>& M) const;
		// array multiplication (element-wise)
	virtual BOOL Times(const CMatrix<Type>& M, CMatrix<Type>& Result) const;
	virtual CMatrix<Type> Times(const CMatrix<Type>& M) const;
		// array division (element-wise)
	virtual BOOL Div(const CMatrix<Type>& M, CMatrix<Type>& Result) const;
	virtual CMatrix<Type> Div(const CMatrix<Type>& M) const;
		// calculate the mean for each row
	virtual BOOL RowMean(CMatrix<Type>& M) const;
	virtual CMatrix<Type> RowMean() const;
		// calculate the mean for each column
	virtual BOOL ColMean(CMatrix<Type>& M) const;
	virtual CMatrix<Type> ColMean() const;
		// calculate the mean for the whole matrix
	virtual double Mean() const;
		// clean up the memory and set size to zero
	virtual void CleanUp();

	// operators
		// access M(r,c) without range checking
	inline const Type& operator () (int r, int c) const {return m_pData[r*m_nCols+c];}
	inline Type& operator () (int r, int c) {return m_pData[r*m_nCols+c];}
		// access M(i) as a vector without range checking
	inline const Type& operator () (int i) const {return m_pData[i];}
	inline Type& operator () (int i) {return m_pData[i];}
		// return the submatrix (r:rr,c:cc)
	virtual CMatrix<Type> operator () (int r, int rr, int c, int cc) const;

		// negation
	virtual CMatrix<Type> operator - () const;

		// arithmetic operations with matrix M
	virtual CMatrix<Type> operator + (const CMatrix<Type>& M) const;
	virtual CMatrix<Type> operator - (const CMatrix<Type>& M) const;
	virtual CMatrix<Type> operator * (const CMatrix<Type>& M) const;
		
		// arithmetic operations with scalar T
	virtual CMatrix<Type> operator + (const Type& T) const;
	virtual CMatrix<Type> operator - (const Type& T) const;
	virtual CMatrix<Type> operator * (const Type& T) const;
	virtual CMatrix<Type> operator / (const Type& T) const;
		
		// assigment operations with matrix arithmetics (all prefix const eliminated by jerry)
	virtual CMatrix<Type>& operator  = (const CMatrix<Type>& M);
	virtual CMatrix<Type>& operator += (const CMatrix<Type>& M);
	virtual CMatrix<Type>& operator -= (const CMatrix<Type>& M);
	virtual CMatrix<Type>& operator *= (const CMatrix<Type>& M);
		
		// assigment operations with scalar arithmetics (all prefix const eliminated by jerry)
	virtual CMatrix<Type>& operator += (const Type& T);
	virtual CMatrix<Type>& operator -= (const Type& T);	
	virtual CMatrix<Type>& operator *= (const Type& T);
	virtual CMatrix<Type>& operator /= (const Type& T);

	// serializations
	virtual void SerializeTxt(CArchive& ar);
	virtual void Serialize(CArchive& ar);

	// debug purpose
	virtual void DumpMatrix();

protected:
	int   m_nRows;
	int   m_nCols;
	Type* m_pData;

	// helper functions
	double Pythag(double a,double b) const;
};

template <class Type>
CMatrix<Type> operator + (const Type& T, const CMatrix<Type>& M);

template <class Type>
CMatrix<Type> operator - (const Type& T, const CMatrix<Type>& M);

template <class Type>
CMatrix<Type> operator * (const Type& T, const CMatrix<Type>& M);

template <class Type>
CMatrix<Type> operator / (const Type& T, const CMatrix<Type>& M);

template <class Type>
Type dot(const CMatrix<Type>& A, const CMatrix<Type>& B);

template <class Type>
CMatrix<Type> cross(const CMatrix<Type>& A, const CMatrix<Type>& B);

// macro for instantiating different types of the matrix class
#define Instantiate_Matrix(Type)\
	template class CMatrix<Type>;\
	template CMatrix<Type> operator + (const Type& T, const CMatrix<Type>& M);\
	template CMatrix<Type> operator - (const Type& T, const CMatrix<Type>& M);\
	template CMatrix<Type> operator * (const Type& T, const CMatrix<Type>& M);\
	template CMatrix<Type> operator / (const Type& T, const CMatrix<Type>& M);\
	template Type dot(const CMatrix<Type>& A, const CMatrix<Type>& B);\
	template CMatrix<Type> cross(const CMatrix<Type>& A, const CMatrix<Type>& B);

typedef CTypedPtrList<CObList, CMatrix<BYTE>*> CByteMatrixList;
typedef CTypedPtrList<CObList, CMatrix<INT>*> CIntMatrixList;
typedef CTypedPtrList<CObList, CMatrix<FLOAT>*> CFloatMatrixList;
typedef CTypedPtrList<CObList, CMatrix<double>*> CDoubleMatrixList;
