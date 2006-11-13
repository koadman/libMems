#ifndef _Matrix_h_
#define _Matrix_h_

#include "gn/gnSetup.h"
#include <iostream>

template<class T>
class Matrix
{
public:
	Matrix();
	Matrix(unsigned nrows, unsigned ncols);
	// Throws a BadSize object if either size is zero
	class BadSize { };

	// Based on the Law Of The Big Three:
	~Matrix();	
	Matrix(const Matrix<T>& m);
	Matrix<T>& operator= (const Matrix<T>& m);
   	// Access methods to get the (i,j) element:	
	T& operator() (unsigned i, unsigned j);
	const T& operator() (unsigned i, unsigned j) const;
	// These throw a BoundsViolation object if i or j is too big
	class BoundsViolation { };
	// Support for initializing each matrix element to a value
	void init( const T& init_val );
	
	void print( ostream& os ) const;
	unsigned rows();
	unsigned cols();
protected:
	T* data_;
	unsigned nrows_, ncols_;
};
   
template<class T>
inline Matrix<T>::Matrix()
{
	data_ = NULL;
}

template<class T>
inline unsigned Matrix<T>::rows()
{
	return nrows_;
}

template<class T>
inline unsigned Matrix<T>::cols()
{
	return ncols_;
}

template<class T>
inline T& Matrix<T>::operator() (unsigned row, unsigned col)
{
	if (row >= nrows_ || col >= ncols_) 
		throw BoundsViolation();
	return data_[row*ncols_ + col];
}
   
template<class T>
inline const T& Matrix<T>::operator() (unsigned row, unsigned col) const
{
	if (row >= nrows_ || col >= ncols_) {
		cout << "debug me ";
		throw BoundsViolation();
	}
	return data_[row*ncols_ + col];
}
   
template<class T>
inline Matrix<T>::Matrix(unsigned nrows, unsigned ncols)
	: data_  (new T[nrows * ncols]),
	  nrows_ (nrows),
	  ncols_ (ncols)
{
	if (nrows == 0 || ncols == 0)
		throw BadSize();

}
template<class T>
inline Matrix<T>::Matrix(const Matrix<T>& m){
	*this = m;
}

template<class T>
inline Matrix<T>& Matrix<T>::operator=( const Matrix<T>& m )
{
	data_ = new T[m.nrows_ * m.ncols_];
	nrows_ = m.nrows_;
	ncols_ = m.ncols_;
	memcpy( data_, m.data_, nrows_ * ncols_ * sizeof( T ) );
	return *this;
}

template<class T>
inline Matrix<T>::~Matrix()
{
	if( data_ != NULL )
		delete[] data_;
}

template<class T>
inline void Matrix<T>::init( const T& init_val )
{
	for( unsigned rowI = 0; rowI < nrows_; rowI++ )
		for( unsigned colI = 0; colI < ncols_; colI++ )
			data_[ rowI * ncols_ + colI ] = init_val;
}

template<class T>
inline void Matrix<T>::print( ostream& os ) const{
	for( unsigned rowI = 0; rowI < nrows_; rowI++ ){
		for( unsigned colI = 0; colI < ncols_; colI++ ){
			if( colI > 0 )
				os << '\t';
			os << data_[ rowI * ncols_ + colI ];
		}
		os << endl;
	}
}

#endif // _Matrix_h_
