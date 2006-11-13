#ifndef _NumericMatrix_h_
#define _NumericMatrix_h_

#include "Matrix.h"

template<class T>  // See section on templates for more
class NumericMatrix : public Matrix<T>
{
public:
	NumericMatrix(){};
	NumericMatrix(unsigned nrows, unsigned ncols);
   
      // Based on the Law Of The Big Three:
	~NumericMatrix();
	NumericMatrix(const NumericMatrix<T>& m);
	NumericMatrix<T>& operator= (const NumericMatrix<T>& m);
	
	// define some arithmetic operators 
	NumericMatrix<T>& operator+= (const NumericMatrix<T>& m);
	NumericMatrix<T>& operator-= (const NumericMatrix<T>& m);
	// not implemented
	NumericMatrix<T>& operator*= (const NumericMatrix<T>& m);
	NumericMatrix<T>& operator*= (const T& m);
	NumericMatrix<T>& operator/= (const T& m);

	// the following 5 are not implemented
	NumericMatrix<T>& operator+ (const NumericMatrix<T>& m ) const;
	const NumericMatrix<T>& operator- (const NumericMatrix<T>& m ) const;
	const NumericMatrix<T>& operator* (const NumericMatrix<T>& m ) const;
	const NumericMatrix<T>& operator* (const T& n) const;
	const NumericMatrix<T>& operator/ (const T& n) const;

};
   
template<class T>
inline NumericMatrix<T>::NumericMatrix(unsigned nrows, unsigned ncols)
	: Matrix<T>( nrows, ncols )
{
}
   
template<class T>
inline NumericMatrix<T>::NumericMatrix(const NumericMatrix<T>& m){
	*this = m;
}

template<class T>
inline NumericMatrix<T>& NumericMatrix<T>::operator= (const NumericMatrix<T>& m)
{
	Matrix<T>::operator=( m );
	return *this;
}

template<class T>
inline NumericMatrix<T>::~NumericMatrix()
{
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator+= (const NumericMatrix<T>& m){
	// make sure matrix dimensions agree
	if (this->nrows_ != m.nrows_ || this->ncols_ != m.ncols_)
		throw BadSize();

	// do the arithmetic on each matrix entry
	for(unsigned i = 0; i < nrows_ * ncols_; i++ )
		this->data_[ i ] += m.data_[ i ];
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator-= (const NumericMatrix<T>& m){
	// make sure matrix dimensions agree
	if (this->nrows_ != m.nrows_ || this->ncols_ != m.ncols_)
		throw BadSize();

	// do the arithmetic on each matrix entry
	for(unsigned i = 0; i < nrows_ * ncols_; i++ )
		this->data_[ i ] -= m.data_[ i ];
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator*= (const NumericMatrix<T>& m){
	// make sure matrix dimensions agree
	if (this->ncols_ != m.nrows_)
		throw BadSize();
	// do a matrix multiply
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator*= (const T& m){
	// do the arithmetic on each matrix entry
	for(unsigned i = 0; i < nrows_ * ncols_; i++ )
		this->data_[ i ] *= m;
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator/= (const T& m){
	// do the arithmetic on each matrix entry
	for(unsigned i = 0; i < nrows_ * ncols_; i++ )
		this->data_[ i ] /= m;
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator+ (const NumericMatrix<T>& m) const {

}
template<class T>
inline
const NumericMatrix<T>& NumericMatrix<T>::operator- (const NumericMatrix<T>& m) const {

}
template<class T>
inline
const NumericMatrix<T>& NumericMatrix<T>::operator* (const NumericMatrix<T>& m) const {

}
template<class T>
inline
const NumericMatrix<T>& NumericMatrix<T>::operator* (const T& n) const {

}
template<class T>
inline
const NumericMatrix<T>& NumericMatrix<T>::operator/ (const T& n) const {

}


#endif // _NumericMatrix_h_
