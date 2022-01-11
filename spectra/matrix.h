/*============================================================
 *
 *   Copyright 2011 Lex Kemper <kemper.af@gmail.com>
 *
 *============================================================
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>
#include <string.h>
#include <complex>
	using namespace std;

typedef complex<double> cdouble;
	

template <class dtype>
class matrix
{
	public:

	explicit matrix(const int dim1, const int dim2)
	{
		this-> dim1 = dim1;
		this-> dim2 = dim2;
		mat = new dtype[dim1*dim2];
		this->zero();

		memCount += dim1*dim2 * sizeof(dtype);
	}

	~matrix()
	{
		delete[] mat;
		memCount -= dim1*dim2 * sizeof(dtype);
	}

	inline int size() const
		{	return dim1*dim2; }
	inline int d1() const
		{   return dim1; }
	inline int d2() const
		{   return dim2; }
	
	inline dtype& at(const int &i, const int &j) const
	{
		assert(i < dim1 && j < dim2);
		return mat[i + dim1*j];
	}

	inline dtype& operator()(const int &i, const int &j) const
	{
		assert(i < dim1 && j < dim2);
		return mat[i + dim1*j];
	}
		
	inline dtype* ptr() const
		{	return mat; }

	inline dtype* col(const int &j) const
		{	assert(j < dim2);
			return mat + dim1 * j; }

	inline void set(const int &i, const int &j, const dtype &val)
		{	
			assert(i < dim1 && j < dim2);
			mat[i + dim1*j] = val; }

	inline void zero()
	{	
		for(int j=0; j < dim2; j++)
		{
			for(int i=0; i < dim1; i++)
				mat[i + dim1*j] = 0;
		}
	}

	void operator*= (const cdouble scalar)
	{
		for(int i=0; i < dim1*dim2; i++)
			mat[i] *= scalar;
	}

	void operator/= (const cdouble scalar)
	{
		for(int i=0; i < dim1*dim2; i++)
			mat[i] /= scalar;
	}

	void operator= (const matrix<dtype>& param)
	{	
		if(dim1 == param.dim1 && dim2 == param.dim2 )
			memcpy(mat, param.mat, dim1*dim2*sizeof(dtype));
	}


	void operator+= (const matrix<dtype>& rhs)
	{
		assert(rhs.dim1 == dim1 && rhs.dim2 == dim2);

		for(size_t i=0; i < dim1*dim2; i++)
			mat[i] += rhs.mat[i];
	}


	static long memCount;


	protected:
	int dim1, dim2;
	dtype* mat;


	// Allows for inheritance without being called
	protected:
	explicit matrix()
	{ ; }
};


template<typename cdouble> long matrix<cdouble>::memCount=0;

#endif
