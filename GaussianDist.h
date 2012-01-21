/********************************************************************
	Year:      2009
	Author:    Johannes Feulner
*********************************************************************/

#ifndef GaussianDist_h__
#define GaussianDist_h__

#include <vector>

#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/algo/vnl_cholesky.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include "randnumgen.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include <iostream>

/*! \brief Class stores parameters of a normal distribution
	supports learning, evaluating, drawing and axis analysis
*/

template<int dimension>
class GaussianDist
{
public:
	typedef double ScalarType;
	typedef vnl_vector_fixed<ScalarType, dimension> DomainType;
	typedef vnl_matrix_fixed<ScalarType, dimension, dimension> CovarianceMatrixType;
	typedef std::vector<DomainType> DomainVectorType;

	ScalarType Eval(const DomainType& vec)const { return Eval(vec.data_block()); }
	ScalarType Eval(const ScalarType* vec)const;

	//! Draw a sample from the distribution
	void Draw(DomainType& v)const;

	//! Learn from samples
	void Learn(const DomainVectorType& trainingdata);

	const DomainType& GetMean()const { return m_mean; }
	const CovarianceMatrixType& GetCovarianceMatrix()const { return m_cov; }

	std::ostream& WriteAscii(std::ostream& s)const;
	std::istream& ReadAscii (std::istream& s);

	/*! \brief Computes PCA (eigensystem of covariance matrix)
		\param axis eigenvectors will be stored in columns. Same order as eigenvalues
		\param eigenvalues in order of descending absolute value, i.e. longest axis comes first
	*/
	void GetPrincipalComponents(CovarianceMatrixType& axis, DomainType& eigenvalues)const;

private:
	void ComputeInvCovLowTriAndPrefac();

	DomainType m_mean;
	CovarianceMatrixType m_cov;
	CovarianceMatrixType m_inv_cov;
	ScalarType m_prefac;
	CovarianceMatrixType m_lowertri;

	// used for sorting
	struct EigenvalueMappingWeakOrdering
	{
		int *mapping;
		vnl_symmetric_eigensystem< ScalarType > *eigs;


		bool operator ()(int a, int b)const
		{
			return fabs(eigs->get_eigenvalue(mapping[a])) > fabs(eigs->get_eigenvalue(mapping[b]));
		}
	};
};

template<int dimension>
std::ostream& GaussianDist<dimension>::WriteAscii(std::ostream &s) const
{
	s << m_mean << "\n" << m_cov;
	return s;
}

template<int dimension>
std::istream& GaussianDist<dimension>::ReadAscii(std::istream &s)
{
	s >> m_mean >> m_cov;
	if (s)
	{
		ComputeInvCovLowTriAndPrefac();
	}
	return s;
}



template<int dimension>
void GaussianDist<dimension>::GetPrincipalComponents(
	CovarianceMatrixType& axis, DomainType& eigenvalues)const
{
	vnl_symmetric_eigensystem< ScalarType > eigsys(m_cov);
	// sort eigensystem:
	int mapping[dimension];
	for (int i=0; i < dimension; ++i)
		mapping[i] = i;
	EigenvalueMappingWeakOrdering comp;
	comp.eigs = &eigsys;
	comp.mapping = mapping;
	std::sort(mapping, mapping+dimension, comp);
	for (int c=0; c < dimension; ++c)
	{
		int src = mapping[c];
		eigenvalues[c] = eigsys.get_eigenvalue(src);
		vnl_vector<ScalarType> vec = eigsys.get_eigenvector(src);
		for (int r=0; r < dimension; ++r)
			axis[r][c] = vec[r];
	}
}



template<int dimension>
void GaussianDist<dimension>::Learn(const DomainVectorType& trainingdata)
{
	m_mean.fill(0);
	m_cov.fill(0);

	// compute mean
	for (size_t i=0; i < trainingdata.size(); ++i)
	{
		const DomainType& v = trainingdata[i];
		m_mean += v;
	}
	ScalarType n = ScalarType(trainingdata.size());
	ScalarType fac = ScalarType(1)/n;
	m_mean *= fac;

	// make data zero-mean for better numerical behavior
	for (size_t i=0; i < trainingdata.size(); ++i)
	{
		DomainType v = trainingdata[i];
		v -= m_mean;
		for (int r=0; r < dimension; ++r)
		{
			for (int c=0; c < dimension; ++c)
			{
				m_cov[r][c] += v[r]*v[c];
			}
		}
	}
	m_cov *= fac;
	ComputeInvCovLowTriAndPrefac();
}



template<int dimension>
void GaussianDist<dimension>::ComputeInvCovLowTriAndPrefac()
{
	m_inv_cov = vnl_matrix_inverse<ScalarType>(m_cov);
	ScalarType det = abs(vnl_determinant(m_cov));
	m_prefac = ScalarType(1)/(pow(2.0*M_PI, double(dimension)*0.5)*sqrt(det));
	assert(m_prefac > 0);
	vnl_cholesky chol(m_cov);
	m_lowertri = chol.lower_triangle();
}

template<int dimension>
void GaussianDist<dimension>::Draw(DomainType& v)const
{
	DomainType gauss;
	for (int i=0; i < dimension; ++i)
		gauss[i] =  stdNormal();
	v = m_lowertri*gauss;
	v += m_mean;
}

template<int dimension>
typename GaussianDist<dimension>::ScalarType
GaussianDist<dimension>::Eval(const ScalarType* vec)const
{
	vnl_vector_ref<ScalarType> x(unsigned(dimension), const_cast<ScalarType*>(vec));
	DomainType dev = x-m_mean;
	ScalarType e = exp(-0.5*dot_product(dev, m_inv_cov*dev));
	assert(m_prefac > 0);
	return e*m_prefac;
}

template<int dimension>
std::ostream& operator << (std::ostream& s, const GaussianDist<dimension>& g)
{
	return g.WriteAscii(s);
}

template<int dimension>
std::istream& operator >> (std::istream& s, GaussianDist<dimension>& g)
{
	return g.ReadAscii(s);
}

#endif // GaussianDist_h__
