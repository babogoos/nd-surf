/********************************************************************
	Year:      2009
	Author:    Johannes Feulner
	Siemens AG, CT SE SCR2
*********************************************************************/

#ifndef Surf_h__
#define Surf_h__

#include <vector>
#include <assert.h>
#include "IntegralImage.h"
#include "GaussianDist.h"
#include "Ipoint.h"

namespace nsurf

{

class SurfBase
{
public:
	enum RotInvarMethod
	{
		upright,
		PCA,
		MeanGradPCA
	};
};

/*!	Computes SURF descriptor at arbitrary dimension
*/

template<typename TPixel, int dim>
class Surf : public SurfBase
{
public:
	typedef TPixel PixelType;
	typedef typename IntegralType<PixelType>::PixelType IntPixelType;
	static const unsigned int Dimension = dim;

	typedef itk::Image<PixelType, dim> ImageType;
	typedef IntegralImage<PixelType, dim> IntegralImageType;
	typedef itk::Index<dim> IndexType;
	typedef itk::ImageRegion<dim> RegionType;

	//! Standard Constructor (img is an integral image)
	Surf(const IntegralImageType& img, std::vector<Ipoint<dim> > *ipts = 0);

	//! Describe all features in the supplied vector
	void getDescriptors(RotInvarMethod rotinv)const;

	void getDescriptor(Ipoint<dim>& ipoint, RotInvarMethod rotinv)const;
	void getOrientation(Ipoint<dim>& ipoint, RotInvarMethod rotinv)const;

	//! Checks whether computing a descriptor here is allowed
	bool ValidIndex(const IndexType& idx, double scale)const;

	static int SupportRegionWidth(double scale);

	//! Minimum allowed image index
	static int MinIndex(double scale);

	/*! Maximum allowed image index
		\param imageSize number of pixels in this dimension
	*/
	static int MaxIndex(double scale, int imageSize);



private:

	void getOrientationFromGradients_MeanGrad(
		const std::vector<vnl_vector_fixed<double, dim> >& gradients,
		const vnl_vector_fixed<double, dim>& meanGrad,
		typename Ipoint<dim>::RotMatrixType& rot)const;

	void getOrientationFromGradients_PCA(
		const std::vector<vnl_vector_fixed<double, dim> >& gradients,
		const vnl_vector_fixed<double, dim>& meanGrad,
		typename Ipoint<dim>::RotMatrixType& rot)const;


	//! Calculate Haar wavelet response
	IntPixelType haar(int dimension, const IndexType& pos, int size)const;

	void MakeIndexImagesAndGaussian();

	//---------------- Private Variables -----------------//

	//! Integral image where Ipoints have been detected
	IntegralImageType m_img;

	//! Size of integral image (for convenience)
	RegionType m_lpr;

	//! Ipoints vector
	std::vector<Ipoint<dim> > *m_ipts;

	typedef itk::Image<IndexType, dim> IndexImageType;

	//! Precomputed image with bin offset indices
	typename IndexImageType::Pointer m_binImage;

	//! Precomputed container with indices for orientation detection
	std::vector<IndexType> m_orientationIndices;

	//! Precomputed image with intra-bin-indices
	typename IndexImageType::Pointer m_binSampleImage;

	//! Precomputed gaussian kernel
	typedef itk::Image<double, dim> FloatImageType;
	typename FloatImageType::Pointer m_gauss33;
};

}

// ------------------------ Implementation ------------------------- //

#include <stdexcept>

#include <itkImageRegionConstIterator.h>
#include <vnl/algo/vnl_determinant.h>

namespace nsurf
{

//! Round float to nearest integer
inline int fRound(double flt)
{
	return (int) floor(flt+0.5);
}


//! Round float vector to nearest integer
template<typename TVec>
void vfRound(TVec& v)
{
	typedef typename TVec::value_type value_type;
	for (unsigned d=0; d < v.size(); ++d)
		v[d] = value_type(floor(flt+0.5f))
}


template<typename TPixel, int dim>
int Surf<TPixel, dim>::
SupportRegionWidth(double scale)
{
	int BinsPerDim = DescriptorDimensions<Dimension>::NumSubBinsPerDimension;
	int halfSpan = DescriptorDimensions<Dimension>::NumSamplesPerDimAndBin*BinsPerDim / 2;

	int minIndex = fRound(-scale*halfSpan);
	int maxIndex = fRound( scale*(halfSpan - 1));
	return maxIndex - minIndex + 1;
}


template<typename TPixel, int dim>
bool Surf<TPixel, dim>::ValidIndex(const IndexType& idx, double scale)const
{
	const itk::Size<dim>& size = m_lpr.GetSize();
	for (int d=0; d < dim; ++d)
	{
		if (idx[d] < MinIndex(scale))
			return false;
		if (idx[d] > MaxIndex(scale, size[d]))
			return false;
	}
	return true;
}


template<typename TPixel, int dim>
int Surf<TPixel, dim>::MinIndex(double scale)
{
	int BinsPerDim = DescriptorDimensions<Dimension>::NumSubBinsPerDimension;
	int halfSpan = DescriptorDimensions<Dimension>::NumSamplesPerDimAndBin*BinsPerDim / 2;
	int minIndex = fRound(-scale*halfSpan);
	return -minIndex;
}

template<typename TPixel, int dim>
int Surf<TPixel, dim>::MaxIndex(double scale, int imageSize)
{
	int BinsPerDim = DescriptorDimensions<Dimension>::NumSubBinsPerDimension;
	int halfSpan = DescriptorDimensions<Dimension>::NumSamplesPerDimAndBin*BinsPerDim / 2;
	int maxIndex = fRound( scale*(halfSpan - 1));
	return imageSize-1-maxIndex;
}

template<typename TPixel, int dim>
Surf<TPixel, dim>::
Surf(const IntegralImageType& img, std::vector<Ipoint<dim> > *ipts)
	: m_ipts(ipts)

{
	m_img = img;
	m_lpr = m_img.GetIntegralImage()->GetLargestPossibleRegion();
	MakeIndexImagesAndGaussian();
}

template<typename TPixel, int dim>
void Surf<TPixel, dim>::
MakeIndexImagesAndGaussian()
{
	int BinsPerDim = DescriptorDimensions<Dimension>::NumSubBinsPerDimension;

	int halfSpan = DescriptorDimensions<Dimension>::NumSamplesPerDimAndBin*BinsPerDim / 2;
	itk::Size<dim> size;

	typedef itk::ImageRegionIteratorWithIndex<IndexImageType> IndexIterType;
	typedef itk::ImageRegionIteratorWithIndex<FloatImageType> FloatIterType;

	{
		size.Fill(BinsPerDim);
		m_binImage = IndexImageType::New();
		m_binImage->SetRegions(size);
		m_binImage->Allocate();

		for (IndexIterType i(m_binImage, m_binImage->GetLargestPossibleRegion()); !i.IsAtEnd(); ++i)
		{
			const itk::Index<dim>& idx = i.GetIndex();
			itk::Index<dim> pixel;
			for (unsigned d=0; d < Dimension; ++d)
				pixel[d] = idx[d]*DescriptorDimensions<Dimension>::NumSamplesPerDimAndBin - halfSpan;
			i.Set(pixel);
		}
	}
	{
		size.Fill(DescriptorDimensions<Dimension>::NumSamplesPerDimAndBin);
		m_binSampleImage = IndexImageType::New();
		m_binSampleImage->SetRegions(size);
		m_binSampleImage->Allocate();

		for (IndexIterType i(m_binSampleImage, m_binSampleImage->GetLargestPossibleRegion()); !i.IsAtEnd(); ++i)
		{
			const itk::Index<dim>& idx = i.GetIndex();
			i.Set(idx);
		}
	}
	{
		int radius = 6;
		int step = 1;

		int maxssq = radius*radius;
		size.Fill(2*radius+1);
		typename ImageType::Pointer myimage = ImageType::New();
		myimage->SetRegions(size);
		typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IterType;
		for (IterType i(myimage, myimage->GetLargestPossibleRegion()); !i.IsAtEnd(); ++i)
		{
			itk::Index<dim> idx = i.GetIndex();
			int ssq = 0;
			for (int d=0; d < Dimension; ++d)
			{
				idx[d] -= radius;
				ssq += idx[d]*idx[d];
				idx[d] *= step;
			}
			if (ssq < maxssq)
				m_orientationIndices.push_back(idx);
		}
	}
	{
		//double sigma = 3.3;
		double sigma = 10.0;
		double var = sigma*sigma;
		static const double pi = 3.1415926535897932384626433832795;

		double prefac = pow(2*pi, Dimension * 0.5 * -1)*pow(var, Dimension*0.5*-1);
		double isigsq = 1.0/var;
		size.Fill(halfSpan + 1);
		m_gauss33 = FloatImageType::New();
		m_gauss33->SetRegions(size);
		m_gauss33->Allocate();

		for (FloatIterType i(m_gauss33, m_gauss33->GetLargestPossibleRegion()); !i.IsAtEnd(); ++i)
		{
			const itk::Index<dim>& idx = i.GetIndex();
			double exarg = 0;
			for (int d=0; d < Dimension; ++d)
				exarg += idx[d]*idx[d];
			exarg = -0.5*exarg*isigsq;
			double val = prefac * exp(exarg);
			i.Set(double(val));
		}
	}

}


template<typename TPixel, int dim>
typename Surf<TPixel, dim>::IntPixelType
Surf<TPixel, dim>::haar(int dimension, const IndexType& pos, int size)const
{
	itk::Size<dim> nsize;
	nsize.Fill(size);
	int halfsize = size/2;
	itk::Index<dim> index = pos;
	for (int d=0; d < dim; ++d)
		index[d] -= halfsize;
	RegionType reg(index, nsize);
	if (!m_lpr.IsInside(reg))
		return 0;
	reg.SetSize(dimension, halfsize);

	IntPixelType result = - m_img.GetIntegral(reg);
	reg.SetIndex(dimension, pos[dimension]);
	result += m_img.GetIntegral(reg);
	return result;
}


template<typename TPixel, int dim>
void Surf<TPixel, dim>::getDescriptors(RotInvarMethod rotinv)const
{
	int ipts_size = int(m_ipts->size());
	// U-SURF loop just gets descriptors
	for (int i = 0; i < ipts_size; ++i)
	{
		if (rotinv != upright)
			getOrientation((*m_ipts)[i], rotinv);
		getDescriptor((*m_ipts)[i], rotinv);
	}
}


template<typename TPixel, int dim>
void Surf<TPixel, dim>::getOrientationFromGradients_PCA(
	const std::vector<vnl_vector_fixed<double, dim> >& gradients,
	const vnl_vector_fixed<double, dim>& meanGrad,
	typename Ipoint<dim>::RotMatrixType& rot)const
{
	typedef vnl_vector_fixed<double, dim> vec_t;
	// find PCA
	GaussianDist<dim> gaussian;
	gaussian.Learn(gradients);
	vec_t eigenvalues;
	gaussian.GetPrincipalComponents(rot, eigenvalues);
	// normalize columns
	rot.normalize_columns();


	// make colums point into mean gradient direction
	for (int c=0; c < dim-1; ++c)
	{
		double sum=0;
		for (int r=0; r < dim; ++r)
			sum += rot[r][c]*meanGrad[r];
		if (sum < 0)
			rot.scale_column(c, -1);
	}

	// compute average gradient for pos/neg axis
	// make sure that det(A) = 1
	if (vnl_determinant(rot) < 0)
		rot.scale_column(Dimension-1, -1);
}

template<int dim>
void MakeBasis(const vnl_vector_fixed<double, dim>& firstAxis, vnl_matrix_fixed<double, dim, dim>& m)
{
	typedef vnl_vector_fixed<double, dim> vec_t;
	vec_t one(1);
	vec_t n = firstAxis;
	n.normalize();
	vec_t basis[dim];
	basis[0] = n;


	for (int i=1; i < dim; ++i)
	{
		int dimrot = 0;
		double nvlen;
		vec_t nv;
		do
		{
			nv.fill(0);
			nv[dimrot] = 1;


			for (int n=0; n < i; ++n)
			{
				const vec_t& lr = basis[n];
				nv -= dot_product(nv, lr)*lr;
				assert(abs(dot_product(nv, lr)) < 1e-5);
			}
			nvlen = nv.magnitude();
			++dimrot;
		}while(nvlen < 0.01 && dimrot < dim);
		assert(!(nvlen < 0.01));

		nv /= nvlen;
		basis[i] = nv;
#ifndef NDEBUG
		for (int n=0; n < i; ++n)
			assert(abs(dot_product(nv, basis[n])) < 1e-5);
#endif
	}
	for (int i=0; i < dim; ++i)
		m.set_column(i, basis[i]);

	assert(abs(vnl_determinant(m)) > 0.01);
}


template<typename TPixel, int dim>
void Surf<TPixel, dim>::getOrientationFromGradients_MeanGrad(
	const std::vector<vnl_vector_fixed<double, dim> >& gradients,
	const vnl_vector_fixed<double, dim>& meanGrad,
	typename Ipoint<dim>::RotMatrixType& rot)const
{
	typedef vnl_vector_fixed<double, dim> vec_t;
	typedef vnl_matrix_fixed<double, dim, dim> mat_t;

	vec_t n0 = meanGrad;
	double meanGradLen = meanGrad.magnitude();
	if (meanGradLen < 1e-5)
	{	// bad conditioned. Return identity
		rot.fill(0.0);
		rot.fill_diagonal(1.0);
		return;
	}
	n0 /= meanGradLen;

	mat_t basis, toBasis;
	MakeBasis(n0, basis);
	toBasis = basis.transpose();

#ifndef NDEBUG
	vec_t e0(0.0);
	e0[0] = 1;
	assert((toBasis*n0 - e0).magnitude() < 1e-5);
#endif

	typedef vnl_vector_fixed<double, dim-1> vecl_t;

	std::vector<vecl_t> projgrad;
	projgrad.resize(gradients.size());


	for (size_t i=0; i < gradients.size(); ++i)
	{
		const vec_t& v = gradients[i];
		// project
		vec_t pv = v - dot_product(n0,v)*n0;
		assert(abs(dot_product(pv,meanGrad)) < 1e-5);
		// transform
		pv = toBasis*pv;
		assert(abs(pv[0]) < 1e-5);
		vecl_t lv;
		for (int d=0; d < dim-1; ++d)
			lv[d] = pv[d+1];
		projgrad[i] = lv;
	}

	// find PCA
	GaussianDist<dim-1> gaussian;
	gaussian.Learn(projgrad);
	typedef vnl_matrix_fixed<double, dim-1, dim-1> matl_t;
	vecl_t eigenvalues;
	matl_t eigenvectors;
	gaussian.GetPrincipalComponents(eigenvectors, eigenvalues);
	// normalize columns
	eigenvectors.normalize_columns();

	// mean should be zero
#ifndef NDEBUG
	vecl_t sum(0.0);
	for (size_t i=0; i < projgrad.size(); ++i)
		sum += projgrad[i];
	sum /= double(projgrad.size());
	assert(sum.magnitude() < 1e-4);
#endif

	// compute second order moment

	vecl_t m2(0.0);
	for (size_t i=0; i < projgrad.size(); ++i)
		m2 += projgrad[i]*projgrad[i].squared_magnitude();

	// make colums point into second order moment direction
	for (int c=0; c < dim-2; ++c)
	{
		double sum=0;
		for (int r=0; r < dim; ++r)
			sum += eigenvectors[r][c]*m2[r];
		if (sum < 0)
			eigenvectors.scale_column(c, -1);
	}

	// backtransform
	rot.set_column(0, n0);
	for (int i=1; i < dim; ++i)
	{
		vec_t v;
		v[0] = 0;
		vecl_t lv = eigenvectors.get_column(i-1);
		for (int d=0; d < dim-1; ++d)
			v[d+1] = lv[d];
		v = basis*v;
		rot.set_column(i,v);
	}

	// compute average gradient for pos/neg axis
	// make sure that det(A) = 1
	if (vnl_determinant(rot) < 0)
		rot.scale_column(Dimension-1, -1);
}


template<typename TPixel, int dim>
void Surf<TPixel, dim>::getOrientation(Ipoint<dim>& ipoint, RotInvarMethod rotinv)const
{
	typedef vnl_vector_fixed<double, dim> vec_t;
	typedef Ipoint<dim>::RotMatrixType mat_t;
	std::vector<vec_t> gradients;
	vec_t meanGradient(0.0);
	int haarSize = fRound(ipoint.scale*2);
	int minIndex = haarSize/2;
	itk::Index<3> maxIndex;
	for (int d=0; d < Dimension; ++d)
		maxIndex[d] = m_img.GetIntegralImage()->GetLargestPossibleRegion().GetSize(d)-haarSize;

	// extract gradients in circular region
	gradients.reserve(m_orientationIndices.size());
	for (size_t i=0; i < m_orientationIndices.size(); ++i)
	{
		itk::Index<dim> idx = ipoint.index;
		for (int d=0; d < Dimension; ++d)
		{
			idx[d] += ipoint.scale*m_orientationIndices[i][d];
			if (idx[d] < minIndex || idx[d] > maxIndex[d])
				continue;
		}
		vec_t grad;
		for (int d=0; d < Dimension; ++d)
			grad[d] = double(haar(d, idx, haarSize));
		meanGradient += grad;
		gradients.push_back(grad);
	}

	switch (rotinv)
	{
		case SurfBase::MeanGradPCA:
			getOrientationFromGradients_MeanGrad(gradients, meanGradient, ipoint.rotMatrix);
			break;

		case SurfBase::PCA:
			getOrientationFromGradients_PCA(gradients, meanGradient, ipoint.rotMatrix);
			break;

		default:
			assert(0);
	}

	bool plotgradients = false;
	if (plotgradients)
	{

		static bool firstrun = true;


		{
			mat_t rotmt = ipoint.rotMatrix.transpose();
			char name[128];
			{
				sprintf(name, "gradients_%dD.txt", dim);
				std::ofstream s(name, firstrun ? std::ios::out : std::ios::out | std::ios::app);
				for (size_t i=0; i < gradients.size(); ++i)
					s << gradients[i] << "\n";
				s << "\n\n";
			}
			{
				sprintf(name, "rot_%dD.txt", dim);
				std::ofstream s(name, firstrun ? std::ios::out : std::ios::out | std::ios::app);
				// output transposed matrix
				//s << rotmt << "\n\n";

				for (unsigned r=0; r < dim; ++r)
				{
					for (unsigned c=0; c < dim; ++c)
						s << 0 << " ";
					s << "\n";
					for (unsigned c=0; c < dim; ++c)
					{
						double fac = (r==0 ? 200 : 100);
						s << rotmt(r,c)*100 << " ";
					}
					s << "\n";
				}
				s << "\n\n";
			}
			if (dim == 3)
			{
				sprintf(name, "meangradplane_%dD.txt", dim);
				std::ofstream s(name, firstrun ? std::ios::out : std::ios::out | std::ios::app);
				vec_t a1 = rotmt.get_row(1);
				vec_t a2 = rotmt.get_row(2);
				vec_t c[4];
				c[0] = -a1-a2;
				c[1] = -a1+a2;
				c[2] = a1+a2;
				c[3] = a1-a2;

				for (int i=0; i < 5; ++i)
					s << c[i%4] << "\n";
				s << "\n\n";
			}
			{
				sprintf(name, "meangrad_%dD.txt", dim);
				std::ofstream s(name, firstrun ? std::ios::out : std::ios::out | std::ios::app);
				vec_t origin(0.0);
				// draw mean grad:
				s << origin << "\n" << rotmt.get_row(0) << "\n";
				s << "\n\n";
			}
			firstrun = false;
		}
	}



#ifndef NDEBUG
	double deter = vnl_determinant(ipoint.rotMatrix);
	assert(abs(deter-1) < 1e-5);
#endif
}


//! Get the upright descriptor vector of the provided Ipoint
template<typename TPixel, int dim>
void Surf<TPixel, dim>::getDescriptor(Ipoint<dim>& ipoint, RotInvarMethod rotinv)const
{
	assert(rotinv == upright || abs(vnl_determinant(ipoint.rotMatrix) - 1.0) < 1e-4);
	const IndexType& idx = ipoint.index; // was: x, y
	int count=0;
	typedef vnl_vector_fixed<double, dim> vec_t;
	vec_t delta, mdelta, r; // was: dx, dy, mdx, mdy, rx, ry
	double scale;
	float *desc;
	double len = 0;

	Ipoint<dim> *ipt = &ipoint;
	scale = ipt->scale;

	int haarSize = 2*fRound(scale);
	IndexType pos = ipt->index; // was. x, y
	Ipoint<dim>::RotMatrixType invRotMat = ipt->rotMatrix;
	invRotMat.inplace_transpose();

	// check bounds
/*
	if (!ValidIndex(pos, scale))
	{
		ipt->valid = false;
		return;
	}
*/
	desc = ipt->descriptor;

#define WEIGHTPERSAMPLE

	typedef itk::ImageRegionConstIterator<IndexImageType> IndexIterType;
	// iterate over (typically) 4^d bins
	for (IndexIterType i(m_binImage, m_binImage->GetLargestPossibleRegion()); !i.IsAtEnd(); ++i)
	{
		const IndexType& binIdx = i.Get();
		delta.fill(0);
		mdelta.fill(0);

#ifndef WEIGHTPERSAMPLE
		IndexType binCenterIdx = binIdx;// -8 .. 7
		for (unsigned d=0; d < Dimension; ++d)
			binCenterIdx[d] += DescriptorDimensions<Dimension>::NumSamplesPerDimAndBin/2;
		IndexType binCenterGaussIdx = binCenterIdx;
		for (unsigned d=0; d < Dimension; ++d)
			binCenterGaussIdx[d] = binCenterGaussIdx[d] >= 0 ? binCenterGaussIdx[d]+1 : -binCenterGaussIdx[d];
#endif

		// iterate over NumSamplesPerDimAndBin^d sample points within this bin
		for (IndexIterType j(m_binSampleImage, m_binSampleImage->GetLargestPossibleRegion());
			!j.IsAtEnd(); ++j)
		{
			const IndexType& sampleIdx = j.Get(); // 0...4
			IndexType absIdx;
			for (unsigned d=0; d < Dimension; ++d)
				absIdx[d] = sampleIdx[d] + binIdx[d]; // -10 ... 9
#ifdef WEIGHTPERSAMPLE
			IndexType gaussIdx = absIdx;
			for (unsigned d=0; d < Dimension; ++d)
				//gaussIdx[d] = abs(gaussIdx[0]);
				gaussIdx[d] = gaussIdx[d] >= 0 ? gaussIdx[d]+1 : -gaussIdx[d];

			// compute haar wavelet responses for each dimension
			// weight with gaussian kernel centered at pos
			double gauss = m_gauss33->GetPixel(gaussIdx);
#endif
			if (rotinv == upright)
			{
				IndexType haarIndex;
				for (unsigned d=0; d < Dimension; ++d)
					haarIndex[d] = idx[d] + fRound(scale*absIdx[d]);
				for (unsigned d=0; d < Dimension; ++d)
				{
					r[d] = haar(d, haarIndex, haarSize);
#ifdef WEIGHTPERSAMPLE
					r[d] *= gauss;
#endif
				}
			}
			else
			{
				IndexType haarIndex;
				typedef vnl_vector_fixed<double, dim> vecd_t;
				vecd_t rhaarIdx;
				for (unsigned d=0; d < Dimension; ++d)
					rhaarIdx[d] = scale*absIdx[d];
				rhaarIdx = ipt->rotMatrix*rhaarIdx;
				for (unsigned d=0; d < Dimension; ++d)
					haarIndex[d] = int(rhaarIdx[d]+0.5) + idx[d];
				for (unsigned d=0; d < Dimension; ++d)
				{
					r[d] = haar(d, haarIndex, haarSize);
#ifdef WEIGHTPERSAMPLE
					r[d] *= gauss;
#endif
				}
				r = invRotMat*r;
			}
			// compute statistics
			delta += r;
			for (unsigned d=0; d < Dimension; ++d)
				mdelta[d] += fabs(r[d]);
		}
#ifndef WEIGHTPERSAMPLE
		double gauss = m_gauss33->GetPixel(binCenterGaussIdx);
		delta *= gauss;
		mdelta *= gauss;
#endif

		for (unsigned d=0; d < Dimension; ++d)
		{
			assert(delta[d] < 0 || delta[d] >= 0);
			desc[count++] = float(delta[d]);
		}
		for (unsigned d=0; d < Dimension; ++d)
		{
			assert(mdelta[d] < 0 || mdelta[d] >= 0);
			desc[count++] = float(mdelta[d]);
		}

		len += delta.squared_magnitude() + mdelta.squared_magnitude();
	}
	assert(count == DescriptorDimensions<Dimension>::DescriptorSize);
	len = sqrt(len);
	assert(abs(len) >= 1e-6);
	double ilen = double(1.0/len);
	for(int i = 0; i < DescriptorDimensions<Dimension>::DescriptorSize; i++)
		desc[i] *= ilen;
	ipt->valid = true;
}
#ifdef WEIGHTPERSAMPLE
#undef WEIGHTPERSAMPLE
#endif
}

#endif // Surf_h__
