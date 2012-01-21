/********************************************************************
	Year:      2009
	Author:    Johannes Feulner
*********************************************************************/

#ifndef IntegralImage_h__
#define IntegralImage_h__


/*!	Class for integral images of arbitrary dimension
*/

#include <vnl/vnl_vector_fixed.h>
#include <itkImage.h>
#include "IntegralType.h"
#include <vector>
#include <algorithm>
#include <assert.h>

template<int dim>
struct PermutTypes
{
	typedef vnl_vector_fixed<int, dim> TupleType;
	typedef std::vector<TupleType> TupleVectorType;
	typedef std::vector<TupleVectorType> TupleVectorVectorType;
};

template<int dim>
struct IndexPermutations : public PermutTypes<dim>
{
public:
	IndexPermutations();

	//! GetPermutations()[0][0] is (1 1 ... 1)
	//! GetPermutations()[dim][0] is (0 0 ... 0)
	const TupleVectorVectorType& GetPermutations() { return permutations; }

private:
	void generate();
	TupleVectorVectorType permutations;
};


//template<int dim> typename IndexPermutations<dim>::TupleVectorVectorType IndexPermutations<dim>::permutations;

template<int dim>
IndexPermutations<dim>::IndexPermutations()
{
	//if (permutations.empty())
	generate();
}

#ifndef NDEBUG

inline unsigned choose_k_out_of_n(unsigned n, unsigned k) {
	if (k > n)
		return 0;

	if (k > n/2)
		k = n-k; // Take advantage of symmetry

	long double accum = 1;
	for (unsigned i = 1; i <= k; i++)
		accum = accum * (n-k+i) / i;

	return unsigned(accum + 0.5); // avoid rounding error
}

#endif


template<int dim>
void IndexPermutations<dim>::generate()
{
	permutations.resize(dim+1);
	for (size_t d=0; d < permutations.size(); ++d)
	{
		TupleType tuple(1);
		for (size_t n=0; n < d; ++n)
			tuple[n] = 0;
		do
		{
			permutations[d].push_back(tuple);
		} while (std::next_permutation(tuple.begin(), tuple.end()));
		assert(permutations[d].size() == choose_k_out_of_n(dim, d));
	}
}

template<typename TSourcePixel, int dim>
class IntegralImage : public PermutTypes<dim>
{
public:
	typedef TSourcePixel SourcePixelType;
	typedef typename IntegralType<SourcePixelType>::PixelType IntPixelType;
	static const unsigned int Dimension = dim;

	typedef itk::Image<SourcePixelType, dim> SourceType;
	typedef itk::Image<IntPixelType, dim> IntType;
	typedef typename IntType::Pointer IntPointer;
	typedef typename IntType::ConstPointer IntConstPointer;
	typedef itk::ImageRegion<dim> RegionType;
	typedef itk::Index<dim> IndexType;

	void Compute(const SourceType* source);
	void clear();

	IntPointer      GetIntegralImage()      { return m_integralImage; }
	IntConstPointer GetIntegralImage()const { return m_integralImage.GetPointer(); }

	//! No boundary checks
	IntPixelType GetIntegral(const RegionType& reg)const;

	//! Crops reg to fit image if necessary
	IntPixelType GetIntegralSafe(const RegionType& reg)const;

private:
	IntPointer m_integralImage;
	static IndexPermutations<dim> permut;
};

template<typename TSourcePixel, int dim> IndexPermutations<dim> IntegralImage<TSourcePixel,dim>::permut;

#include <itkCastImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>


template<typename TSourcePixel, int dim>
void IntegralImage<TSourcePixel, dim>::
clear()
{
	m_integralImage = 0;
}


template<typename TSourcePixel, int dim>
void IntegralImage<TSourcePixel, dim>::
Compute(const SourceType* source)
{
	typedef itk::CastImageFilter<SourceType, IntType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(source);
	caster->Update();
	m_integralImage = caster->GetOutput();
	RegionType lpr = m_integralImage->GetLargestPossibleRegion();

	typedef itk::ImageRegionIteratorWithIndex<IntType> IterType;

	for (unsigned d=0; d < Dimension; ++d)
	{
		for (IterType i(m_integralImage, lpr); !i.IsAtEnd(); ++i)
		{
			const IndexType& idx = i.GetIndex();
			if (idx[d] == 0)
				continue;
			IndexType leftidx = idx;
			--leftidx[d];
			i.Set(i.Get() + m_integralImage->GetPixel(leftidx));
		}
	}
}

template<typename TSourcePixel, int dim>
typename IntegralImage<TSourcePixel, dim>::IntPixelType
IntegralImage<TSourcePixel, dim>::
GetIntegral(const RegionType& reg)const
{
	assert(m_integralImage->GetLargestPossibleRegion().IsInside(reg));

	// stores for each of the 2^d sub-regions starting at the origin the last index
	IndexType lowhigh[2];

	for (unsigned d=0; d < Dimension; ++d)
	{
		lowhigh[1][d] = reg.GetIndex(d)+reg.GetSize(d)-1;
		lowhigh[0][d]  = reg.GetIndex(d)-1;
	}
	IntPixelType result = 0;
	IntPixelType sign = 1;

	for (unsigned b=0; b <= Dimension; ++b)
	{
		const TupleVectorType& tuples = permut.GetPermutations()[b];

		for (size_t t=0; t < tuples.size(); ++t)
		{
			IndexType idx;
			const TupleType& tuple = tuples[t];
			bool validbox = true;
			assert(int(validbox) == 1);
			for (unsigned d=0; d < Dimension; ++d)
			{
				idx[d] = lowhigh[tuple[d]][d];
				validbox &= (idx[d] >= 0);
				assert(int(idx[d] >= 0) == 1 || (idx[d] >= 0) == 0);
			}
			if (!validbox) // happens at the border of the volume
				continue;
			IntPixelType boxint = m_integralImage->GetPixel(idx);
			result += boxint*sign;
		}
		sign = -sign;
	}
	return result;
}

template<typename TSourcePixel, int dim>
typename IntegralImage<TSourcePixel, dim>::IntPixelType
IntegralImage<TSourcePixel, dim>::
GetIntegralSafe(const RegionType& reg)const
{
	RegionType safereg = reg;
	safereg.Crop(m_integralImage->GetLargestPossibleRegion());
	return GetIntegral(safereg);
}

#endif // IntegralImage_h__
