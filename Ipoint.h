/********************************************************************
	Year:      2009
	Author:    Johannes Feulner
*********************************************************************/

#ifndef Ipoint_h__
#define Ipoint_h__

namespace nsurf
{

template<int dim>
struct DescriptorDimensions
{};

template<>
struct DescriptorDimensions<2>
{
	static const int NumSubBinsPerDimension = 4;
	static const int DescriptorSize = 64;
	static const int NumSamplesPerDimAndBin = 5;
};


#define NUSUBIPEDI 3

#if NUSUBIPEDI == 2

template<>
struct DescriptorDimensions<3>
{
	static const int NumSubBinsPerDimension = 2;
	static const int DescriptorSize = 48;
	static const int NumSamplesPerDimAndBin = 5;
};

#elif NUSUBIPEDI == 3

template<>
struct DescriptorDimensions<3>
{
	static const int NumSubBinsPerDimension = 3;
	static const int DescriptorSize = 162;
	static const int NumSamplesPerDimAndBin = 4;
};

#elif NUSUBIPEDI == 4

template<>
struct DescriptorDimensions<3>
{
	static const int NumSubBinsPerDimension = 4;
	static const int DescriptorSize = 384;
	static const int NumSamplesPerDimAndBin = 4;
};

#endif


template<int dim>
class Ipoint
{

public:
	Ipoint() : valid(false){}

	typedef itk::Index<dim> IndexType;

	//! Coordinates of the detected interest point
	IndexType index;

	//! Detected scale
	double scale;

	typedef vnl_matrix_fixed<double, dim, dim> RotMatrixType;
	// contains rotated unit vectors
	RotMatrixType rotMatrix;

	//! Vector of descriptor components
	float descriptor[DescriptorDimensions<dim>::DescriptorSize];

	bool valid;
};

}

#endif // Ipoint_h__
