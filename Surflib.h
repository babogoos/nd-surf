/********************************************************************
	Year:      2009
	Author:    Johannes Feulner
*********************************************************************/

#ifndef Surflib_h__
#define Surflib_h__

#include "Surf.h"

namespace nsurf
{
	//! Library function describes interest points in vector
	template<typename TPixel, int dim>
	void surfDes(const itk::Image<TPixel, dim>* img,
		std::vector<Ipoint<dim> > &ipts, /* reference to vector of Ipoints */
		SurfBase::RotInvarMethod rotinvar) /* run in rotation invariant mode? */
	{
		// Create integral image representation of the image
		IntegralImage<TPixel, dim> int_img;
		int_img.Compute(img);
		Surf<TPixel, dim> surfgen(int_img, &ipts);

		// Extract the descriptors for the ipts
		surfgen.getDescriptors(rotinvar);
	}
}

#endif // Surflib_h__
