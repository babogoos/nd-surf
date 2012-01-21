/********************************************************************
	Year:      2009
	Author:    Johannes Feulner
*********************************************************************/


#ifndef IntegralType_h__
#define IntegralType_h__

/*!	Use IntegralType<TSourcePixel>::PixelType to get a reasonable
	pixel type of an integral image for a source image with pixel type
	TSourcePixel.
*/

template<typename T>
struct IntegralType
{};

template<>
struct IntegralType<short>
{ typedef int PixelType; };

template<>
struct IntegralType<unsigned short>
{ typedef unsigned int PixelType; };

template<>
struct IntegralType<char>
{ typedef int PixelType; };

template<>
struct IntegralType<unsigned char>
{ typedef unsigned int PixelType; };

template<>
struct IntegralType<int>
{ typedef int PixelType; };

template<>
struct IntegralType<unsigned int>
{ typedef unsigned int PixelType; };

template<>
struct IntegralType<float>
{ typedef double PixelType; };

template<>
struct IntegralType<double>
{ typedef double PixelType; };

#endif //IntegralType_h__
