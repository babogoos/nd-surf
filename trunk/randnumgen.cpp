/********************************************************************
	Year:      2009
	Author:    Johannes Feulner
*********************************************************************/

#include "randnumgen.h"


#ifdef USE_STOCC
 #include "lib/stocc/stocc.h"
#else
 #include <cstdlib>
#endif

#ifndef _USE_MATH_DEFINES
 #define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <iostream>

using namespace std;


#ifdef USE_STOCC
StochasticLib1 g_randomGen(0);

void randSeed(int seed)
{
	g_randomGen.RandomInit(seed);
}

float randf()
{
	return float(g_randomGen.Random());
}

double randd()
{
	return g_randomGen.Random();
}

double stdNormal()
{
	return g_randomGen.Normal(0,1);
}

#else

void randSeed(int seed)
{
	srand(seed);
}

float randf()
{
	return float(rand())/float(RAND_MAX);
}

double randd()
{
	return double(rand())/double(RAND_MAX);
}

double stdNormal()
{
         double x1, x2, w;

         do {
                 x1 = 2.0 * randd() - 1.0;
                 x2 = 2.0 * randd() - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );
		 return x1 * w;
}

#endif //USE_STOCC
