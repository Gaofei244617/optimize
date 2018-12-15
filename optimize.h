#ifndef _OPTIMIZE_ALGORITHM_
#define _OPTIMIZE_ALGORITHM_

#include "ga\ga.h"
#include "pso\pso.h"

namespace opt
{
	/******************* Helper Functions *******************/
	// Create a GA group
	template<class R, class... Args>
	GAGroup<R(Args...)> createGAGroup(R(*func)(Args...), const int N = 1000)
	{
		return GAGroup<R(Args...)>(func, N);
	}

	// Create a PSO group
	template<class R, class... Args>
	PSO<R(Args...)> createPSO(R(*func)(Args...), const int N = 1000)
	{
		return PSO<R(Args...)>(func, N);
	}
}

#endif
