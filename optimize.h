#ifndef _OPTIMIZE_ALGORITHM_
#define _OPTIMIZE_ALGORITHM_

#include "ga\ga.h"

namespace opt
{
	// Helper Function: Create a GA group
	template<class R, class... Args>
	GAGroup<R(Args...)> createGAGroup(R(*func)(Args...), const int N = 1000)
	{
		return GAGroup<R(Args...)>(func, N);
	}
}

#endif
