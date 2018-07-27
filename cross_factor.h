#ifndef _CROSS_FACTOR_H_
#define _CROSS_FACTOR_H_

#include <utility>
#include <cmath>
#include "Range_Random.h"

namespace opt
{
	// SBX交叉算法(模拟二进制单点交叉)
	// 子代满足 (1) p1 + p2 = c1 + c2, (2) beta = |(c2 - c1) / (p2 - p1)|
	// p1、p2为父代, c1、c2为子代, 一对父代交叉产生两个子代

	std::pair<double, double> cross_SBX(const double m, const double f)
	{
		double beta = 0;                              // 均匀分布因子
		double nc = 1.0;                              // 交叉分布指数(大于0)，推荐为1; nc越大：子代个体离父代越远

		double rand_u = random_real(0, 1);            // 随机数缓存 
		if (rand_u <= 0.5)
		{
			beta = std::pow(2.0 * rand_u, 1.0 / (nc + 1.0));
		}
		else
		{
			beta = std::pow(1.0 / (2.0 * (1.0 - rand_u)), 1.0 / (nc + 1));
		}

		// 基因交叉产生子代基因
		double mm = 0.5 * ((1 + beta)*m + (1 - beta)*f);
		double ff = 0.5 * ((1 - beta)*m + (1 + beta)*f);

		return std::make_pair(mm, ff);
	}
}

#endif // !_CROSS_FACTOR_H_

