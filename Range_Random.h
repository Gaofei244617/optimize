/*
指定区间的随机数生成函数
*/
#ifndef _RANGE_RANDOM_
#define _RANGE_RANDOM_

namespace opt
{
	// 生成区间[min, max)内的 double 类型的随机数
	double random_real(const double min, const double max);

	// 生成区间[min, max]内的 int 类型的随机数
	int random_int(const int min, const int max);
}

#endif