#include "Range_Random.h"
#include <random>
#include <ctime>

namespace opt
{
	// 随机数引擎，用于生成一个随机的unsigned整数
	// 该引擎对象重载了()运算符
	static std::default_random_engine rand_eng((unsigned int)time(NULL));

	// 生成区间[min, max)内的 double 类型的随机数
	double random_real(const double min, const double max)
	{
		return std::uniform_real_distribution<double>(min, max)(rand_eng);
	}

	// 生成区间[min, max]内的 int 类型的随机数
	int random_int(const int min, const int max)
	{
		return std::uniform_int_distribution<int>(min, max)(rand_eng);
	}
}