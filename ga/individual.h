#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_

#include <initializer_list>

namespace opt
{
	// 个体类
	struct Individual
	{
		const int nVars;                                           // 变量(基因)个数
		double* vars;                                              // 个体各个变量(基因)
		double fitness;                                            // 个体对环境的适应度

		Individual(const int n);                                   // 构造函数
		Individual(const std::initializer_list<double>& list);          // 构造函数
		Individual(const Individual& other);                       // 拷贝构造
		Individual(Individual&& other);                            // 拷贝构造
		Individual& operator=(const Individual& other);            // 赋值函数
		Individual& operator=(Individual&& other);                 // 移动赋值

		~Individual();                                             // 析构
	};
}

#endif