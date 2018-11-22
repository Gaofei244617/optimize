#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

#include <initializer_list>

namespace opt
{
	// 个体类
	struct Individual
	{
		std::size_t nVars;                                   // 变量(基因)个数
		double* vars;                                              // 个体各个变量(基因)
		double fitness;                                            // 个体对环境的适应度

		Individual();                                              // 构造函数
		Individual(const std::size_t n);                           // 构造函数
		Individual(const std::initializer_list<double>& list);     // 构造函数
		Individual(const Individual& other);                       // 拷贝构造
		Individual(Individual&& other)noexcept;                    // 移动构造
		Individual& operator=(const Individual& other)noexcept;    // 赋值函数
		Individual& operator=(Individual&& other)noexcept;         // 移动赋值

		~Individual();                                             // 析构
	};
}

#endif