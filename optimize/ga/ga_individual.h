#ifndef _GA_INDIVIDUAL_H_
#define _GA_INDIVIDUAL_H_

#include <initializer_list>

namespace opt
{
	// 个体类
	struct GA_Individual
	{
		std::size_t nVars;                                               // 变量(基因)个数
		double* vars;                                                    // 个体各个变量(基因)
		double fitness;                                                  // 个体对环境的适应度

		GA_Individual();                                                 // 构造函数
		GA_Individual(const std::size_t n);                              // 构造函数
		GA_Individual(const std::initializer_list<double>& list);        // 构造函数
		GA_Individual(const GA_Individual& other);                       // 拷贝构造
		GA_Individual(GA_Individual&& other)noexcept;                    // 移动构造
		GA_Individual& operator=(const GA_Individual& other)noexcept;    // 赋值函数
		GA_Individual& operator=(GA_Individual&& other)noexcept;         // 移动赋值

		~GA_Individual();                                                // 析构
	};
}

#endif