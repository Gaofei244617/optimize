#ifndef _PSO_INDIVIDUAL_H_
#define _PSO_INDIVIDUAL_H_

#include <initializer_list>

namespace opt
{
	// 粒子个体
	struct PSO_Individual
	{
		std::size_t nVars;                                                 // 变量个数
		double* xs;                                                        // 粒子各个变量
		double* vs;                                                        // 粒子速度向量
		double* best_xs;                                                   // 粒子所经过的最好位置
		double fitness;                                                    // 个体对环境的适应度

		PSO_Individual();                                                  // 构造函数
		PSO_Individual(const std::size_t n);                               // 构造函数
		PSO_Individual(const std::initializer_list<double>& list);         // 构造函数
		PSO_Individual(const PSO_Individual& other);                       // 拷贝构造
		PSO_Individual(PSO_Individual&& other)noexcept;                    // 移动构造
		PSO_Individual& operator=(const PSO_Individual& other)noexcept;    // 赋值函数
		PSO_Individual& operator=(PSO_Individual&& other)noexcept;         // 移动赋值

		~PSO_Individual();                                                 // 析构
	};
}

#endif