#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_

#include <string>
#include <memory>

namespace opt
{
	// 个体类
	struct Individual
	{
		//std::string groupName;                           // 隶属的种群
		const int nVars;                                   // 变量(基因)个数
		double* vars;                                      // 个体各个变量(基因)
		double fitness;                                    // 个体对环境的适应度

		Individual(int n);                                 // 构造函数
		Individual(const Individual& other);               // 拷贝构造
		Individual& operator=(const Individual& other);    // 赋值构造
		~Individual();                                     // 析构
	};
}

#endif