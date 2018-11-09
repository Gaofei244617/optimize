#include "Individual.h"
#include <string>

namespace opt
{
	// 构造函数，n:变量个数，str:隶属的种群
	Individual::Individual(const int n) : nVars(n), fitness(0)
	{
		vars = new double[n]();
	}

	// 复制构造
	Individual::Individual(const Individual& other) : nVars(other.nVars), fitness(other.fitness)
	{
		vars = new double[other.nVars];
		for (int i = 0; i < nVars; i++)
		{
			vars[i] = (other.vars[i]);
		}
	}

	// 赋值构造
	Individual& Individual::operator=(const Individual& other)
	{
		// 避免自赋值
		if (this != &other)
		{
			// 若变量个数不相等，抛出异常
			if (nVars != other.nVars)
			{
				throw std::string("Two individuals do not match.");
			}

			// 拷贝各个变量的值
			for (int i = 0; i < nVars; i++)
			{
				vars[i] = (other.vars[i]);
			}
			fitness = other.fitness;
		}
		return *this;
	}

	// 析构函数
	Individual::~Individual()
	{
		delete[] vars;
		vars = nullptr;
	}
}