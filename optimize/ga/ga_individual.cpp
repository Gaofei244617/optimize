#include "ga_individual.h"
#include <string>

namespace opt
{
	// 构造函数，n:变量个数，str:隶属的种群
	GA_Individual::GA_Individual() : nVars(0), vars(nullptr), fitness(0)
	{
	}

	// 构造函数，n:变量个数，str:隶属的种群
	GA_Individual::GA_Individual(const std::size_t n) : nVars(n), fitness(0)
	{
		vars = new double[n]();
	}

	// 构造函数
	GA_Individual::GA_Individual(const std::initializer_list<double>& list) : nVars(list.size()), fitness(0)
	{
		vars = new double[nVars]();
		for (std::size_t i = 0; i < nVars; i++)
		{
			vars[i] = *(list.begin() + i);
		}
	}

	// 复制构造
	GA_Individual::GA_Individual(const GA_Individual& other) : nVars(other.nVars), fitness(other.fitness)
	{
		vars = new double[other.nVars];
		for (std::size_t i = 0; i < nVars; i++)
		{
			vars[i] = (other.vars[i]);
		}
	}

	// 移动构造
	GA_Individual::GA_Individual(GA_Individual&& other)noexcept
		:nVars(other.nVars),
		vars(other.vars),
		fitness(other.fitness)
	{
		other.vars = nullptr;
	}

	// 赋值函数
	GA_Individual& GA_Individual::operator=(const GA_Individual& other)noexcept
	{
		// 避免自赋值
		if (this != &other)
		{
			// 若变量个数不相等，抛出异常
			if (nVars != other.nVars)
			{
				delete[] vars;
				vars = new double[other.nVars]();
				nVars = other.nVars;
			}

			// 拷贝各个变量的值
			for (std::size_t i = 0; i < nVars; i++)
			{
				vars[i] = (other.vars)[i];
			}
			fitness = other.fitness;
		}
		return *this;
	}

	// 移动赋值
	GA_Individual& GA_Individual::operator=(GA_Individual&& other)noexcept
	{
		// 避免自赋值
		if (this != &other)
		{
			nVars = other.nVars;
			vars = other.vars;
			other.vars = nullptr;
			other.nVars = 0;
			fitness = other.fitness;
			other.fitness = 0;
		}
		return *this;
	}

	// 析构函数
	GA_Individual::~GA_Individual()
	{
		delete[] vars;
	}
}