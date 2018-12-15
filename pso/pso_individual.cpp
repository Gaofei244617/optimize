#include "pso_individual.h"
#include <string>

namespace opt
{
	// 构造函数，n:变量个数，str:隶属的种群
	PSO_Individual::PSO_Individual() 
		: nVars(0), 
		xs(nullptr), 
		vs(nullptr),
		best_xs(nullptr),
		fitness(0)
	{
	}

	// 构造函数，n:变量个数，str:隶属的种群
	PSO_Individual::PSO_Individual(const std::size_t n) 
		: nVars(n), 
		fitness(0)
	{
		xs = new double[n]();
		vs = new double[n]();
		best_xs = new double[n]();
	}

	// 构造函数
	PSO_Individual::PSO_Individual(const std::initializer_list<double>& list) 
		: nVars(list.size()), 
		fitness(0)
	{
		xs = new double[nVars]();
		vs = new double[nVars]();
		best_xs = new double[nVars]();

		for (std::size_t i = 0; i < nVars; i++)
		{
			xs[i] = *(list.begin() + i);
			best_xs[i] = xs[i];
		}
	}

	// 复制构造
	PSO_Individual::PSO_Individual(const PSO_Individual& other) : nVars(other.nVars), fitness(other.fitness)
	{
		xs = new double[other.nVars];
		vs = new double[other.nVars];
		best_xs = new double[other.nVars];

		for (std::size_t i = 0; i < nVars; i++)
		{
			xs[i] = (other.xs[i]);
			vs[i] = (other.vs[i]);
			best_xs[i] = (other.best_xs[i]);
		}
	}

	// 移动构造
	PSO_Individual::PSO_Individual(PSO_Individual&& other)noexcept
		:nVars(other.nVars),
		xs(other.xs),
		vs(other.vs),
		best_xs(other.best_xs),
		fitness(other.fitness)
	{
		other.xs = nullptr;
		other.vs = nullptr;
		other.best_xs = nullptr;
	}

	// 赋值函数
	PSO_Individual& PSO_Individual::operator=(const PSO_Individual& other)noexcept
	{
		// 避免自赋值
		if (this != &other)
		{
			// 若变量个数不相等，抛出异常
			if (nVars != other.nVars)
			{
				delete[] xs;
				delete[] vs;
				delete[] best_xs;

				xs = new double[other.nVars]();
				vs = new double[other.nVars]();
				best_xs = new double[other.nVars]();

				nVars = other.nVars;
			}

			// 拷贝各个变量的值
			for (std::size_t i = 0; i < nVars; i++)
			{
				xs[i] = (other.xs)[i];
				vs[i] = (other.vs)[i];
				best_xs[i] = (other.best_xs)[i];
			}
			fitness = other.fitness;
		}
		return *this;
	}

	// 移动赋值
	PSO_Individual& PSO_Individual::operator=(PSO_Individual&& other)noexcept
	{
		// 避免自赋值
		if (this != &other)
		{
			nVars = other.nVars;
			xs = other.xs;
			vs = other.vs;
			best_xs = other.best_xs;
			other.xs = nullptr;
			other.vs = nullptr;
			other.best_xs = nullptr;
			other.nVars = 0;
			fitness = other.fitness;
			other.fitness = 0;
		}
		return *this;
	}

	// 析构函数
	PSO_Individual::~PSO_Individual()
	{
		delete[] xs;
		delete[] vs;
		delete[] best_xs;
	}
}