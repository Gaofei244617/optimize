#include "Individual.h"
#include <string>

namespace opt
{
	// ���캯����n:����������str:��������Ⱥ
	Individual::Individual() : nVars(0), vars(nullptr), fitness(0)
	{
	}

	// ���캯����n:����������str:��������Ⱥ
	Individual::Individual(const std::size_t n) : nVars(n), fitness(0)
	{
		vars = new double[n]();
	}

	// ���캯��
	Individual::Individual(const std::initializer_list<double>& list) : nVars(list.size()), fitness(0)
	{
		vars = new double[nVars]();
		for (std::size_t i = 0; i < nVars; i++)
		{
			vars[i] = *(list.begin() + i);
		}
	}

	// ���ƹ���
	Individual::Individual(const Individual& other) : nVars(other.nVars), fitness(other.fitness)
	{
		vars = new double[other.nVars];
		for (std::size_t i = 0; i < nVars; i++)
		{
			vars[i] = (other.vars[i]);
		}
	}

	// �ƶ�����
	Individual::Individual(Individual&& other)noexcept
		:nVars(other.nVars),
		vars(other.vars),
		fitness(other.fitness)
	{
		other.vars = nullptr;
	}

	// ��ֵ����
	Individual& Individual::operator=(const Individual& other)noexcept
	{
		// �����Ը�ֵ
		if (this != &other)
		{
			// ��������������ȣ��׳��쳣
			if (nVars != other.nVars)
			{
				delete[] vars;
				vars = new double[other.nVars]();
				nVars = other.nVars;
			}

			// ��������������ֵ
			for (std::size_t i = 0; i < nVars; i++)
			{
				vars[i] = (other.vars)[i];
			}
			fitness = other.fitness;
		}
		return *this;
	}

	// �ƶ���ֵ
	Individual& Individual::operator=(Individual&& other)noexcept
	{
		// �����Ը�ֵ
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

	// ��������
	Individual::~Individual()
	{
		delete[] vars;
	}
}