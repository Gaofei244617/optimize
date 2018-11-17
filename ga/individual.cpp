#include "Individual.h"
#include <string>

namespace opt
{
	// ���캯����n:����������str:��������Ⱥ
	Individual::Individual(const std::size_t n) : nVars(n), fitness(0)
	{
		vars = new double[n]();
	}

	// ���캯��
	Individual::Individual(const std::initializer_list<double>& list)
		: nVars(list.size()),
		fitness(0)
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
		: nVars(other.nVars),
		vars(other.vars),
		fitness(other.fitness)
	{
		other.vars = nullptr;
	}

	// ��ֵ����
	Individual& Individual::operator=(const Individual& other)
	{
		// �����Ը�ֵ
		if (this != &other)
		{
			// ��������������ȣ��׳��쳣
			if (nVars != other.nVars)
			{
				throw std::string("Two individuals do not match.");
			}

			// ��������������ֵ
			for (std::size_t i = 0; i < nVars; i++)
			{
				vars[i] = (other.vars[i]);
			}
			fitness = other.fitness;
		}
		return *this;
	}

	// �ƶ���ֵ
	Individual& Individual::operator=(Individual&& other)
	{
		// �����Ը�ֵ
		if (this != &other)
		{
			// ��������������ȣ��׳��쳣
			if (nVars != other.nVars)
			{
				throw std::string("Two individuals do not match.");
			}
			vars = other.vars;
			other.vars = nullptr;
			fitness = other.fitness;
		}
		return *this;
	}

	// ��������
	Individual::~Individual()
	{
		delete[] vars;
		vars = nullptr;
	}
}