#include "Individual.h"
#include <string>

namespace opt
{
	// ���캯����n:����������str:��������Ⱥ
	Individual::Individual(int n) : nVars(n), fitness(0)
	{
		vars = new double[n]();
	}

	// ���ƹ���
	Individual::Individual(const Individual& other) : nVars(other.nVars), fitness(other.fitness)
	{
		vars = new double[other.nVars];
		for (int i = 0; i < nVars; i++)
		{
			vars[i] = (other.vars[i]);
		}
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
			for (int i = 0; i < nVars; i++)
			{
				vars[i] = (other.vars[i]);
			}
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