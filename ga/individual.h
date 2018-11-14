#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_

#include <initializer_list>

namespace opt
{
	// ������
	struct Individual
	{
		const int nVars;                                           // ����(����)����
		double* vars;                                              // �����������(����)
		double fitness;                                            // ����Ի�������Ӧ��

		Individual(const int n);                                   // ���캯��
		Individual(const std::initializer_list<double>& list);          // ���캯��
		Individual(const Individual& other);                       // ��������
		Individual(Individual&& other);                            // ��������
		Individual& operator=(const Individual& other);            // ��ֵ����
		Individual& operator=(Individual&& other);                 // �ƶ���ֵ

		~Individual();                                             // ����
	};
}

#endif