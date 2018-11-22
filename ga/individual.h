#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

#include <initializer_list>

namespace opt
{
	// ������
	struct Individual
	{
		std::size_t nVars;                                   // ����(����)����
		double* vars;                                              // �����������(����)
		double fitness;                                            // ����Ի�������Ӧ��

		Individual();                                              // ���캯��
		Individual(const std::size_t n);                           // ���캯��
		Individual(const std::initializer_list<double>& list);     // ���캯��
		Individual(const Individual& other);                       // ��������
		Individual(Individual&& other)noexcept;                    // �ƶ�����
		Individual& operator=(const Individual& other)noexcept;    // ��ֵ����
		Individual& operator=(Individual&& other)noexcept;         // �ƶ���ֵ

		~Individual();                                             // ����
	};
}

#endif