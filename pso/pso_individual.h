#ifndef _PSO_INDIVIDUAL_H_
#define _PSO_INDIVIDUAL_H_

#include <initializer_list>

namespace opt
{
	// ���Ӹ���
	struct PSO_Individual
	{
		std::size_t nVars;                                                 // ��������
		double* xs;                                                        // ���Ӹ�������
		double* vs;                                                        // �����ٶ�����
		double* best_xs;                                                   // ���������������λ��
		double fitness;                                                    // ����Ի�������Ӧ��

		PSO_Individual();                                                  // ���캯��
		PSO_Individual(const std::size_t n);                               // ���캯��
		PSO_Individual(const std::initializer_list<double>& list);         // ���캯��
		PSO_Individual(const PSO_Individual& other);                       // ��������
		PSO_Individual(PSO_Individual&& other)noexcept;                    // �ƶ�����
		PSO_Individual& operator=(const PSO_Individual& other)noexcept;    // ��ֵ����
		PSO_Individual& operator=(PSO_Individual&& other)noexcept;         // �ƶ���ֵ

		~PSO_Individual();                                                 // ����
	};
}

#endif