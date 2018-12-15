#ifndef _GA_INDIVIDUAL_H_
#define _GA_INDIVIDUAL_H_

#include <initializer_list>

namespace opt
{
	// ������
	struct GA_Individual
	{
		std::size_t nVars;                                               // ����(����)����
		double* vars;                                                    // �����������(����)
		double fitness;                                                  // ����Ի�������Ӧ��

		GA_Individual();                                                 // ���캯��
		GA_Individual(const std::size_t n);                              // ���캯��
		GA_Individual(const std::initializer_list<double>& list);        // ���캯��
		GA_Individual(const GA_Individual& other);                       // ��������
		GA_Individual(GA_Individual&& other)noexcept;                    // �ƶ�����
		GA_Individual& operator=(const GA_Individual& other)noexcept;    // ��ֵ����
		GA_Individual& operator=(GA_Individual&& other)noexcept;         // �ƶ���ֵ

		~GA_Individual();                                                // ����
	};
}

#endif