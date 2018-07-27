#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_

#include <string>
#include <memory>

namespace opt
{
	// ������
	struct Individual
	{
		//std::string groupName;                           // ��������Ⱥ
		const int nVars;                                   // ����(����)����
		double* vars;                                      // �����������(����)
		double fitness;                                    // ����Ի�������Ӧ��

		Individual(int n);                                 // ���캯��
		Individual(const Individual& other);               // ��������
		Individual& operator=(const Individual& other);    // ��ֵ����
		~Individual();                                     // ����
	};
}

#endif