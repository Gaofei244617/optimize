/* ָ���������������ɺ��� */
#ifndef _RANGE_RANDOM_
#define _RANGE_RANDOM_

namespace opt
{
	// ��������[min, max)�ڵ� double ���͵������
	double random_real(const double min, const double max);

	// ��������[min, max]�ڵ� int ���͵������
	int random_int(const int min, const int max);
}

#endif
