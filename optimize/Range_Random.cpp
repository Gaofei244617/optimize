#include "range_random.h"
#include <random>
#include <ctime>

namespace opt
{
	// ��������棬��������һ�������unsigned����
	// ���������������()�����
	static std::default_random_engine rand_eng((unsigned int)time(NULL));

	// ��������[min, max)�ڵ� double ���͵������
	double random_real(const double min, const double max)
	{
		return std::uniform_real_distribution<double>(min, max)(rand_eng);
	}

	// ��������[min, max]�ڵ� int ���͵������
	int random_int(const int min, const int max)
	{
		return std::uniform_int_distribution<int>(min, max)(rand_eng);
	}
}