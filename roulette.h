#ifndef _ROULETTE_H_
#define _ROULETTE_H_

#include "range_random.h"

namespace opt
{
	// 轮盘赌封装类
	template<class T>
	class Roulette
	{
	private:
		int mark_num;                                           // 刻度线数量
		T* tick_mark;                                           // 轮盘刻度线(需保证数值单调递增)

	public:
		Roulette();                                             // 构造函数 
		Roulette(const int N);                                  // 构造函数 
		Roulette(const Roulette<T>& other);                     // 拷贝构造
		Roulette(Roulette<T>&& other);                          // 移动构造
		~Roulette();                                            // 析构函数

		int roll();                                             // 轮盘随机转动
		void resetMark(const int INDEX, const T& NUM);          // 重设刻度线的值
	};

	// 构造函数
	template<class T>
	Roulette<T>::Roulette() 
		:mark_num(0),
		tick_mark(nullptr)
	{}
	
	// 构造函数
	template<class T>
	Roulette<T>::Roulette(const int N)
		: mark_num(N),
		tick_mark(new T[N]())
	{}

	// 复制构造
	template<class T>
	Roulette<T>::Roulette(const Roulette<T>& other)
		:mark_num(other.mark_num),
		tick_mark(new T[mark_num]())
	{
		for (int i = 0; i < mark_num; i++)
		{
			tick_mark[i] = other.tick_mark[i];
		}
	}

	// 移动构造
	template<class T>
	Roulette<T>::Roulette(Roulette<T>&& other)
		:mark_num(other.mark_num),
		tick_mark(other.tick_mark)
	{
		other.tick_mark = nullptr;
	}

	// 析构函数
	template<class T>
	Roulette<T>::~Roulette()
	{
		delete[] tick_mark;
	}

	// 轮盘随机转动指向一个区间
	template<class T>
	int Roulette<T>::roll()
	{
		if (mark_num > 0)
		{
			double temp = random_real(tick_mark[0], tick_mark[mark_num - 1]);

			int lowIndex = 0;
			int upIndex = mark_num - 1;
			int tempIndex;

			// 二分法查找,时间复杂度log_n
			// 通过逐步移动lowIndex和upIndex代表的下标位置，每次将搜索区间减小一半
			while (upIndex - lowIndex > 1)
			{
				tempIndex = (lowIndex + upIndex) / 2 + (lowIndex + upIndex) % 2;
				if (temp <= fitArrayCache[tempIndex])
				{
					upIndex = tempIndex;
				}
				else
				{
					lowIndex = tempIndex;
				}
			}
			return lowIndex;  // up = 1, low = 0 时随机个体是indivs[0]
		}
		else
		{
			return 0;
		}
	}

	// 重设刻度线的值
	template<class T>
	inline void Roulette<T>::resetMark(const int INDEX, const T& NUM)
	{
		tick_mark[INDEX] = NUM;
	}
}

#endif