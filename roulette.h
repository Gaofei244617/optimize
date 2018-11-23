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
		std::size_t mark_num;                                   // 刻度线数量
		T* tick_mark;                                           // 轮盘刻度线(需保证数值单调递增)

	public:
		Roulette();                                             // 构造函数
		Roulette(const std::size_t N);                          // 构造函数
		Roulette(const Roulette<T>& other);                     // 拷贝构造
		Roulette(Roulette<T>&& other);                          // 移动构造
		Roulette<T>& operator=(const Roulette<T>& other);       // 赋值函数
		Roulette<T>& operator=(Roulette<T>&& other);            // 移动赋值
		~Roulette();                                            // 析构函数

		int roll();                                             // 轮盘随机转动
		T& operator[](const std::size_t& i);
		void reset(const std::size_t N);                        // 重新划分刻度线
	};

	// 构造函数
	template<class T>
	Roulette<T>::Roulette()
		:mark_num(0),
		tick_mark(nullptr)
	{}

	// 构造函数
	template<class T>
	Roulette<T>::Roulette(const std::size_t N)
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

	// 赋值函数
	template<class T>
	Roulette<T>& Roulette<T>::operator=(const Roulette<T>& other)
	{
		if (this != &other)
		{
			mark_num = other.mark_num;
			tick_mark = new T[mark_num]();
			for (int i = 0; i < mark_num; i++)
			{
				tick_mark[i] = other.tick_mark[i];
			}
		}
		return *this;
	}

	// 移动赋值
	template<class T>
	Roulette<T>& Roulette<T>::operator=(Roulette<T>&& other)
	{
		if (this != &other)
		{
			mark_num = other.mark_num;
			tick_mark = other.tick_mark;
			other.tick_mark = nullptr;
		}
		return *this;
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
				if (temp <= tick_mark[tempIndex])
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

	// 获取轮盘赌刻度线
	template<class T>
	T& Roulette<T>::operator[](const std::size_t& i)
	{
		return tick_mark[i];
	}

	//
	template<class T>
	void Roulette<T>::reset(const std::size_t N)
	{
		if (N != mark_num)
		{
			delete[] tick_mark;
			tick_mark = new T[N]();
			mark_num = N;
		}
	}
}

#endif