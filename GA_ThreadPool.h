#ifndef _GA_THREADPOOL_H_
#define _GA_THREADPOOL_H_

#include <thread>
#include <vector>
#include <utility>
#include <condition_variable>
#include <atomic>

namespace opt
{
	/******************************************************************************************/
	// bool型数组容器
	class bool_array
	{
	private:
		bool* array;
		int Length;
		int Capacity;

	public:
		// 构造函数，所有元素初始化为false
		bool_array():array(nullptr),Length(0),Capacity(0){}

		bool_array(const int N)
			:array(new bool[N]),
			Length(N),
			Capacity(N)
		{
			for (int i = 0; i < N; i++)
			{
				array[i] = false;
			}
		}

		// 复制构造
		bool_array(const bool_array& other)
			:array(new bool[other.Length]),
			Length(other.Length)
		{
			for (int i = 0; i < Length; i++)
			{
				array[i] = (other.array)[i];
			}
		}

		// 判断是否所有元素为true
		bool is_all_true()const
		{
			if (Length > 0)
			{
				bool temp = true;
				for (int i = 0; i < Length; i++)
				{
					temp = temp && array[i];
				}
				return temp;
			}
			return false;
		}

		// 为所有元素设置相同值
		void set_all(const bool b)
		{
			for (int i = 0; i < Length; i++)
			{
				array[i] = b;
			}
		}
		
		// 设置bool容器长度
		void set_length(const int N)
		{
			/*** 分三种情况：(1) N > Capacity; (2.1) Length < N <= Capacity; (2.2) N <= Length ***/
			// (1) N > Capacity的情况,需要需要重新分配内存
			if (N > Capacity)
			{
				bool* temp = new bool[N];
				for (int i = 0; i < Capacity; i++)
				{
					temp[i] = array[i];
				}
				for (int i = Capacity; i < N; i++)
				{
					temp[i] = false;
				}

				Capacity = N;
				Length = N;

				bool* temp2 = temp;
				temp = array;
				array = temp2;
				delete[] temp;
				temp = temp2 = nullptr;
			}
			else
			{
				// (2.1) Length < N <= Capacity; (2.2) N <= Length, 不需要重新分配内存
				if (N > Length)
				{
					for (int i = Length; i < N; i++)
					{
						array[i] = false;
					}
				}
				Length = N;				
			}
		}

		// 下标索引运算符
		bool& operator[](const int n)
		{
			return array[n];
		}

		// 析构函数
		~bool_array()
		{
			delete[] array;
		}
	};

    /*******************************************************************************************/
	// GA算法并行计算相关状态标志
	struct GA_ThreadState
	{
		std::mutex mtx;  // 互斥量
		std::condition_variable cv; // 条件变量, 用于线程同步
		bool crossReady;  // 交叉操作标志
		bool mutReady;    // 变异操作标志
		bool selectReady; // 选择操作标志
		int threadNum;    // 并行计算使用的线程数
		bool_array cross_flag; // 并行计算: 交叉线程状态标志
		bool_array mut_flag;   // 并行计算: 变异线程状态标志

		// 构造函数
		GA_ThreadState()
			:crossReady(false),
			mutReady(false),
			selectReady(false),
			threadNum(1),
			cross_flag(threadNum),
			mut_flag(threadNum)
		{}

		// 复制构造
		GA_ThreadState(const GA_ThreadState& other)
			:crossReady(other.crossReady),
			mutReady(other.mutReady),
			selectReady(other.selectReady),
			threadNum(other.threadNum),
			cross_flag(other.cross_flag),
			mut_flag(other.mut_flag)
		{}
	};
}

#endif
