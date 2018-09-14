#ifndef _GA_THREAD_SYNC_MACHINE_
#define _GA_THREAD_SYNC_MACHINE_

#include <thread>
#include <vector>
#include <utility>
#include <condition_variable>
#include "bool_array.h"
#include <iostream>

namespace opt
{
	template<class F>
	class GAGroup;

	// GA�߳�ͬ����
	template<class R, class... Args>
	class GAThreadSyncMachine
	{
	private:
		friend class GAGroup<R(Args...)>;

		std::mutex mtx;  // ������
		std::condition_variable cv; // ��������, �����߳�ͬ��
		bool crossReady;  // ���������־
		bool mutReady;    // ���������־
		bool selectReady; // ѡ�������־
		int threadNum;    // ���м���ʹ�õ��߳���
		bool_array cross_flag; // ���м���: �����߳�״̬��־
		bool_array mut_flag;   // ���м���: �����߳�״̬��־
		GAGroup<R(Args...)>* ga;

	public:
		// ���湹��
		GAThreadSyncMachine(GAGroup<R(Args...)>* ptr) 
			:crossReady(false),
			mutReady(false),
			selectReady(false),
			threadNum(1),
			cross_flag(threadNum),
			mut_flag(threadNum),
			ga(ptr)
		{}

		// ���ƹ���
		GAThreadSyncMachine(const GAThreadSyncMachine& other, GAGroup<R(Args...)>* ptr)
			:crossReady(other.crossReady),
			mutReady(other.mutReady),
			selectReady(other.selectReady),
			threadNum(other.threadNum),
			cross_flag(other.cross_flag),
			mut_flag(other.mut_flag),
			ga(ptr)
		{}

		// �����桱�߳�ͬ��
		void cross_sync(const int thread_seq)
		{
			cross_flag[thread_seq] = true;
			selectReady = false;
			if (cross_flag.is_all_true())
			{
				cross_flag.set_all(false);
				crossReady = false;
				mutReady = true;
			}
			else
			{
				mutReady = false;
			}
		}

		// �����족�߳�ͬ��
		void mutate_sync(const int thread_seq)
		{
			mut_flag[thread_seq] = true;
			crossReady = false;
			if (mut_flag.is_all_true())
			{
				mut_flag.set_all(false);
				mutReady = false;
				selectReady = true;
			}
			else
			{
				selectReady = false;
			}
		}

		// ��ѡ���߳�ͬ��
		void select_sync()
		{
			selectReady = false;
			crossReady = true;
			mutReady = false;
		}
	};
}

#endif
