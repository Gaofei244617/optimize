#ifndef _GA_THREAD_SYNC_
#define _GA_THREAD_SYNC_

#include <thread>
#include <vector>
#include <utility>
#include <condition_variable>
#include <tuple>
#include "bool_array.h"

namespace opt
{
	template<class F> class GAGroup;

	// GA�߳�ͬ����
	template<class R, class... Args>
	class GAThreadSync
	{
		friend class GAGroup<R(Args...)>;

	private:
		std::mutex mtx;                                   // ������
		std::condition_variable cv;                       // ��������, �����߳�ͬ��
		bool crossReady;                                  // ���������־
		bool mutReady;                                    // ���������־
		bool selectReady;                                 // ѡ�������־
		bool sleep;                                       // 
		int threadNum;                                    // ���м���ʹ�õ��߳���
		bool_array cross_flag;                            // ���м���: �����߳�״̬��־
		bool_array mut_flag;                              // ���м���: �����߳�״̬��־
		bool_array sel_flag;                              // ���м���: ѡ���߳�״̬��־
		GAGroup<R(Args...)>* ga;                          // ͬ����������GA��Ⱥ

	public:
		// ���湹��
		GAThreadSync(GAGroup<R(Args...)>* ptr) 
			:crossReady(false),
			mutReady(false),
			selectReady(false),
			sleep(false),
			threadNum(1),
			cross_flag(threadNum),
			mut_flag(threadNum),
			sel_flag(threadNum),
			ga(ptr)
		{}

		// ���ƹ���
		GAThreadSync(const GAThreadSync& other, GAGroup<R(Args...)>* ptr)
			:crossReady(other.crossReady),
			mutReady(other.mutReady),
			selectReady(other.selectReady),
			sleep(other.sleep),
			threadNum(other.threadNum),
			cross_flag(other.cross_flag),
			mut_flag(other.mut_flag),
			sel_flag(other.sel_flag),
			ga(ptr)
		{}

		void setThreadNum(const int N);	            // �����߳�����

		void cross_sync(const int thread_seq);		// "����"�߳�ͬ��
		void mutate_sync(const int thread_seq);		// "����"�߳�ͬ��
		void select_sync(const int thread_seq);		// "ѡ��"�߳�ͬ��
	};
    
	///////////////////////////////////////////////////////////////////////////////////////////
	// �����߳�����
	template<class R, class... Args>
	void GAThreadSync<R, Args...>::setThreadNum(const int N)
	{
		threadNum = N;
		cross_flag.set_length(N);
		mut_flag.set_length(N);
		sel_flag.set_length(N);
	}

	// �����桱�߳�ͬ��
	template<class R, class... Args>
	void GAThreadSync<R, Args...>::cross_sync(const int thread_seq)
	{
		cross_flag[thread_seq] = true;
		selectReady = false;
		
		// �������߳���ȫ��������
		if (cross_flag.is_all_true())
		{
			// ������������
			ga->switchIndivArray();

			cross_flag.set_all(false);
			crossReady = false;
			mutReady = true;
		}
		else
		{
			mutReady = false;
		}
	}

	// "����"�߳�ͬ��
	template<class R, class... Args>
	void GAThreadSync<R, Args...>::mutate_sync(const int thread_seq)
	{
		mut_flag[thread_seq] = true;
		crossReady = false;
		
		// �������߳���ȫ��������
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

	// "ѡ��"�߳�ͬ��
	template<class R, class... Args>
	void GAThreadSync<R, Args...>::select_sync(const int thread_seq)
	{
		sel_flag[thread_seq] = true;
		mutReady = false;
		
		// ��ѡ���߳���ȫ��������
		if (sel_flag.is_all_true())
		{
			// ��̭�Ӵ�������, �ø������Ÿ���ȡ��
			ga->indivs[ga->group_state.worstIndex] = ga->bestIndivs.back();

			// ����������ָ��Ÿ���
			if (ga->indivs[ga->group_state.bestIndex].fitness >= ga->bestIndivs.back().fitness)
			{
				ga->bestIndivs.push_back(ga->indivs[ga->group_state.bestIndex]);
			}
			else
			{
				ga->bestIndivs.push_back(ga->bestIndivs.back());
			}
			
			// �������̶Ŀ̶���
			ga->updateRoulette();

			// ����GA��Ⱥֹͣ״̬
			ga->updateStopState();

			ga->flushStopFlag();

			sel_flag.set_all(false);
			selectReady = false;
			crossReady = true;
		}
		else
		{
			crossReady = false;
		}
	}
}

#endif
