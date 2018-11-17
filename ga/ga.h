#ifndef _Genetic_Algorithm_
#define _Genetic_Algorithm_

#include <string>
#include <chrono>
#include <tuple>
#include <vector>
#include <cmath>
#include <thread>
#include <functional>

#include "individual.h"
#include "range_random.h"
#include "roulette.h"
#include "opt_time.h"
#include "ga_GroupState.h"
#include "cross_factor.h"
#include "mutate_factor.h"
#include "ga_thread_sync.h"
#include "index_seq.h"
#include "ga_info.h"

/////////////////////////////////////////////////
#include <iostream>

namespace opt
{
	template<class F> class GAGroup;

	// GA��Ⱥ�࣬���ṩ��Ӧ�Ⱥ�������
	template<class R, class... Args>
	class GAGroup<R(Args...)>
	{
		using GenBound = std::initializer_list<std::initializer_list<double>>;
		friend class GAThreadSync<R, Args...>;

	private:
		std::string name;                                                     // ��Ⱥ����
		std::size_t groupSize;                                                // ��ʼ��Ⱥ����������default = 1000��
		const int nVars;                                                      // ��Ӧ�Ⱥ��������ı�������
		Individual* indivs;                                                   // ��Ⱥ����, ָ��groupSize����������(���һ��������Ÿ���)
		Individual* tempIndivs;                                               // �Ӵ����建����
		R(*fitFunc)(Args...);                                                 // ��Ӧ�Ⱥ���ָ��
		std::function<void(const GA_Info&)> monitor;                          // �ⲿ������
		double(*bound)[2];                                                    // ÿ������(����)������, ������ָ���ʾ
		Roulette<double> roulette;                                            // ���̶Ķ���
		double mutateProb;                                                    // ����������,Ĭ��p = 0.1
		double crossProb;                                                     // ���彻�����, Ĭ��p = 0.6
		std::vector<Individual> bestIndivs;                                   // ��¼ÿ�ε��������Ÿ���

		std::vector<std::thread> vec_run;                                     // �����߳���
		std::thread single_thread;                                            // ���߳�

		GA_GroupState group_state;                                            // ��Ⱥ״̬
		std::unique_ptr< GAThreadSync<R, Args...> > thread_sync;              // �߳�ͬ����

	public:
		GAGroup(R(*f)(Args...), const std::size_t size = 1000);               // ���캯��������һ����Ⱥ�����ṩ��Ӧ�Ⱥ�������Ⱥ����
		GAGroup(GAGroup<R(Args...)>& other);                                  // ��������
		GAGroup(GAGroup<R(Args...)>&& other);                                 // �ƶ�����
		GAGroup<R(Args...)>& operator=(const GAGroup<R(Args...)>& other) = delete;
		GAGroup<R(Args...)>& operator=(GAGroup<R(Args...)>&& other) = delete;
		~GAGroup();

		void setName(const std::string& str);                                 // ������Ⱥ����
		void setBoundary(double(*b)[2]);                                      // ���ñ�������
		void setBoundary(const GenBound& b);                                  // ���ñ�������
		void setMaxGeneration(const unsigned int N);                          // ��������������
		void setMaxRuntime(const Second& time);                               // �����������ʱ�䣨�룩
		void setStopTol(const long double t, const unsigned int N = 5);       // �������Ž�ֹͣ���
		void setMutateProb(const double p);                                   // ���û���������
		void setCrossProb(const double p);                                    // ���ý������
		void setThreadNum(const int NUM);                                     // ���ò��м�����߳�����Ĭ��Ϊ1
		void setMonitor(const std::function<void(const GA_Info&)>&);          // �����ⲿ��������

		const std::string getName()const;                                     // ��ȡ��Ⱥ����
		int getNVars()const;                                                  // ��ȡ��Ⱥ��������
		int getGeneration()const;                                             // ��õ�ǰ��Ⱥ����
		int getGroupSize()const;                                              // ��õ�ǰ��Ⱥ��������
		std::vector<Individual> getBestIndivs();                              // ��ȡ���ε��������Ž�
		int getStopCode();                                                    // ��ȡStop Code

		void initGroup(const std::vector<Individual>& indivs = std::vector<Individual>());    // ��ʼ����Ⱥ����
		bool start();                                                         // ��ʼ��������
		void pause();                                                         // ��ͣ(Ϊ��֤����һ���ԣ�����һ������������pause)
		void proceed();                                                       // ��������
		void kill();                                                          // ��������
		void wait_result();                                                   // ������ǰ�߳�,�ȴ��Ż����
		GAGroup<R(Args...)> clone();                                          // ��¡��ǰ��Ⱥ

	private:
		void crossover(const int seq);                                        // ����(���߳�ģʽ)
		void mutate(const int seq);                                           // ����(���߳�ģʽ)
		void select(const int seq);                                           // ѡ��(���߳�ģʽ)

		void run();                                                           // ���̵߳���
		void run_parallel(const int seq);                                     // ���е���

		void updateRoulette();                                                // �������̶Ŀ̶���
		void updateStopState();                                               // ������Ⱥֹͣ״̬
		void switchIndivArray();

		std::pair<int, int> selectPolarIndivs(const int seq, const int interval);   // Ѱ��������õĸ���λ��,����<worst, best>
		template<std::size_t... I>
		R callFitFunc(double* args, const opt::index_seq<I...>&);               // ������Ӧ�Ⱥ���
		bool flushStopFlag();                                                   // �ж��Ƿ��������
	};

	/******************************************* ���������� ***********************************************************/
	// ���캯�������ṩ����Ӧ�Ⱥ�������Ⱥ����
	template<class R, class... Args>
	GAGroup<R(Args...)>::GAGroup(R(*f)(Args...), const std::size_t size)
		:name("None"),
		groupSize(size + size % 2),
		nVars(sizeof...(Args)),
		indivs(nullptr),
		tempIndivs(nullptr),
		fitFunc(f),
		bound(nullptr),
		roulette(groupSize + 1),
		mutateProb(0.1),
		crossProb(0.6),
		group_state(),
		thread_sync(new GAThreadSync<R, Args...>(this))
	{
		// ��������ڴ棬�����������(Indivadual��Ĭ�Ϲ���),���һ��λ�����ڴ�����Ÿ���
		this->indivs = static_cast<Individual*>(::operator new(sizeof(Individual) * groupSize));
		this->tempIndivs = static_cast<Individual*>(::operator new(sizeof(Individual) * groupSize));

		// ���ѷ�����ڴ��Ϲ���������
		for (std::size_t i = 0; i < groupSize; i++)
		{
			new(indivs + i) Individual(nVars);
			new(tempIndivs + i) Individual(nVars);
		}

		// ��������
		bound = new double[nVars][2];
	}

	// ���ƹ���
	template<class R, class... Args>
	GAGroup<R(Args...)>::GAGroup(GAGroup<R(Args...)>& other)
		:name(other.name),
		groupSize(other.groupSize),
		nVars(other.nVars),
		indivs(nullptr),
		tempIndivs(nullptr),
		fitFunc(other.fitFunc),
		bound(nullptr),
		roulette(),
		mutateProb(other.mutateProb),
		crossProb(other.crossProb),
		bestIndivs(),
		group_state(),
		thread_sync(nullptr)
	{
		other.pause();

		// ��Ա����
		this->indivs = static_cast<Individual*>(::operator new(sizeof(Individual) * groupSize));
		this->tempIndivs = static_cast<Individual*>(::operator new(sizeof(Individual) * groupSize));

		for (int i = 0; i < groupSize; i++)
		{
			new(indivs + i) Individual(*(other.indivs + i));
			new(tempIndivs + i) Individual(*(other.tempIndivs + i));
		}

		// ��������
		bound = new double[nVars][2];
		for (int i = 0; i < nVars; i++)
		{
			bound[i][0] = (other.bound)[i][0];
			bound[i][1] = (other.bound)[i][1];
		}

		this->roulette = other.roulette;
		this->bestIndivs = other.bestIndivs;

		this->group_state = other.group_state;
		this->group_state.sleep.signal = false;
		this->group_state.sleep.result = false;
		this->group_state.stopCode = -1;   // δ��ʼ����

		this->thread_sync.reset(new GAThreadSync<R, Args...>(*(other.thread_sync), this));

		other.proceed();
	}

	// �ƶ�����
	template<class R, class... Args>
	GAGroup<R(Args...)>::GAGroup(GAGroup<R(Args...)>&& other)
		:name(other.name),
		groupSize(other.groupSize),
		nVars(other.nVars),
		indivs(other.indivs),
		tempIndivs(other.tempIndivs),
		fitFunc(other.fitFunc),
		bound(other.bound),
		roulette(),
		mutateProb(other.mutateProb),
		crossProb(other.crossProb),
		bestIndivs(),
		group_state(other.group_state),
		thread_sync(nullptr)
	{
		other.kill();        // ��ֹ��Ⱥ����

		other.indivs = nullptr;
		other.tempIndivs = nullptr;
		other.bound = nullptr;

		this->roulette = std::move(other.roulette);
		this->bestIndivs = std::move(other.bestIndivs);

		this->group_state = other.group_state;
		this->group_state.sleep.signal = false;
		this->group_state.sleep.result = false;
		this->group_state.stopCode = -1;

		thread_sync.reset(new GAThreadSync<R, Args...>(*(other.thread_sync), this));
	}

	// ��������
	template<class R, class... Args>
	GAGroup<R(Args...)>::~GAGroup()
	{
		// ����indivs����
		if (indivs != nullptr)
		{
			// ������Ⱥ�е�ÿһ�����壨����groupSize�����壩
			for (std::size_t i = 0; i < groupSize; i++)
			{
				(indivs + i) -> ~Individual();
			}
			// �ͷŸ������ָ��
			::operator delete(indivs);
		}

		// ����tempIndivs����
		if (tempIndivs != nullptr)
		{
			// ������Ⱥ�е�ÿһ�����壨����groupSize�����壩
			for (std::size_t i = 0; i < groupSize; i++)
			{
				(tempIndivs + i) -> ~Individual();
			}
			// �ͷŸ������ָ��
			::operator delete(tempIndivs);
		}

		// �ͷ�Boundary����ָ��
		delete[] bound;
	}
	/*****************************************************************************************************************/

	/**************************************** Setter ************************************************************/
	// ������Ⱥ����
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setName(const std::string& str)
	{
		this->name = str;
	}

	// �������б�������,�����ʼ���б�
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setBoundary(const GenBound& b)
	{
		const std::size_t len = b.size();
		for (std::size_t i = 0; i < len; i++)
		{
			bound[i][0] = *((*(b.begin() + i)).begin());         // �±߽�
			bound[i][1] = *((*(b.begin() + i)).begin() + 1);     // �ϱ߽�
		}
		group_state.setBoundFlag = true;
	}

	// �������б�������,��������ָ��
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setBoundary(double(*b)[2])
	{
		for (int i = 0; i < nVars; i++)
		{
			bound[i][0] = b[i][0];       // �±߽�
			bound[i][1] = b[i][1];       // �ϱ߽�
		}
		group_state.setBoundFlag = true;
	}

	// ��������������
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMaxGeneration(const unsigned int N)
	{
		group_state.setMaxGeneFlag = true;
		group_state.maxGene = N;
	}

	// �����������ʱ�䣨�룩
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMaxRuntime(const Second& time)
	{
		group_state.setRuntimeFlag = true;
		group_state.maxRuntime = time;
	}

	// �������Ž�ֹͣ���
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setStopTol(const long double t, const unsigned int N)
	{
		group_state.converCount = N;
		group_state.setStopTolFlag = true;
		group_state.stopTol = t;
	}

	// ���û���������
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMutateProb(const double p)
	{
		mutateProb = p;
	}

	// ���û��򽻲����
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setCrossProb(const double p)
	{
		crossProb = p;
	}

	// ���ò��м�����߳�����Ĭ��Ϊ1
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setThreadNum(const int NUM)
	{
		if (NUM > 1)
		{
			thread_sync->setThreadNum(NUM);
		}
		if (NUM < 1)
		{
			throw std::string("Thread number error.");
		}
	}

	// �����ⲿ��������
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMonitor(const std::function<void(const GA_Info&)>& func)
	{
		this->monitor = func;
	}

	/*****************************************************************************************************************/

	/**************************************Getter*********************************************************************/
	// ��ȡ��Ⱥ����
	template<class R, class... Args>
	const std::string GAGroup<R(Args...)>::getName()const
	{
		return this->name;
	}

	// ��ȡ��Ⱥ��������
	template<class R, class... Args>
	int GAGroup<R(Args...)>::getNVars()const
	{
		return nVars;
	}

	// ��ȡǰ��Ⱥ����
	template<class R, class... Args>
	int GAGroup<R(Args...)>::getGeneration()const
	{
		return group_state.nGene;
	}

	// ��ȡ��ǰ��Ⱥ��������
	template<class R, class... Args>
	int GAGroup<R(Args...)>::getGroupSize()const
	{
		return groupSize;
	}

	// ��ȡ���ε��������Ž�
	template<class R, class... Args>
	std::vector<Individual> GAGroup<R(Args...)>::getBestIndivs()
	{
		return bestIndivs;
	}

	// ��ȡֹͣ����
	template<class R, class... Args>
	int GAGroup<R(Args...)>::getStopCode()
	{
		return group_state.stopCode;
	}

	/******************************************************************************************************************/
	// ��ʼ����
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::start()
	{
		// ��ʼ����Ⱥ
		// �����Ⱥ�ѱ���ʼ�������ظ���ʼ��
		if (group_state.initFlag == false)
		{
			initGroup();
		}

		// ����Ⱥ�����������
		if (group_state.runable())
		{
			group_state.startTime = std::chrono::steady_clock::now();
			group_state.stopCode = 0; // ���ڵ���

			if (thread_sync->threadNum > 1) // ���߳�ģʽ
			{
				// ������Ⱥ�����߳���
				for (int i = 0; i < thread_sync->threadNum; i++)
				{
					vec_run.emplace_back(&GAGroup<R(Args...)>::run_parallel, this, i);
				}
			}
			else // ���߳�ģʽ
			{
				single_thread = std::thread(&GAGroup<R(Args...)>::run, this);
			}
			return true;
		}
		return false;
	}

	// ��������(���߳�ģʽ)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::run()
	{
		while (true)
		{
			// �Ƿ���ͣ����(Ϊ��֤����һ���ԣ�����һ������������pause)
			{
				std::unique_lock<std::mutex> lck(thread_sync->mtx);
				thread_sync->cv.wait(lck, [this]() { return !(this->group_state.sleep.signal); });
			}

			updateRoulette();                 // ���¸������̶Ŀ̶���
			crossover(0);                     // ����
			switchIndivArray();               // ������������
			mutate(0);                        // ����
			select(0);                        // ����ѡ��
			updateStopState();                // ����ֹͣ״̬

			// �����ⲿ��������
			if (monitor)
			{
				GA_Info ga_info(group_state.time, group_state.nGene, bestIndivs.back());
				monitor(ga_info);
			}

			// �ж��Ƿ񵽴�ֹͣ����
			if (flushStopFlag())
			{
				return;
			}

			if (group_state.sleep.signal == true)
			{
				group_state.sleep.result = true;
				thread_sync->cv.notify_all();
			}
		}
	}

	// ��������(���߳�ģʽ)
	template<class R, class ...Args>
	void GAGroup<R(Args...)>::run_parallel(const int seq)
	{
		while (true)
		{
			/////////////////////////////// Crossover ///////////////////////////////
			// Ϊ��֤����һ���ԣ�����һ������������pause
			{
				std::unique_lock<std::mutex> lck(thread_sync->mtx);
				thread_sync->cv.wait(lck, [this, &seq]() {
					return !(this->thread_sync->cross_flag[seq]) && this->thread_sync->crossReady && !(this->group_state.sleep.signal);
				});
			}

			if (group_state.stopFlag == true) { return; }

			crossover(seq);

			// �߳�ͬ��
			thread_sync->mtx.lock();
			thread_sync->cross_sync(seq);
			thread_sync->mtx.unlock();
			thread_sync->cv.notify_all();

			/////////////////////////////// Mutate ///////////////////////////////
			{
				std::unique_lock<std::mutex> lck(thread_sync->mtx);
				thread_sync->cv.wait(lck, [this, &seq]() {
					return !(this->thread_sync->mut_flag[seq]) && this->thread_sync->mutReady;
				});
			}

			mutate(seq);

			// �߳�ͬ��
			thread_sync->mtx.lock();
			thread_sync->mutate_sync(seq);
			thread_sync->mtx.unlock();
			thread_sync->cv.notify_all();

			/////////////////////////////// Select ///////////////////////////////
			{
				std::unique_lock<std::mutex> lck(thread_sync->mtx);
				thread_sync->cv.wait(lck, [this, &seq]() {
					return !(this->thread_sync->sel_flag[seq]) && this->thread_sync->selectReady;
				});
			}

			int worst_temp = 0;
			int best_temp = 0;

			int thread_num = thread_sync->threadNum;

			// Ѱ���Ӵ��������Ÿ���
			std::tie(worst_temp, best_temp) = this->selectPolarIndivs(seq, thread_num);

			// �߳�ͬ��
			thread_sync->mtx.lock();

			// ����������������λ��;
			if (indivs[worst_temp].fitness < indivs[group_state.worstIndex].fitness)
			{
				group_state.worstIndex = worst_temp;
			}
			if (indivs[best_temp].fitness > indivs[group_state.bestIndex].fitness)
			{
				group_state.bestIndex = best_temp;
			}

			thread_sync->select_sync(seq);
			thread_sync->mtx.unlock();
			thread_sync->cv.notify_all();
		}
	}

	// ������ǰ�߳�, �ȴ�������
	template<class R, class... Args>
	void GAGroup<R(Args...)>::wait_result()
	{
		if (thread_sync->threadNum == 1) // ���߳�ģʽ
		{
			if (single_thread.joinable())
			{
				single_thread.join();
			}
		}
		else  // ���߳�ģʽ
		{
			for (int i = 0; i < thread_sync->threadNum; i++)
			{
				if (vec_run[i].joinable())
				{
					vec_run[i].join();
				}
			}
		}
	}

	// ��¡��ǰ״̬��Ⱥ
	template<class R, class... Args>
	GAGroup<R(Args...)> GAGroup<R(Args...)>::clone()
	{
		return GAGroup<R(Args...)>(*this);
	}

	// ��ͣ(Ϊ��֤����һ���ԣ�����һ������������pause)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::pause()
	{
		// �����Ⱥδ���ڵ�������״̬
		if (group_state.stopCode != 0)
		{
			return;
		}

		this->group_state.sleep.signal = true;
		std::unique_lock<std::mutex> lck(thread_sync->mtx);
		thread_sync->cv.wait(lck, [this]() {return this->group_state.sleep.result; });
	}

	// ��������
	template<class R, class... Args>
	void GAGroup<R(Args...)>::proceed()
	{
		// �����Ⱥδ���ڵ�������״̬
		if (group_state.stopCode != 0)
		{
			return;
		}

		this->group_state.sleep.signal = false;
		this->group_state.sleep.result = false;
		this->thread_sync->cv.notify_all();
	}

	// ֹͣ����
	template<class R, class... Args>
	void GAGroup<R(Args...)>::kill()
	{
		// �����Ⱥδ���ڵ�������״̬
		if (group_state.stopCode != 0)
		{
			return;
		}

		this->group_state.stopCode = 4;
		this->wait_result();
	}

	/******************************************* Private Functions *****************************************************/
	// ��ʼ����Ⱥ
	template<class R, class... Args>
	void GAGroup<R(Args...)>::initGroup(const std::vector<Individual>& init_indivs)
	{
		// ��ʼ��ǰ��Ҫ��֤�����û����������
		if (group_state.setBoundFlag)
		{
			// 1.��ʼ��ÿһ������
			for (std::size_t i = 0; i < init_indivs.size(); i++)
			{
				indivs[i] = init_indivs[i];
			}
			for (std::size_t i = init_indivs.size(); i < groupSize; i++)
			{
				// ���ø����ʼ����
				for (int j = 0; j < nVars; j++)
				{
					indivs[i].vars[j] = random_real(bound[j][0], bound[j][1]);
				}
				// ������Ӧ��
				indivs[i].fitness = callFitFunc(indivs[i].vars, opt::make_index_seq<sizeof...(Args)>());
			}

			// 2.Ѱ�����������Ӧ�ȸ���
			std::tie(group_state.worstIndex, group_state.bestIndex) = this->selectPolarIndivs(0, 1);
			bestIndivs.push_back(indivs[group_state.bestIndex]);                     // ��¼���Ž�

			// 3.�������̶Ŀ̶���
			updateRoulette();

			// 4.���ó�ʼ����־λ
			group_state.initFlag = true;

			// 5.��ʼ���߳�ͬ����
			thread_sync->crossReady = true;
			thread_sync->mutReady = false;
			thread_sync->selectReady = false;
		}
		else
		{
			throw std::string("Fitness function or genes boundary is not set.");
		}
	}

	// ����(���߳�ģʽ)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::crossover(const int seq)
	{
		int Index_M = 0;             // ������
		int Index_F = 0;             // ĸ����
		double rand_cross = 0;       // ���������

		int thread_num = thread_sync->threadNum;          // �߳���

		// ���ѡȡ���彻��, fitnessԽ�󣬱�ѡ�еĸ���Խ��
		for (std::size_t i = seq * 2; i < groupSize; i += 2 * thread_num)
		{
			// ���ѡȡ����������Ϊ����
			Index_M = roulette.roll();
			Index_F = roulette.roll();

			// �����Ӵ�����, ����ڻ�������tempIndivs��
			for (int j = 0; j < nVars; j++)
			{
				// ���ݽ�����ʾ��������Ƿ���н���
				rand_cross = random_real(0, 1);
				if (rand_cross <= crossProb)   // ����
				{
					// ���򽻲�����Ӵ�����
					std::tie(tempIndivs[i].vars[j], tempIndivs[i + 1].vars[j]) = cross_SBX(indivs[Index_M].vars[j], indivs[Index_F].vars[j]);

					// �Ӵ������Ƿ񳬹��߽�
					if (tempIndivs[i].vars[j] < bound[j][0] || tempIndivs[i].vars[j] > bound[j][1])
					{
						tempIndivs[i].vars[j] = random_real(bound[j][0], bound[j][1]);
					}
					if (tempIndivs[i + 1].vars[j] < bound[j][0] || tempIndivs[i + 1].vars[j] > bound[j][1])
					{
						tempIndivs[i + 1].vars[j] = random_real(bound[j][0], bound[j][1]);
					}
				}
				else  // ������
				{
					// �̳и������򣬻��򲻷�������
					tempIndivs[i].vars[j] = indivs[Index_M].vars[j];
					tempIndivs[i + 1].vars[j] = indivs[Index_F].vars[j];
				}
			}
		}
	}

	// ����(���߳�ģʽ)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::mutate(const int seq)
	{
		double rand_num = 0; // �����

		int thread_num = thread_sync->threadNum;  // �߳���

		// ÿ������ı���
		for (std::size_t i = seq; i < groupSize; i += thread_num)
		{
			// ÿ������ı���
			for (int j = 0; j < nVars; j++)
			{
				rand_num = random_real(0, 1);
				if (rand_num < mutateProb)
				{
					// �����������
					indivs[i].vars[j] = mutate_PM(indivs[i].vars[j], bound[j][0], bound[j][1]);
				}
			}

			// ����ÿ��������Ӧ��
			indivs[i].fitness = callFitFunc(indivs[i].vars, opt::make_index_seq<sizeof...(Args)>());
		}
	}

	// ѡ��(���߳�ģʽ)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::select(const int seq)
	{
		int thread_num = thread_sync->threadNum; // �߳���

		// Ѱ���Ӵ��������Ÿ���
		std::tie(group_state.worstIndex, group_state.bestIndex) = this->selectPolarIndivs(seq, thread_num);

		// ��̭�Ӵ�������, �ø������Ÿ���ȡ��
		indivs[group_state.worstIndex] = bestIndivs.back();

		// ����Ӵ����ָ��Ÿ���
		if (indivs[group_state.bestIndex].fitness >= bestIndivs.back().fitness)
		{
			bestIndivs.push_back(indivs[group_state.bestIndex]);
		}
		else
		{
			bestIndivs.push_back(bestIndivs.back());
		}
	}

	// �������̶Ŀ̶���
	template<class R, class... Args>
	void GAGroup<R(Args...)>::updateRoulette()
	{
		int worst = 0;

		// ����indivs��fitness��Ѱ��fitness����Сֵ,��¼��Сֵλ��
		for (std::size_t i = 1; i < groupSize; i++)
		{
			if (indivs[i].fitness < indivs[worst].fitness)
			{
				worst = i;
			}
		}

		// �������̶Ŀ̶���
		for (std::size_t i = 1; i < groupSize + 1; i++)
		{
			roulette[i] = roulette[i - 1] + (indivs[i - 1].fitness - indivs[worst].fitness);
		}
	}

	// ������Ⱥֹͣ״̬
	template<class R, class... Args>
	void GAGroup<R(Args...)>::updateStopState()
	{
		// ����������һ
		group_state.nGene++;

		// ��¼��ǰʱ��
		group_state.nowTime = std::chrono::steady_clock::now();

		// �ж����Ÿ���fitnessֵ����һ���Ĳ������
		if (abs(bestIndivs[group_state.nGene - 1].fitness - bestIndivs[group_state.nGene].fitness) <= group_state.stopTol)
		{
			group_state.count++;
		}
		else
		{
			group_state.count = 0;
		}
	}

	// ���¸�������ָ��
	template<class R, class... Args>
	void GAGroup<R(Args...)>::switchIndivArray()
	{
		// �����Ӵ����建�����͸�������������ָ��
		Individual* temp = indivs;
		indivs = tempIndivs;
		tempIndivs = temp;
	}

	// Ѱ��������õĸ���λ��, �β�Ϊ�߳�ID(0 ~ N)�������߳�����, ����pair<worst, best>
	template<class R, class... Args>
	std::pair<int, int> GAGroup<R(Args...)>::selectPolarIndivs(const int seq, const int interval)
	{
		int worst = seq;
		int best = seq;

		// ����indivs��fitness��Ѱ��fitness�ļ�ֵ����¼��ֵλ��
		for (std::size_t i = seq + interval; i < groupSize; i += interval)
		{
			if (indivs[i].fitness < indivs[worst].fitness)
			{
				worst = i;
			}

			if (indivs[i].fitness > indivs[best].fitness)
			{
				best = i;
			}
		}

		return std::make_pair(worst, best);
	}

	// ������Ӧ�Ⱥ���
	template<class R, class... Args>
	template<std::size_t... I>
	inline R GAGroup<R(Args...)>::callFitFunc(double* args, const opt::index_seq<I...>&)
	{
		return fitFunc(args[I]...);
	}

	// �ж��Ƿ��������
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::flushStopFlag()
	{
		// Stop Code : -1-δ��ʼ����; 0-���ڵ���; 1-�ﵽ����������; 2-�ﵽ������ʱ��; 3-���Ž��������ȶ�ֵ; 4-��Ϊֹͣ����

		// �Ƿ�ﵽ����������
		if (group_state.setMaxGeneFlag && group_state.nGene >= group_state.maxGene)
		{
			group_state.stopCode = 1;
		}

		// �Ƿ�ﵽ������ʱ��
		std::chrono::duration<double> evolTime = group_state.nowTime - group_state.startTime;
		group_state.time = evolTime.count(); // ��
		if (group_state.setRuntimeFlag && group_state.time >= group_state.maxRuntime.value)
		{
			group_state.stopCode = 2;
		}

		// �Ƿ�����5�����Ž�Ĳ���С��ֹͣ���
		if (group_state.setStopTolFlag && group_state.count == group_state.converCount)
		{
			group_state.stopCode = 3;
		}

		// ����stop code�ж���Ⱥ�Ƿ�ֹͣ����
		if (group_state.stopCode != -1 && group_state.stopCode != 0)
		{
			group_state.stopFlag = true;
		}

		return group_state.stopFlag;
	}
}

#endif