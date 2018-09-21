#ifndef _Genetic_Algorithm_
#define _Genetic_Algorithm_

#include <string>
#include <utility>
#include <chrono>
#include <tuple>
#include <vector>
#include <cmath>
#include <thread>
#include "Individual.h"
#include "Range_Random.h"
#include "GA_Time.h"
#include "GA_GroupState.h"
#include "cross_factor.h"
#include "mutate_factor.h"
#include "ga_thread_sync.h"

/*************Test*******************/
#include <iostream>
#include "out_bool.h"
#include <iomanip>
//#include <crtdbg.h>
/************************************/

namespace opt
{
	template<class F> class GAGroup;

	// ��Ⱥ�࣬���ṩ��Ӧ�Ⱥ�������
	template<class R, class... Args>
	class GAGroup<R(Args...)>
	{
		using GenBound = std::initializer_list<std::initializer_list<double>>;
		friend void GAThreadSync<R, Args...>::select_sync(const int thread_seq);

	private:
		std::string name;                                                     // ��Ⱥ����
		int groupSize;                                                        // ��ʼ��Ⱥ����������default = 1000��
		const int nVars;                                                      // ��Ӧ�Ⱥ��������ı�������    
		Individual* indivs;                                                   // ��Ⱥ����, ָ��groupSize + 1����������(���һ��������Ÿ���)
		Individual* tempIndivs;                                               // �Ӵ����建����
		R(*fitFunc)(Args...);                                                 // ��Ӧ�Ⱥ���ָ��
		double(*bound)[2];                                                    // ÿ������(����)������, ������ָ���ʾ
		double* fitArrayCache;                                                // ���̶Ŀ̶���: fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		double mutateProb;                                                    // ����������,Ĭ��p = 0.1
		double crossProb;                                                     // ���彻�����, Ĭ��p = 0.6
		std::vector<Individual> bestIndivs;                                   // ��¼ÿ�ε��������Ÿ���	

		std::vector<std::thread> vec_cross;                                   // �����߳���
		std::vector<std::thread> vec_mut;                                     // �����߳���
		std::vector<std::thread> vec_select;                                  // ����ѡ���߳���
		std::thread single_thread;                                            // ���߳�
		
		GA_GroupState group_state;                                            // ��Ⱥ״̬
		GAThreadSync<R, Args...> thread_sync;                                 // �߳�ͬ����  

	public:
		GAGroup(R(*f)(Args...), const int size = 1000);                       // ���캯��������һ����Ⱥ�����ṩ��Ӧ�Ⱥ�������Ⱥ����
		GAGroup(const GAGroup<R(Args...)>& other);                            // ��������
		GAGroup(GAGroup<R(Args...)>&& other);                                 // �ƶ�����
		~GAGroup();

		void setName(const std::string& str);                                 // ������Ⱥ����
		void setBoundary(double(*b)[2]);                                      // ���ñ�������
		void setBoundary(const GenBound& b);                                  // ���ñ�������
		void setMaxGeneration(const unsigned int N);                          // ��������������
		void setMaxRuntime(const Second& time);                               // �����������ʱ�䣨�룩
		void setStopTol(long double t, unsigned int N = 5);                   // �������Ž�ֹͣ���
		void setMutateProb(double p);                                         // ���û���������
		void setCrossProb(double p);                                          // ���ý������
		void setThreadNum(const int NUM);                                     // ���ò��м�����߳�����Ĭ��Ϊ1

		const std::string getName()const;                                     // ��ȡ��Ⱥ����
		int getNVars()const;                                                  // ��ȡ��Ⱥ�������� 
		int getGeneration()const;                                             // ��õ�ǰ��Ⱥ����
		int getGroupSize()const;                                              // ��õ�ǰ��Ⱥ��������	
		const Individual& getIndivByIndex(const int index)const;              // ��ȡ��Ⱥ���壨ͨ���±꣩

		void test();                                                       ///////////////////////////////

		bool start();                                                         // ��ʼ����
		void wait_result();                                                   // �ȴ��������
		bool pause();                                                         // ֹͣ����
		bool proceed();                                                       // ��������                                                    

	private:
		void initGroup();                                                     // ��ʼ����Ⱥ����

		void crossover();                                                     // ����(���߳�ģʽ)
		void mutate();                                                        // ����(���߳�ģʽ)
		void select();                                                        // ѡ��(���߳�ģʽ)
		void run();                                                           // ��������(���߳�ģʽ)
		
		void cross_thread(const int seq);                                     // ����(���߳�ģʽ)
		void mut_thread(const int seq);                                       // ����(���߳�ģʽ)
		void sel_thread(const int seq);                                       // ѡ��(���߳�ģʽ)
		
		void updateFitArrayCache();                                           // �������̶Ŀ̶���
		void updateStopState();                                               // ������Ⱥֹͣ״̬

		std::pair<int, int> selectPolarIndivs(const int seq, const int interval);   // Ѱ��������õĸ���λ��,����<worst, best>
		int randomPickIndiv();                                                // ���ݸ����fitness���ѡȡһ�����壬���ظ���λ��
		template<std::size_t... I>
		R callFitFunc(double* args, const std::index_sequence<I...>&);        // ������Ӧ�Ⱥ���
		bool getStopFlag();                                                   // �ж��Ƿ��������
	};

	/******************************************* Helper Functions*****************************************************/
	template<class R, class... Args>
	GAGroup<R(Args...)> createGAGroup(R(*func)(Args...), const int N = 1000)
	{
		return GAGroup<R(Args...)>(func, N);
	}
	/*****************************************************************************************************************/

	/******************************************* ���������� ***********************************************************/
	// ���캯�������ṩ����Ӧ�Ⱥ�������Ⱥ����
	template<class R, class... Args>
	GAGroup<R(Args...)>::GAGroup(R(*f)(Args...), const int size)
		:name("None"),
		groupSize(size + size % 2),
		nVars(sizeof...(Args)),
		fitFunc(f),
		indivs(nullptr),
		tempIndivs(nullptr),
		bound(nullptr),
		fitArrayCache(nullptr),
		mutateProb(0.1),
		crossProb(0.6),
		thread_sync(this)
	{
		// ��������ڴ棬�����������(Indivadual��Ĭ�Ϲ���),���һ��λ�����ڴ�����Ÿ���
		this->indivs = static_cast<Individual*>(::operator new(sizeof(Individual) * (groupSize + 1)));
		this->tempIndivs = static_cast<Individual*>(::operator new(sizeof(Individual) * (groupSize + 1)));

		// ���ѷ�����ڴ��Ϲ���������
		for (int i = 0; i < groupSize + 1; i++)
		{
			new(indivs + i) Individual(nVars);
			new(tempIndivs + i) Individual(nVars);
		}

		// ��������
		bound = new double[nVars][2];

		// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		fitArrayCache = new double[groupSize + 1]();
	}

	// ���ƹ���
	template<class R, class... Args>
	GAGroup<R(Args...)>::GAGroup(const GAGroup<R(Args...)>& other)
		:name(other.name),
		groupSize(other.groupSize),
		nVars(other.nVars),
		fitFunc(other.fitFunc),
		mutateProb(other.mutateProb),
		crossProb(other.crossProb),
		indivs(nullptr),
		tempIndivs(nullptr),
		bound(nullptr),
		fitArrayCache(nullptr),
		group_state(other.group_state),
		thread_sync(other.thread_sync, this)
	{
		// ��������ڴ棬�����������(Indivadual��Ĭ�Ϲ���)
		this->indivs = static_cast<Individual*>(::operator new(sizeof(Individual) * (groupSize + 1)));
		this->tempIndivs = static_cast<Individual*>(::operator new(sizeof(Individual) * (groupSize + 1)));

		// ���ѷ�����ڴ��Ϲ���������
		for (int i = 0; i < groupSize + 1; i++)
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

		// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		fitArrayCache = new double[groupSize + 1]();
		for (int i = 0; i < groupSize + 1; i++)
		{
			fitArrayCache[i] = (other.fitArrayCache)[i];
		}
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
		fitArrayCache(other.fitArrayCache),
		mutateProb(other.mutateProb),
		crossProb(other.crossProb),
		group_state(other.group_state)
	{
		other.indivs = nullptr;
		other.tempIndivs = nullptr;
		other.bound = nullptr;
		other.fitArrayCache = nullptr;
	}

	// ��������
	template<class R, class... Args>
	GAGroup<R(Args...)>::~GAGroup()
	{
		// ����indivs����
		if (indivs != nullptr)
		{
			// ������Ⱥ�е�ÿһ�����壨����groupSize + 1�����壩
			for (int i = 0; i < groupSize + 1; i++)
			{
				(indivs + i) -> ~Individual();
			}
			// �ͷŸ������ָ��
			::operator delete(indivs);
		}

		// ����tempIndivs����
		if (tempIndivs != nullptr)
		{
			// ������Ⱥ�е�ÿһ�����壨����groupSize + 1�����壩
			for (int i = 0; i < groupSize + 1; i++)
			{
				(tempIndivs + i) -> ~Individual();
			}
			// �ͷŸ������ָ��
			::operator delete(tempIndivs);
		}

		// �ͷ�Boundary����ָ��
		delete[] bound;

		// �ͷŻ�������
		delete[] fitArrayCache;
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
		const size_t len = b.size();
		for (size_t i = 0; i < len; i++)
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
	void GAGroup<R(Args...)>::setStopTol(long double t, unsigned int N)
	{
		group_state.converCount = N;
		group_state.setStopTolFlag = true;
		group_state.stopTol = t;
	}

	// ���û���������
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMutateProb(double p)
	{
		mutateProb = p;
	}

	// ���û��򽻲����
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setCrossProb(double p)
	{
		crossProb = p;
	}

	// ���ò��м�����߳�����Ĭ��Ϊ1
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setThreadNum(const int NUM)
	{
		if (NUM >= 1)
		{
			thread_sync.setThreadNum(NUM);
		}
		else
		{
			throw std::string("Thread number error.");
		}
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

	// ��ȡ��Ⱥ���壨ͨ���±꣩
	template<class R, class... Args>
	const Individual& GAGroup<R(Args...)>::getIndivByIndex(const int index)const
	{
		return *(indivs + index);
	}
	
	/************************************************************************/
	/************************************************************************/
	/************************************************************************/
	// ��ȡ���Ÿ���
	template<class R, class... Args>
	void GAGroup<R(Args...)>::test()
	{
		using namespace std;

		// �Ż�����Ļ�����fitness
		std::vector<double> vec;
		for (int i = 0; i < nVars; i++)
		{
			vec.push_back(indivs[groupSize].vars[i]);
		}
		vec.push_back(indivs[groupSize].fitness);

		// ����Ż����
		std::cout << "Thread number: " << thread_sync.threadNum << std::endl;
		
		// ���ֹͣ����
		// Stop Code : -1-δֹͣ; 0-���Ž��������ȶ�ֵ; 1-�ﵽ����������; 2-�ﵽ������ʱ��; 3-��Ϊֹͣ����
		switch (group_state.stopCode)
		{
		case 0:
			std::cout << "Stop condition: reach the convergency." << std::endl;
			break;
		case 1:
			std::cout << "Stop condition: reach the max generation." << std::endl;
			break;
		case 2:
			std::cout << "Stop condition: reach the max time." << std::endl;
			break;
		}

		cout << vec[0] << endl;
		cout << vec[1] << endl;
		cout << "fitness: " << vec[2] << endl;

		// ����Ӵ����Ž��������
		cout << "\n*********************************" << endl;
		cout << "�Ӵ����Ž�������̣�" << endl;
		cout << bestIndivs[0].fitness << endl;
		cout << endl;
		for (size_t i = 1; i < bestIndivs.size(); i++)
		{
			if (bestIndivs[i].fitness != bestIndivs[i - 1].fitness)
			{
				cout << bestIndivs[i].fitness << endl;
				cout << endl;
			}
		}
	}

	/******************************************************************************************************************/

	// ��ʼ����
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::start()
	{
		// ��ʼ����Ⱥ
		initGroup();

		// ����Ⱥ�����������
		if (group_state.runable())
		{
			group_state.startTime = std::chrono::steady_clock::now();
			// ���߳�ģʽ
			if (thread_sync.threadNum > 1)
			{
				// ������Ⱥ�����߳���
				for (int i = 0; i < thread_sync.threadNum; i++)
				{
					vec_cross.emplace_back(&GAGroup<R(Args...)>::cross_thread, this, i);
				}

				// ������Ⱥ�����߳���
				for (int i = 0; i < thread_sync.threadNum; i++)
				{
					vec_mut.emplace_back(&GAGroup<R(Args...)>::mut_thread, this, i);
				}

				// ���컷��ѡ���߳���
				for (int i = 0; i < thread_sync.threadNum; i++)
				{
					vec_select.emplace_back(&GAGroup<R(Args...)>::sel_thread, this, i);
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

	// ������ǰ�߳�, �ȴ�������
	template<class R, class... Args>
	void GAGroup<R(Args...)>::wait_result()
	{
		// ���߳�ģʽ
		if (thread_sync.threadNum == 1)
		{
			single_thread.join();
		}

		// ���߳�ģʽ
		if (thread_sync.threadNum > 1)
		{
			for (int i = 0; i < thread_sync.threadNum; i++)
			{
				vec_cross[i].join();
				vec_mut[i].join();
				vec_select[i].join();
			}
		}
	}

	/******************************************* Private Functions *****************************************************/
	// ��ʼ����Ⱥ     
	template<class R, class... Args>
	void GAGroup<R(Args...)>::initGroup()
	{
		// ��ʼ��ǰ��Ҫ��֤�����û����������
		if (group_state.setBoundFlag)
		{
			// 1.��ʼ��ÿһ������
			for (int i = 0; i < groupSize; i++)
			{
				// ���ø����ʼ����
				for (int j = 0; j < nVars; j++)
				{
					indivs[i].vars[j] = random_real(bound[j][0], bound[j][1]);
				}
				// ������Ӧ��
				indivs[i].fitness = callFitFunc(indivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
			}

			// 2.Ѱ�����������Ӧ�ȸ���
			std::tie(group_state.worstIndex, group_state.bestIndex) = this->selectPolarIndivs(0, 1);
			indivs[groupSize] = indivs[group_state.bestIndex];                 // indivs����ĩβλ�����ڴ洢���Ž�
			bestIndivs.push_back(indivs[groupSize]);                     // ��¼���Ž�

			// 3.��ʼ��fitArrayCache����
			updateFitArrayCache();

			//4.���ó�ʼ����־λ
			group_state.initFlag = true;

			thread_sync.crossReady = true;
			thread_sync.mutReady = false;
			thread_sync.selectReady = false;
		}
		else
		{
			throw std::string("Fitness function or genes boundary is not set.");
		}
	}

	// ����(���߳�ģʽ)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::crossover()
	{
		int Index_M = 0;             // ������
		int Index_F = 0;             // ĸ����

		double rand_cross = 0;       // ���������    

		// ���ѡȡgroupSize/2�Ը��彻��, fitnessԽ�󣬱�ѡ�еĸ���Խ��
		for (int i = 0; i < groupSize; i += 2)
		{
			// ���ѡȡ����������Ϊ����
			Index_M = randomPickIndiv();
			Index_F = randomPickIndiv();

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

		// �������Ÿ����Ŵ����Ӵ�
		tempIndivs[groupSize] = indivs[groupSize];

		// �����Ӵ����建�����͸�������������ָ��
		Individual* temp = indivs;
		indivs = tempIndivs;
		tempIndivs = temp;
	}

	// ����(���߳�ģʽ)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::mutate()
	{
		double rand_num = 0; // �����

		// ÿ������ı���
		for (int i = 0; i < groupSize; i++)
		{
			// ÿ������ı���
			for (int j = 0; j < nVars; j++)
			{
				rand_num = random_real(0, 1);
				if (rand_num < mutateProb)
				{
					// �������
					indivs[i].vars[j] = mutate_PM(indivs[i].vars[j], bound[j][0], bound[j][1]);
				}
			}
			
			// ����ÿ��������Ӧ��
			indivs[i].fitness = callFitFunc(indivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
		}
	}

	// ѡ��(���߳�ģʽ)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::select()
	{
		// Ѱ���Ӵ��������Ÿ���
		std::tie(group_state.worstIndex, group_state.bestIndex) = this->selectPolarIndivs(0, 1);
		
		// ��̭�Ӵ�������, �ø������Ÿ���ȡ��
		indivs[group_state.worstIndex] = indivs[groupSize];
		// ����Ӵ����ָ��Ÿ���
		if (indivs[group_state.bestIndex].fitness >= indivs[groupSize].fitness)
		{
			indivs[groupSize] = indivs[group_state.bestIndex];
		}

		// �ٴ�Ѱ���Ӵ��������Ÿ���
		std::tie(group_state.worstIndex, group_state.bestIndex) = this->selectPolarIndivs(0, 1);

		// ��¼���Ž�
		bestIndivs.push_back(indivs[groupSize]);
	}

	// ��������(���߳�ģʽ)
	template<class R, class... Args>
	void GAGroup<R(Args...)>::run()
	{
		while (true)
		{
			// ���¸������̶Ŀ̶���
			updateFitArrayCache();

			crossover();   // ����
			mutate();      // ����
			select();      // ����ѡ��

			// ����ֹͣ״̬
			updateStopState();

			// �ж��Ƿ񵽴�ֹͣ����
			if (getStopFlag())
			{
				return;
			}
		}
	}

	// �������̶Ŀ̶���
	template<class R, class... Args>
	void GAGroup<R(Args...)>::updateFitArrayCache()
	{
        // fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		for (int i = 1; i < groupSize; i++)
		{
			fitArrayCache[i] = fitArrayCache[i - 1] + (indivs[i - 1].fitness - indivs[group_state.worstIndex].fitness);
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

	// Ѱ��������õĸ���λ��, �β�Ϊ�߳�ID(0 ~ N)�������߳�����, ����pair<worst, best>
	template<class R, class... Args>
	std::pair<int, int> GAGroup<R(Args...)>::selectPolarIndivs(const int seq, const int interval)
	{
		int worst = seq;
		int best = seq;

		// ����indivs��fitness��Ѱ��fitness�ļ�ֵ����¼��ֵλ��
		for (int i = seq + interval; i < groupSize; i += interval)
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

	// ���ݸ����fitness���ѡȡһ�����壬���ظ���λ��
	template<class R, class... Args>
	int GAGroup<R(Args...)>::randomPickIndiv()
	{
		// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		double temp = random_real(0, fitArrayCache[groupSize]);

		int lowIndex = 0;
		int upIndex = groupSize;
		int tempIndex;

		// ���ַ������������λ��,ʱ�临�Ӷ�log_n
		// ͨ�����ƶ�lowIndex��upIndex������±�λ�ã�ÿ�ν����������Сһ��
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
		return lowIndex;  // up = 1, low = 0 ʱ���������indivs[0]
	}

	template<class R, class... Args>
	template<std::size_t... I>
	inline R GAGroup<R(Args...)>::callFitFunc(double* args, const std::index_sequence<I...>&)
	{
		return fitFunc(args[I]...);
	}

	// �ж��Ƿ��������
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::getStopFlag()
	{
		// Stop Code : -1-δֹͣ; 0-���Ž��������ȶ�ֵ; 1-�ﵽ����������; 2-�ﵽ������ʱ��; 3-��Ϊֹͣ����

		// �Ƿ�����5�����Ž�Ĳ���С��ֹͣ���
		if (group_state.setStopTolFlag && group_state.count == group_state.converCount)
		{
			group_state.stopCode = 0;
		}

		// �Ƿ�ﵽ����������
		if (group_state.setMaxGeneFlag && group_state.nGene >= group_state.maxGene)
		{
			group_state.stopCode = 1;
		}

		// �Ƿ�ﵽ������ʱ��
		std::chrono::duration<double> evolTime = group_state.nowTime - group_state.startTime;
		if (group_state.setRuntimeFlag && evolTime.count() >= group_state.maxRuntime.value)
		{
			group_state.stopCode = 2;
		}

		// ����stop code�ж���Ⱥ�Ƿ�ֹͣ����
		if (group_state.stopCode != -1)
		{
			group_state.stopFlag = true;
		}

		return group_state.stopFlag;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	template<class R, class... Args>
	void GAGroup<R(Args...)>::cross_thread(const int seq)
	{
		int Index_M = 0;
		int Index_F = 0;
		double rand_cross = 0;

		while (true)
		{
			// �ٽ���
			{
				std::unique_lock<std::mutex> lck(thread_sync.mtx);
				thread_sync.cv.wait(lck, [this, &seq]() {return !(this->thread_sync).cross_flag[seq] && (this->thread_sync).crossReady || this->getStopFlag(); });
			}// �뿪������, gm.thread_state.mtx��unlock()�����Զ�ִ��
			
			 // �Ƿ�ֹͣ����
			if (getStopFlag())
			{
				return;
			}

			// �������彻������Ӵ�
			for (int i = seq; i < groupSize; i += 2 * thread_sync.threadNum)
			{
				// ���ѡȡ����������Ϊ����
				Index_M = randomPickIndiv();
				Index_F = randomPickIndiv();

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
					else    // ������
					{
						// �̳и������򣬻��򲻷�������
						tempIndivs[i].vars[j] = indivs[Index_M].vars[j];
						tempIndivs[i + 1].vars[j] = indivs[Index_F].vars[j];
					}
				}
			}

			// �߳�ͬ��
			thread_sync.mtx.lock();
			thread_sync.cross_sync(seq); // Change this function
			thread_sync.mtx.unlock();

			thread_sync.cv.notify_all();
		}
	}

	////////////////////////////////////////////////////////
	template<class R, class... Args>
	void GAGroup<R(Args...)>::mut_thread(const int seq)
	{
		double rand_num = 0; // �����

		while (true)
		{
			// �ٽ���
			{
				std::unique_lock<std::mutex> lck(thread_sync.mtx);
				thread_sync.cv.wait(lck, [this, &seq]() {return !(this->thread_sync).mut_flag[seq] && (this->thread_sync).mutReady || this->getStopFlag(); });
			}

			// �Ƿ�ֹͣ����
			if (getStopFlag())
			{
				return;
			}

			// ÿ���������(���Ž���岻��������)
			for (int i = seq; i < groupSize; i += thread_sync.threadNum)
			{
				if (i != group_state.bestIndex)
				{
					// ÿ������ı���
					for (int j = 0; j < nVars; j++)
					{
						rand_num = random_real(0, 1);
						if (rand_num < mutateProb)
						{
							// �������
							indivs[i].vars[j] = mutate_PM(indivs[i].vars[j], bound[j][0], bound[j][1]);
						}
					}
					// ���������Ӧ��
					indivs[i].fitness = callFitFunc(indivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
				}
			}

			// �߳�ͬ��
			thread_sync.mtx.lock();
			thread_sync.mutate_sync(seq);
			thread_sync.mtx.unlock();

			thread_sync.cv.notify_all();
		}
	}

	// ѡ��
	template<class R, class ...Args>
	void GAGroup<R(Args...)>::sel_thread(const int seq)
	{
		int worst_temp = 0;
		int best_temp = 0;

		while (true)
		{
			// �ٽ���
			{
				std::unique_lock<std::mutex> lck(thread_sync.mtx);				
				thread_sync.cv.wait(lck, [this, &seq]() {return !(this->thread_sync).sel_flag[seq] && (this->thread_sync).selectReady || this->getStopFlag(); });
			}

			// Ѱ���Ӵ��������Ÿ���
			std::tie(worst_temp, best_temp) = this->selectPolarIndivs(seq, thread_sync.threadNum);

			// �߳�ͬ��
			thread_sync.mtx.lock();			
			//����������������λ��;
			if (indivs[worst_temp].fitness < indivs[group_state.worstIndex].fitness)
			{
				group_state.worstIndex = worst_temp;
			}
			if (indivs[best_temp].fitness > indivs[group_state.bestIndex].fitness)
			{
				group_state.bestIndex = best_temp;
			}

			thread_sync.select_sync(seq);
			thread_sync.mtx.unlock();

			thread_sync.cv.notify_all();

			// �Ƿ�ֹͣ����
			if (getStopFlag())
			{
				return;
			}
		}
	}
}

#endif