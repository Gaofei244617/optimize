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
#include "GA_State.h"
#include "cross_factor.h"
#include "mutate_factor.h"
#include "GA_ThreadPool.h"

/*************Test*******************/
#include <iostream>
#include "out_bool.h"
/************************************/

namespace opt
{
	template<class F>
	class GAGroup;

	// ��Ⱥ�࣬���ṩ��Ӧ�Ⱥ�������
	template<class R, class... Args>
	class GAGroup<R(Args...)>
	{
		using GenBound = std::initializer_list<std::initializer_list<double>>;

	private:
		std::string name;                                                     // ��Ⱥ����
		int groupSize;                                                        // ��ʼ��Ⱥ����������default = 1000��
		const int nVars;                                                      // ��Ӧ�Ⱥ��������ı�������    
		Individual* indivs;                                                   // ��Ⱥ����, ָ��groupSize + 1����������(���һ��������Ÿ���)
		Individual* tempIndivs;                                               // �Ӵ����建����
		R(*fitFunc)(Args...);                                                 // ��Ӧ�Ⱥ���ָ��
		double(*bound)[2];                                                    // ÿ������(����)������, ������ָ���ʾ
		double* fitArrayCache;                                                // fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		double mutateProb;                                                    // ����������,Ĭ��p = 0.1
		double crossProb;                                                     // ���彻�����, Ĭ��p = 0.6
		
		GA_State state;                                                       // ��Ⱥ״̬
		GA_ThreadState thread_state;                                          // 

	public:
		std::vector<Individual> bestIndivs;                                   // ��¼ÿ�ε��������Ÿ���	

		GAGroup(R(*f)(Args...), const int size = 1000);                       // ���캯��������һ����Ⱥ�����ṩ��Ӧ�Ⱥ�������Ⱥ����
		GAGroup(const GAGroup<R(Args...)>& other);                            // ��������
		GAGroup(GAGroup<R(Args...)>&& other);                                 // �ƶ�����
		~GAGroup();

		void setName(const std::string& str);                                 // ������Ⱥ����
		void setBoundary(double(*b)[2]);                                      // ���ñ�������
		void setBoundary(const GenBound& b);                                  // ���ñ�������
		void setMaxGeneration(int N);                                         // ��������������
		void setMaxRuntime(const Second& time);                               // �����������ʱ�䣨�룩
		void setStopTol(long double t);                                       // �������Ž�ֹͣ���
		void setMutateProb(double p);                                         // ���û���������
		void setCrossProb(double p);                                          // ���ý������
		void setThreadNum(const int NUM);                                     // ���ò��м�����߳�����Ĭ��Ϊ1

		const std::string getName()const;                                          // ��ȡ��Ⱥ����
		int getNVars()const;                                                  // ��ȡ��Ⱥ�������� 
		int getGeneration()const;                                             // ��õ�ǰ��Ⱥ����
		int getGroupSize()const;                                              // ��õ�ǰ��Ⱥ��������	
		const Individual& getIndivByIndex(int index)const;                    // ��ȡ��Ⱥ���壨ͨ���±꣩

		void test();                                                       ///////////////////////////////

		bool start();                                                         // ��ʼ����
		bool single_start();
		bool pause();                                                         // ֹͣ����
		bool proceed();                                                       // ��������                                                    

	private:
		void initGroup();                                                     // ��ʼ����Ⱥ����
		bool crossover();                                                     // ����
		void cross_thread(const int seq);
		void mut_thread(const int seq);
		void select();
		bool mutate();                                                        // ����
		void run();                                                           // ��������
		std::pair<int, int> findPolarIndex();                                 // Ѱ��������õĸ���λ��,����<worst, best>
		int randomPickIndiv();                                                // ���ݸ����fitness���ѡȡһ�����壬���ظ���λ��
		template<std::size_t... I>
		R callFitFunc(double* args, const std::index_sequence<I...>&);        // ������Ӧ�Ⱥ���
		bool flushStopFlag();                                                 // �ж��Ƿ��������
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
		crossProb(0.6)
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
		state(other.state)
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
		state(other.state)
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
		state.setBoundFlag = true;
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
		state.setBoundFlag = true;
	}

	// ��������������
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMaxGeneration(int N)
	{
		state.setMaxGeneFlag = true;
		state.maxGene = N;
	}

	// �����������ʱ�䣨�룩
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setMaxRuntime(const Second& time)
	{
		state.setRuntimeFlag = true;
		state.maxRuntime = time;
	}

	// �������Ž�ֹͣ���
	template<class R, class... Args>
	void GAGroup<R(Args...)>::setStopTol(long double t)
	{
		state.setStopTolFlag = true;
		state.stopTol = t;
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
			thread_state.threadNum = NUM;
			thread_state.cross_flag.set_length(NUM);
			thread_state.mut_flag.set_length(NUM);
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
		return state.nGene;
	}

	// ��ȡ��ǰ��Ⱥ��������
	template<class R, class... Args>
	int GAGroup<R(Args...)>::getGroupSize()const
	{
		return groupSize;
	}

	// ��ȡ��Ⱥ���壨ͨ���±꣩
	template<class R, class... Args>
	const Individual& GAGroup<R(Args...)>::getIndivByIndex(int index)const
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
		std::cout << "Stop code: " << state.stopCode << std::endl;
		cout << vec[0] << endl;
		cout << vec[1] << endl;
		cout << "fitness: " << vec[2] << endl;
		cout << "**************************************" << endl;

		// ����Ӵ����Ž��������
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
		//std::cout << "aaaaa" << std::endl;

		std::vector<std::thread> vec_cross;
		std::vector<std::thread> vec_mut;

		if (state.runable())
		{
			//state.startTime = std::chrono::steady_clock::now();
			//// �������߳��������Ⱥ����
			//std::thread t(&GAGroup<R(Args...)>::run, this);
			//t.join();
			//return true;

			// ������Ⱥ�����̳߳�
			for (int i = 0; i < thread_state.threadNum; i++)
			{
				vec_cross.emplace_back(&GAGroup<R(Args...)>::cross_thread, this, i);
			}

			// ������Ⱥ�����̳߳�
			for (int i = 0; i < thread_state.threadNum; i++)
			{
				vec_mut.emplace_back(&GAGroup<R(Args...)>::mut_thread, this, i);
			}

			// ����ѡ��,  ���߳���ִ��
			select();

			// ������ǰ�߳�
			for (int i = 0; i < thread_state.threadNum; i++)
			{
				vec_cross[i].join();
				vec_mut[i].join();
			}

			return true;
		}
		return false;
	}

	// ���߳�ִ��
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::single_start()
	{
		initGroup();

		if (state.runable())
		{
			state.startTime = std::chrono::steady_clock::now();
			std::thread t(&GAGroup<R(Args...)>::run, this);
			t.join();
			return true;
		}
		return false;
	}

	// ��ͣ����
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::pause()
	{
		state.stopFlag = true;
		state.stopCode = 3;
		return true;
	}

	// ��������
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::proceed()
	{
		state.stopFlag = false;
		state.stopCode = -1;
		std::thread t(&GAGroup<R(Args...)>::run, this);
		t.join();
		return true;
	}

	/******************************************* Private Functions *****************************************************/
	// ��ʼ����Ⱥ
	template<class R, class... Args>
	void GAGroup<R(Args...)>::initGroup()
	{
		// ��ʼ��ǰ��Ҫ��֤�����û����������
		if (state.setBoundFlag)
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
			std::tie(state.worstIndex, state.bestIndex) = this->findPolarIndex();
			indivs[groupSize] = indivs[state.bestIndex];                 // indivs����ĩβλ�����ڴ洢���Ž�
			bestIndivs.push_back(indivs[groupSize]);                     // ��¼���Ž�

			// 3.��ʼ��fitArrayCache����
			// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
			for (int i = 1; i < groupSize; i++)
			{
				fitArrayCache[i] = fitArrayCache[i - 1] + (indivs[i - 1].fitness - indivs[state.worstIndex].fitness);
			}

			//4.���ó�ʼ����־λ
			state.initFlag = true;

			thread_state.crossReady = true;
			thread_state.mutReady = false;
			thread_state.selectReady = false;
		}
		// throw std::string("Fitness function or genes boundary is not set.");
	}

	template<class R, class... Args>
	bool GAGroup<R(Args...)>::crossover()
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
				else    // ������
				{
					// �̳и������򣬻��򲻷�������
					tempIndivs[i].vars[j] = indivs[Index_M].vars[j];
					tempIndivs[i + 1].vars[j] = indivs[Index_F].vars[j];
				}
			}

			// �����Ӵ�������Ӧ��
			tempIndivs[i].fitness = callFitFunc(tempIndivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
			tempIndivs[i + 1].fitness = callFitFunc(tempIndivs[i + 1].vars, std::make_index_sequence<sizeof...(Args)>());
		}

		// �������Ÿ����Ŵ����Ӵ�
		tempIndivs[groupSize] = indivs[groupSize];

		// �����Ӵ����建�����͸�����Ⱥ����������ָ��
		Individual* temp = indivs;
		indivs = tempIndivs;
		tempIndivs = temp;

		// Ѱ���Ӵ��������Ÿ���
		std::tie(state.worstIndex, state.bestIndex) = this->findPolarIndex();

		// ����Ӵ����ָ��Ÿ���
		if (indivs[state.bestIndex].fitness >= indivs[groupSize].fitness)
		{
			indivs[groupSize] = indivs[state.bestIndex];
		}

		return true;
	}

	// ����
	template<class R, class... Args>
	bool GAGroup<R(Args...)>::mutate()
	{
		double rand_num = 0; // �����

		// ÿ������ı���,���Ž���岻��������
		for (int i = 0; i < groupSize; i++)
		{
			if (i != state.bestIndex)
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
				// ������Ӧ��
				indivs[i].fitness = callFitFunc(indivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
			}
		}

		// ����fitArrayCache����
		// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
		for (int i = 1; i < groupSize; i++)
		{
			fitArrayCache[i] = fitArrayCache[i - 1] + (indivs[i - 1].fitness - indivs[state.worstIndex].fitness);
		}

		// Ѱ���Ӵ��������Ÿ���
		std::tie(state.worstIndex, state.bestIndex) = this->findPolarIndex();

		// ����������ָ��Ÿ���
		if (indivs[state.bestIndex].fitness >= indivs[groupSize].fitness)
		{
			indivs[groupSize] = indivs[state.bestIndex];
		}

		return true;
	}

	// ��������
	template<class R, class... Args>
	void GAGroup<R(Args...)>::run()
	{
		while (true)
		{
			crossover();
			mutate();

			// ��¼���Ž�
			bestIndivs.push_back(indivs[groupSize]);        

			// ����������һ
			state.nGene++;
			// ��¼��ǰʱ��
			state.nowTime = std::chrono::steady_clock::now();

			// �ж����Ž��������һ�εĲ���ֵ
			if (abs(bestIndivs[state.nGene - 1].fitness - bestIndivs[state.nGene].fitness) <= state.stopTol)
			{
				state.count++;
			}
			else
			{
				state.count = 0;
			}

			// �ж��Ƿ񵽴�ֹͣ����
			if (flushStopFlag())
			{
				return;
			}
		}
	}

	// Ѱ��������õĸ���λ��,����pair<worst, best>
	template<class R, class... Args>
	std::pair<int, int> GAGroup<R(Args...)>::findPolarIndex()
	{
		int worst = 0;
		int best = 0;

		// ����indivs��fitness��Ѱ��fitness�ļ�ֵ����¼��ֵλ��
		for (int i = 1; i < groupSize; i++)
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
	bool GAGroup<R(Args...)>::flushStopFlag()
	{
		// �Ƿ�ﵽ����������
		if (state.setMaxGeneFlag && state.nGene >= state.maxGene)
		{
			state.stopCode = 1;
		}

		// �Ƿ�ﵽ������ʱ��
		std::chrono::duration<double> evolTime = state.nowTime - state.startTime;
		if (state.setRuntimeFlag && evolTime.count() >= state.maxRuntime.value)
		{
			state.stopCode = 2;
		}

		// �Ƿ�����5�����Ž�Ĳ���С��ֹͣ���
		if (state.setStopTolFlag && state.count == 5)
		{
			state.stopCode = 0;
		}

		// ����stop code�ж���Ⱥ�Ƿ�ֹͣ����
		if (state.stopCode != -1)
		{
			state.stopFlag = true;
		}

		return state.stopFlag;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	template<class R, class... Args>
	void GAGroup<R(Args...)>::cross_thread(const int seq)
	{
		//std::condition_variable cv_cross;
		int Index_M = 0;
		int Index_F = 0;
		double rand_cross = 0;

		while (true)
		{
			{
				std::unique_lock<std::mutex> lck(thread_state.mtx);
				//std::cout << "---* Cross before: "<< seq << std::endl;
				//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
				//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
				//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
				thread_state.cv.wait(lck, [this, &seq]() {return !(this->thread_state).cross_flag[seq] && (this->thread_state).crossReady || this->flushStopFlag(); });
			}// �뿪������, thread_state.mtx��unlock()�����Զ�ִ��
			
			//std::this_thread::sleep_for(1s);

			 // �Ƿ�ֹͣ����
			if (flushStopFlag())
			{
				//cout << "Exit Cross..." << endl;
				return;
			}

			for (int i = seq; i < groupSize; i += 2 * thread_state.threadNum)
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

				// �����Ӵ�������Ӧ��
				tempIndivs[i].fitness = callFitFunc(tempIndivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
				tempIndivs[i + 1].fitness = callFitFunc(tempIndivs[i + 1].vars, std::make_index_sequence<sizeof...(Args)>());
			}

			// �߳�ͬ��
			thread_state.mtx.lock();
			thread_state.cross_flag[seq] = true;
			thread_state.selectReady = false;
			//std::cout << "---* Cross after: " << seq << "\n" << std::endl;
			if (thread_state.cross_flag.is_all_true())
			{
				thread_state.cross_flag.set_all(false);
				thread_state.crossReady = false;
				thread_state.mutReady = true;
				//std::cout << "--------Cross finish once in thread: --------" << seq << std::endl;
			}
			else
			{
				thread_state.mutReady = false;
			}
			//std::cout << "\n++ Change in cross: " << seq << endl;
			//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
			//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
			//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
			thread_state.mtx.unlock();
			thread_state.cv.notify_all();
		}
	}

	////////////////////////////////////////////////////////
	template<class R, class... Args>
	void GAGroup<R(Args...)>::mut_thread(const int seq)
	{
		//std::condition_variable cv_mut;  // ��������
		double rand_num = 0; // �����

		while (true)
		{
			{
				// �ٽ���
				std::unique_lock<std::mutex> lck(thread_state.mtx);
				//std::cout << "---# Mutate before: " << seq << std::endl;
				//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
				//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
				//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
				thread_state.cv.wait(lck, [this, &seq]() {return !(this->thread_state).mut_flag[seq] && (this->thread_state).mutReady || this->flushStopFlag(); });
			}

			//std::this_thread::sleep_for(1s);

			// �Ƿ�ֹͣ����
			if (flushStopFlag())
			{
				//cout << "Exit Mutate..." << endl;
				return;
			}

			// ÿ���������(���Ž���岻��������)
			for (int i = seq; i < groupSize; i += thread_state.threadNum)
			{
				if (i != state.bestIndex)
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
					// ������Ӧ��
					indivs[i].fitness = callFitFunc(indivs[i].vars, std::make_index_sequence<sizeof...(Args)>());
				}
			}

			// �߳�ͬ��
			thread_state.mtx.lock();
			thread_state.mut_flag[seq] = true;
			thread_state.crossReady = false;
			//std::cout << "---# Mutate after: " << seq << "\n" << std::endl;
			if (thread_state.mut_flag.is_all_true())
			{
				thread_state.mut_flag.set_all(false);
				thread_state.mutReady = false;
				thread_state.selectReady = true;
				//std::cout << "--------Mutate finish once--------" << std::endl;
			}
			else
			{
				thread_state.selectReady = false;
			}
			//std::cout << "\n++ Change in mutate: " << seq << endl;
			//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
			//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
			//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
			thread_state.mtx.unlock();
			thread_state.cv.notify_all();
		}
	}

	// ѡ��
	template<class R, class ...Args>
	void GAGroup<R(Args...)>::select()
	{
		//std::condition_variable cv_select;

		while (true)
		{
			// �ٽ���
			{
				std::unique_lock<std::mutex> lck(thread_state.mtx);				
				//std::cout << "---^ Select before. " << std::endl;
				//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
				//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
				//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
				thread_state.cv.wait(lck, [this]() {return (this->thread_state).selectReady; });
			}
			
			//std::this_thread::sleep_for(1s);
			// ����fitArrayCache����
			// fitArrayCache[n] - fitArrayCache[n-1] == Individual[n-1].fitness - minFitness
			// fitArrayCache[0]ʼ��Ϊ0
			for (int i = 1; i < groupSize; i++)
			{
				fitArrayCache[i] = fitArrayCache[i - 1] + (indivs[i - 1].fitness - indivs[state.worstIndex].fitness);
			}

			// Ѱ���Ӵ��������Ÿ���
			std::tie(state.worstIndex, state.bestIndex) = this->findPolarIndex();

			// ����������ָ��Ÿ���
			if (indivs[state.bestIndex].fitness >= indivs[groupSize].fitness)
			{
				indivs[groupSize] = indivs[state.bestIndex];
			}

			// �߳�ͬ��
			thread_state.mtx.lock();
			state.nGene++;
			thread_state.selectReady = false;
			thread_state.crossReady = true;
			thread_state.mutReady = false;
			//std::cout << "---^ Select after." << "\n" << std::endl;
			//std::cout << "--------Select finish once--------" << std::endl;
			//std::cout << "�� " << state.nGene << " ��..." << endl;
			//std::cout << "\n++ Change in select: " << endl;
			//std::cout << "CrossReady: " << out_bool(thread_state.crossReady) << endl;
			//std::cout << "MutateReady: " << out_bool(thread_state.mutReady) << endl;
			//std::cout << "SelecReady: " << out_bool(thread_state.selectReady) << endl;
			thread_state.mtx.unlock();
			thread_state.cv.notify_all();
			// �Ƿ�ֹͣ����
			if (flushStopFlag())
			{
				//cout << "Exit Select..." << endl;
				return;
			}
		}
	}
}

#endif