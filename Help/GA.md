# GA类：GAGroup

## 类定义

```cpp
namespace opt
{
    template<class F> class GAGroup; // 无具体实现
    template<class R, class... Args> class GAGroup<R(Args...)>;  // 使用适应度函数偏特化GAGroup
}
```

## 构造函数

```cpp
// 构造函数，f:适应度函数,size:个体数量
GAGroup(R(*f)(Args...), const int size = 1000);
// 拷贝构造
GAGroup(const GAGroup<R(Args...)>& other);  
// 移动构造
GAGroup(GAGroup<R(Args...)>&& other);
// 禁用赋值函数
GAGroup<R(Args...)>& operator=(const GAGroup<R(Args...)>& other) = delete;
GAGroup<R(Args...)>& operator=(GAGroup<R(Args...)>&& other) = delete;
```

## 析构函数

```cpp
~GAGroup();
```

## 成员函数

### (1) 设置种群参数

```cpp
void setName(const std::string& str);
```

设置种群名称。

```cpp
void setBoundary(double(*b)[2]);
```

设置变量区间。

 ```cpp
void setBoundary(const GenBound& b);
```

设置变量区间。

 ```cpp
void setMaxGeneration(const unsigned int N);
```

设置最大迭代次数。

```cpp
void setMaxRuntime(const Second& time);
```

设置最大运行时间（秒）。

```cpp
void setStopTol(long double t, unsigned int N = 5);
```

设置最优解停止误差。

```cpp
void setMutateProb(double p);
```

设置基因变异概率。

```cpp
void setCrossProb(double p);
```

设置交叉概率。

```cpp
void setThreadNum(const int NUM);
```

设置并行计算的线程数，default=1。

### (2) 获取种群参数

```cpp
const std::string getName()const;
```

获取种群名称。

```cpp
int getNVars()const;
```

获取种群变量个数。

```cpp
int getGeneration()const;
```

获得当前种群代数。

```cpp
int getGroupSize()const;
```

获得当前种群个体数量。

```cpp
std::vector<Individual> getBestIndivs();
```

获取历次迭代的最优解。

```cpp
int getStopCode();
```

获取Stop Code。

### (3) 用户函数

```cpp
bool start();
```

开启迭代进化过程。

```cpp
void wait_result();
```

阻塞当前线程，等待优化结果。

```cpp
void pause();
```

暂停迭代(为保证数据一致性，需在一次完整迭代后pause)。

```cpp
void proceed();
```

继续迭代。

```cpp
void kill();
```

结束迭代。

```cpp
GAGroup<R(Args...)> clone();
```

克隆当前种群