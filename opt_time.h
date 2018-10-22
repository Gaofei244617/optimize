#ifndef _OPT_TIME_H_
#define _OPT_TIME_H_

namespace opt
{
	struct Second
	{
		long double value;
		explicit Second(long double N) :value(N) {}
		Second& operator= (const Second& other)
		{
			if (this != &other)
			{
				value = other.value;
			}
			return *this;
		}
	};

	struct Minute
	{
		long double value;
		explicit Minute(long double N) :value(N) {}
		operator Second()
		{
			return Second(value * 60);
		}
	};

	struct Hour
	{
		long double value;
		explicit Hour(long double N) :value(N) {}
		operator Second()
		{
			return Second(value * 3600);
		}
	};
}

#endif
