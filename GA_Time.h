#ifndef _GA_TIME_H_
#define _GA_TIME_H_

namespace opt
{
	struct Second
	{
		const long double value;
		explicit Second(long double N) :value(N) {}
	};

	struct Minute
	{
		const long double value;
		explicit Minute(long double N) :value(N) {}
		operator Second()
		{
			return Second(value * 60);
		}
	};

	struct Hour
	{
		const long double value;
		explicit Hour(long double N) :value(N) {}
		operator Second()
		{
			return Second(value * 3600);
		}
	};
}

#endif
