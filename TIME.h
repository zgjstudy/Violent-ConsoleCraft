#ifndef TIME_H
#define TIME_H

#include <windows.h>

namespace RENDERER
{
	//精确计时类
	class Timestamp
	{
	protected:
		LARGE_INTEGER _frequency;
		LARGE_INTEGER _startCount;

	public:
		Timestamp()
		{
			QueryPerformanceFrequency(&_frequency);
			QueryPerformanceCounter(&_startCount);
		}

		//获取微秒
		double getElapsedTimeInMicroSec()
		{
			LARGE_INTEGER endCount;
			QueryPerformanceCounter(&endCount);

			double  startTimeInMicroSec = _startCount.QuadPart * (1000000.0 / _frequency.QuadPart);
			double  endTimeInMicroSec = endCount.QuadPart * (1000000.0 / _frequency.QuadPart);

			return  endTimeInMicroSec - startTimeInMicroSec;
		}

		//获得毫秒
		double getElapsedTimeInMilliSec()
		{
			return this->getElapsedTimeInMicroSec() * 0.001;
		}

		//刷新时钟
		void update()
		{
			QueryPerformanceCounter(&_startCount);
		}
	};

	//高速rand函数
	inline int randEx()
	{
		LARGE_INTEGER seed;
		QueryPerformanceFrequency(&seed);
		QueryPerformanceCounter(&seed);
		srand(seed.QuadPart);

		return rand();
	}
}

#endif // !TIME_H