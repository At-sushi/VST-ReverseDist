#define _USE_MATH_DEFINES
#include "PhaseMaker.h"


CPhaseMaker::CPhaseMaker(void)
{
	Reset();

	SinCache.resize(4);
	for (int i = 0; i < SinCache.size(); i++)
		SinCache[i] = (-std::cos(M_PI_4 * i) + 1) / 2;
}


CPhaseMaker::~CPhaseMaker(void)
{
}


void CPhaseMaker::OnSampleRateChanged(float sampleRate)
{
	Reset();

	SinCache.resize( std::round(4 * 44100.0 / sampleRate) );
	for (int i = 0; i < SinCache.size(); i++)
		SinCache[i] = (-std::cos(M_PI_4 * i) + 1) / 2;
}


double CPhaseMaker::Process(void)
{
	double result = MovingState.empty() ? lastValue : MovingState[0].from;

	for (int i = 0; i < MovingState.size(); i++)
	{
		result += SinCache[MovingState[i].count++] * (MovingState[i].to - MovingState[i].from);
		if (MovingState[i].count >= SinCache.size())
		{
			MovingState.pop_front();
			i--;
		}
	}

	return result;
}
