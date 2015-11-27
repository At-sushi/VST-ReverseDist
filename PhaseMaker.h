#pragma once
#include <deque>
#include <vector>


class CPhaseMaker
{
public:
	typedef struct
	{
			double from, to;
			int count;
	} PHASING;

	CPhaseMaker(void);
	~CPhaseMaker(void);

	void MoveTo(double phase)
	{
		if (phase != lastValue)
		{
			PHASING value = {lastValue, phase, 0};
			lastValue = phase;
			MovingState.push_back(value);
		}
	}

	void OnSampleRateChanged(float sampleRate);
	double Process(void);

	void Reset()	{
		lastValue = 0;
		MovingState.clear();
	}

private:
	std::deque<PHASING> MovingState;
	std::vector<double> SinCache;
	double lastValue;
};

