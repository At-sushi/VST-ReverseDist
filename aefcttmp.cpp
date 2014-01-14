#include <stdio.h>
#include <string.h>
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <intrin.h>
#include <xmmintrin.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>

extern "C" {
#include "fft_mmx.h"
}
#include ".\AEfctTmp.h"

ARDistProgram::ARDistProgram()
{
	fGain = fDistortion = fDenoise = fThreshold = 0.5f;

	name[0] = '\0';		// 暫定
}

//-------------------------------------------------------------------------------------------------------
// SSEに対応しているか確認
static bool IsSSESupported()
{
	int cpuInfo[4];

	__cpuid(cpuInfo, 1);
	if (!(cpuInfo[3] & 0x02000000))
	{
		MessageBox(NULL, "Error: SSE isn't supported on your PC.", "Plugin error", MB_ICONERROR | MB_OK);
		return false;
	}

	return true;
}

AudioEffect* createEffectInstance (audioMasterCallback audioMaster)
{
	return (IsSSESupported() ? new AReverseDist (audioMaster) : NULL);
}


AReverseDist::AReverseDist(audioMasterCallback audioMaster)
: AudioEffectX(audioMaster, kNumPrograms, kNumParams)
{
	ffft_init();
	setProgram(0);

	setNumInputs(2);
	setNumOutputs(2);

	setUniqueID('RdsT');
}

AReverseDist::~AReverseDist(void)
{
}

//------------------------------------------------------------------------
void AReverseDist::resume ()
{
	for (int i = 0; i < 2; i++)
		PhaseMaker[i].Reset();

	AudioEffectX::resume ();
}

//------------------------------------------------------------------------
void AReverseDist::setProgram (VstInt32 program)
{
	ARDistProgram* ap = &programs[program];

	curProgram = program;
	setParameter (kGain, ap->fGain);	
	setParameter (kDistortion, ap->fDistortion);
	setParameter (kDenoise, ap->fDenoise);
	setParameter (kThreshold, ap->fThreshold);
}

//------------------------------------------------------------------------
void AReverseDist::setProgramName (char *name)
{
	strcpy (programs[curProgram].name, name);
}

//------------------------------------------------------------------------
void AReverseDist::getProgramName (char *name)
{
	if (!strcmp (programs[curProgram].name, "Init"))
		sprintf (name, "%s %d", programs[curProgram].name, curProgram + 1);
	else
		strcpy (name, programs[curProgram].name);
}

//-----------------------------------------------------------------------------------------
bool AReverseDist::getProgramNameIndexed (VstInt32 category, VstInt32 index, char* text)
{
	if (index < kNumPrograms)
	{
		strcpy (text, programs[index].name);
		return true;
	}
	return false;
}

//------------------------------------------------------------------------
void AReverseDist::setParameter (VstInt32 index, float value)
{
	ARDistProgram* ap = &programs[curProgram];

	switch (index)
	{
	case kGain :		fGain = ap->fGain = value;				break;
	case kDistortion :	fDistortion = ap->fDistortion = value;	break;
	case kDenoise :		fDenoise = ap->fDenoise = value;		break;
	case kThreshold :	fThreshold = ap->fThreshold = value;	break;
	}
}

//------------------------------------------------------------------------
float AReverseDist::getParameter (VstInt32 index)
{
	switch (index)
	{
	case kGain :		return fGain;
	case kDistortion :	return fDistortion;
	case kDenoise :		return fDenoise;
	case kThreshold :	return fThreshold;
	}
	
	return 0;
}

//------------------------------------------------------------------------
void AReverseDist::getParameterName (VstInt32 index, char *label)
{
	switch (index)
	{
	case kGain :		strcpy (label, "Gain");			break;
	case kDistortion :	strcpy (label, "Distortion");	break;
	case kDenoise :		strcpy (label, "Denoise");		break;
	case kThreshold :	strcpy (label, "Threshold");	break;
	}
}

//------------------------------------------------------------------------
void AReverseDist::getParameterDisplay (VstInt32 index, char *text)
{
	switch (index)
	{
	case kGain :		dB2string(fGain, text, kVstMaxParamStrLen);					break;
	case kDistortion :	float2string (fDistortion, text, kVstMaxParamStrLen);		break;
	case kDenoise :		dB2string (fDenoise * 0.1f, text, kVstMaxParamStrLen);		break;
	case kThreshold :	dB2string (fThreshold * 0.005f, text, kVstMaxParamStrLen);	break;
	}
	;
}

//------------------------------------------------------------------------
void AReverseDist::getParameterLabel (VstInt32 index, char *label)
{
	switch (index)
	{
	case kGain :	strcpy (label, "dB");	break;
	case kDenoise :	strcpy (label, "dB");	break;
	case kThreshold :strcpy (label, "dB");	break;
	}
	;
}

//------------------------------------------------------------------------
bool AReverseDist::getEffectName (char* name)
{
	strcpy (name, "Bit-Invert Distortion");
	return false;
}

//------------------------------------------------------------------------
bool AReverseDist::getProductString (char* text)
{
	strcpy (text, "Bit-Invert Distortion");
	return false;
}

//------------------------------------------------------------------------
bool AReverseDist::getVendorString (char* text)
{
	strcpy (text, "Reinherd Corp.");
	return false;
}

//------------------------------------------------------------------------
void AReverseDist::setSampleRate (float sampleRate)
{
	for (int i = 0; i < 2; i++)
		PhaseMaker[i].OnSampleRateChanged(sampleRate);
	AudioEffectX::setSampleRate (sampleRate);
}

//------------------------------------------------------------------------
void AReverseDist::setBlockSize (VstInt32 blockSize)
{
	// hann窓初期化
	int blockSize2 = blockSize * 2;
	int divide = (int)ceil((double)blockSize2 / (double)FFTSIZE);

	while (blockSize2 % divide != 0) divide++;
	int windowSize = blockSize2 / divide;

	assert (windowSize <= FFTSIZE);
	hannWindow.resize(windowSize);
	for (int i = 0; i < windowSize; i++)
		hannWindow[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / windowSize);
	halfWindow[0].clear();
	halfWindow[0].resize(windowSize / 2, 0.0f);
	halfWindow[1].clear();
	halfWindow[1].resize(windowSize / 2, 0.0f);
	endWindow[0].clear();
	endWindow[0].resize(windowSize / 2, 0.0f);
	endWindow[1].clear();
	endWindow[1].resize(windowSize / 2, 0.0f);
	setInitialDelay(windowSize / 2);
	AudioEffectX::setBlockSize (windowSize);
}

//---------------------------------------------------------------------------
inline static float reduct(float source, float strength)
{
	float reducted = source - strength;

	return (source > 0.0f) ? max(reducted, 0.0f) : min(reducted, 0.0f);
}

inline static __m128 reductSSE(__m128 source, __m128 str128)
{
	const __m128 zero = _mm_setzero_ps();
	__m128 isplus = _mm_cmpnle_ps(source, zero);

	source = _mm_sub_ps(source, str128);

	__m128 result = _mm_and_ps(
						_mm_max_ps( source, zero ),
						isplus
					);
	result = _mm_or_ps(result, _mm_andnot_ps(
						isplus,
						_mm_min_ps( source, zero )
					));

	return result;	// = (source > 0.0f) ? max(reducted, 0.0f) : min(reducted, 0.0f)
}

inline static __m128 rsqrt24(__m128 source)
{
	__m128 half = _mm_set1_ps(0.5), third = _mm_set1_ps(3.0),
		rsqrt = _mm_rsqrt_ps(source);

	return _mm_mul_ps(
				_mm_mul_ps(half, rsqrt),
				_mm_sub_ps(
						third,
						_mm_mul_ps(
								_mm_mul_ps(source, rsqrt),
								rsqrt
						)
				)
			);
}

static void denoise(float strength, vcpx_vector& vpx1, vcpx_vector& vpx2)
{
	if (strength <= 0.0f)		// 変換不要
		return;

	for (int i = 0; i < FFTSIZE; i++)
	{
		float *p1re = &vpx1[i/4].re[i%4],
			*p1im = &vpx1[i/4].im[i%4],
			*p2re = &vpx2[i/4].re[i%4],
			*p2im = &vpx2[i/4].im[i%4];

		// 複素平面上の長さに合わせて減算
		float length1 = sqrt(*p1re * *p1re + *p1im * *p1im);
		*p1re = reduct(*p1re, strength * *p1re / length1);
		*p1im = reduct(*p1im, strength * *p1im / length1);

		float length2 = sqrt(*p2re * *p2re + *p2im * *p2im);
		*p2re = reduct(*p2re, strength * *p2re / length2);
		*p2im = reduct(*p2im, strength * *p2im / length2);
	}
}

// せっかくSSEが使えるんだから…
static void denoiseSSE(float strength, vcpx_t* vpx1, vcpx_t* vpx2)
{
	if (strength <= 0.0f)		// 変換不要
		return;

	__m128 str128 = _mm_set1_ps(strength);

	for (int i = 0; i < FFTSIZE / 4; i++)
	{
		__m128 p1re = _mm_load_ps(vpx1->re),
			p1im = _mm_load_ps(vpx1->im);

		// 複素平面上の長さに合わせて減算
		__m128 rlength1 = rsqrt24(
							_mm_add_ps(_mm_mul_ps(p1re, p1re), _mm_mul_ps(p1im, p1im))
						 );		// = 1 / sqrt(p1re * p1re + p1im * p1im)
		_mm_store_ps(vpx1->re, reductSSE(p1re,
			_mm_mul_ps(_mm_mul_ps(p1re, rlength1), str128)		// = str128 * (p1re / (1/rlength1))
			));
		_mm_store_ps((vpx1++)->im, reductSSE(p1im,
			_mm_mul_ps(_mm_mul_ps(p1im, rlength1), str128)		// = str128 * (p1re / (1/rlength1))
			));

		__m128 p2re = _mm_load_ps(vpx2->re),
			p2im = _mm_load_ps(vpx2->im);
		__m128 rlength2 = rsqrt24(
							_mm_add_ps(_mm_mul_ps(p2re, p2re), _mm_mul_ps(p2im, p2im))
						 );		// = 1 / sqrt(p2re * p2re + p2im * p2im)
		_mm_store_ps(vpx2->re, reductSSE(p2re,
			_mm_mul_ps(_mm_mul_ps(p2re, rlength2), str128)		// = str128 * (p2re / (1/rlength2))
			));
		_mm_store_ps((vpx2++)->im, reductSSE(p2im,
			_mm_mul_ps(_mm_mul_ps(p2im, rlength2), str128)		// = str128 * (p2re / (1/rlength2))
			));
	}
}

void AReverseDist::processDistortion(float* in1, float* ptmp1, int channel)
{
	assert(channel == 0 || channel == 1);

	if (fabs(*in1) > fThreshold * 0.005)
	{
		if (*in1 >= 0.0f)
			PhaseMaker[channel].MoveTo(-M_PI_2 / 4);
		else
			PhaseMaker[channel].MoveTo(M_PI_2 / 4);
	}
	else
		PhaseMaker[channel].MoveTo(0);

	*ptmp1 = *in1 + PhaseMaker[channel].Process() * fDistortion;
}

//---------------------------------------------------------------------------
void AReverseDist::processReplacing (float** inputs, float** outputs, VstInt32 sampleFrames)
{
	// FFT_SSE用
	static _MM_ALIGN16 float tmp1[FFTSIZE],			// static宣言はスタックオーバーフロー対策
		  tmp2[FFTSIZE];

	float *in1 = inputs[0],
		  *in2 = inputs[1],
		  *out1 = outputs[0],
		  *out2 = outputs[1],
		  *ptmp1,
		  *ptmp2;

	while (sampleFrames > 0)
	{
		// Distortion処理
		int sf2 = sampleFrames;
		int hannnum = (int)hannWindow.size() / 2;
		assert(min(sf2, FFTSIZE / 2) >= hannnum);
		ptmp1 = tmp1;
		ptmp2 = tmp2;
		for (int i = 0; i < hannnum; i++)
		{
			*ptmp1++ = endWindow[0][i] * hannWindow[i];
			*ptmp2++ = endWindow[1][i] * hannWindow[i];
		}

		for (int i = 0; i < FFTSIZE / 2; i++)
		{
			if (--sf2 >= 0 && hannnum < (int)hannWindow.size())
			{
				processDistortion(in1++, ptmp1, 0);
				endWindow[0][i] = *ptmp1;
				*ptmp1++ *= hannWindow[hannnum];

				processDistortion(in2++, ptmp2, 1);
				endWindow[1][i] = *ptmp2;
				*ptmp2++ *= hannWindow[hannnum++];
			}
			else
			{
				// ぜろー
				*ptmp1++ = 0;
				*ptmp2++ = 0;
			}
		}
		if (fDenoise > 0.0f)
		{
			static _MM_ALIGN16 vcpx_vector vpx1, vpx2;

			// FFT変換
			ffft(tmp1, vpx1);
			ffft(tmp2, vpx2);
			// ノイズ削減処理
			denoiseSSE(fDenoise * 0.1f, vpx1, vpx2);
			// 逆変換
			frfft(vpx1, tmp1);
			frfft(vpx2, tmp2);
		}

		int fftLimit = hannnum / 2;
		std::vector<float>::iterator ihalf1 = halfWindow[0].begin(), ihalf2 = halfWindow[1].begin();
		ptmp1 = tmp1;
		ptmp2 = tmp2;

		// どっちかの制限に到達するまで転送
		while (--fftLimit >= 0 &&
				--sampleFrames >= 0)
		{
			*out1++ = ((*ptmp1++) + (*ihalf1++)) * fGain;
			*out2++ = ((*ptmp2++) + (*ihalf2++)) * fGain;
		}
		assert(ihalf1 == halfWindow[0].end() && ihalf2 == halfWindow[1].end());
		memcpy(&halfWindow[0][0], tmp1 + (int)halfWindow[0].size(), (int)(halfWindow[0].size() * sizeof(float)));
		memcpy(&halfWindow[1][0], tmp2 + (int)halfWindow[1].size(), (int)(halfWindow[1].size() * sizeof(float)));
	}
}
