#include <stdio.h>
#include <string.h>
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <intrin.h>
#include <xmmintrin.h>
#include <math.h>

extern "C" {
#include "fft_mmx.h"
}
#include ".\AEfctTmp.h"

ARDistProgram::ARDistProgram()
{
	fGain = fDistortion = fDenoise = fThreshold = fFeedback = 0.5f;

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
	PreSample[0] = PreSample[1] = 0.f;

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
	setParameter (kFeedback, ap->fFeedback);
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
	case kFeedback :	fFeedback = ap->fFeedback = value;		break;
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
	case kFeedback :	return fFeedback;
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
	case kFeedback :	strcpy (label, "Feedback");		break;
	}
}

//------------------------------------------------------------------------
void AReverseDist::getParameterDisplay (VstInt32 index, char *text)
{
	switch (index)
	{
	case kGain :		dB2string(fGain, text, kVstMaxParamStrLen);					break;
	case kDistortion :	float2string (fDistortion, text, kVstMaxParamStrLen);		break;
	case kDenoise :		dB2string (fDenoise * 0.05f, text, kVstMaxParamStrLen);		break;
	case kThreshold :	dB2string (fThreshold * 0.005f, text, kVstMaxParamStrLen);	break;
	case kFeedback:		dB2string (fFeedback, text, kVstMaxParamStrLen);			break;
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
	case kFeedback :strcpy (label, "dB");	break;
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
	AudioEffectX::setSampleRate (sampleRate);
}

//------------------------------------------------------------------------
void AReverseDist::setBlockSize (VstInt32 blockSize)
{
	AudioEffectX::setBlockSize (blockSize);
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
		__m128 rlength1 = _mm_rsqrt_ps(
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
		__m128 rlength2 = _mm_rsqrt_ps(
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

//---------------------------------------------------------------------------
void AReverseDist::processReplacing (float** inputs, float** outputs, VstInt32 sampleFrames)
{
	// FFT_SSE用
	_MM_ALIGN16 float tmp1[FFTSIZE],
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
		ptmp1 = tmp1;
		ptmp2 = tmp2;

		for (int i = 0; i < FFTSIZE; i++)
		{
			if (--sf2 >= 0)
			{
				if (fabs(*in1) > fThreshold * 0.005)
					*ptmp1 = *in1 >= 0.0f ? *in1++ - fDistortion : *in1++ + fDistortion;
				else
					*ptmp1 = *in1++;
				// フィードバック
				*ptmp1 = PreSample[0] = *ptmp1 - PreSample[0] * fFeedback;
				ptmp1++;

				if (fabs(*in2) > fThreshold * 0.005)
					*ptmp2 = *in2 >= 0.0f ? *in2++ - fDistortion : *in2++ + fDistortion;
				else
					*ptmp2 = *in2++;
				// フィードバック
				*ptmp2 = PreSample[1] = *ptmp2 - PreSample[1] * fFeedback;
				ptmp2++;
			}
			else
			{
				// 前の値をコピー
				*ptmp1++ = tmp1[i - 1];
				*ptmp2++ = tmp2[i - 1];
			}
		}

		if (fDenoise > 0.0f)
		{
			_MM_ALIGN16 vcpx_vector vpx1, vpx2;

			// FFT変換
			ffft(tmp1, vpx1);
			ffft(tmp2, vpx2);
			// ノイズ削減処理
			denoiseSSE(fDenoise * 0.05f, vpx1, vpx2);
			// 逆変換
			frfft(vpx1, tmp1);
			frfft(vpx2, tmp2);
		}

		int fftLimit = FFTSIZE;
		ptmp1 = tmp1;
		ptmp2 = tmp2;

		// どっちかの制限に到達するまで転送
		while (--fftLimit >= 0 &&
				--sampleFrames >= 0)
		{
			*out1++ = *ptmp1++ * fGain;
			*out2++ = *ptmp2++ * fGain;
		}
	}
}
