#pragma once
#include "public.sdk\source\vst2.x\audioeffectx.h"
#include "PhaseMaker.h"

enum
{
	// Global
	kNumPrograms = 16,

	// Parameters Tags
	kGain = 0,
	kDistortion,
	kDenoise,
	kThreshold,
	kFeedback,

	kNumParams
};

class ARDistProgram
{
	friend class AReverseDist;
public:
	ARDistProgram();
	~ARDistProgram() {}

private:
	float fGain;
	float fDistortion;
	float fDenoise;
	float fThreshold;
	float fFeedback;
	char name[24];
};

class AReverseDist :
	public AudioEffectX
{
public:
	AReverseDist(audioMasterCallback audioMaster);
	~AReverseDist(void);

	//---from AudioEffect-----------------------
	virtual void processReplacing (float** inputs, float** outputs, VstInt32 sampleFrames);

	virtual void setProgram (VstInt32 program);
	virtual void setProgramName (char* name);
	virtual void getProgramName (char* name);
	virtual bool getProgramNameIndexed (VstInt32 category, VstInt32 index, char* text);
	
	virtual void setParameter (VstInt32 index, float value);
	virtual float getParameter (VstInt32 index);
	virtual void getParameterLabel (VstInt32 index, char* label);
	virtual void getParameterDisplay (VstInt32 index, char* text);
	virtual void getParameterName (VstInt32 index, char* text);

	virtual void setSampleRate (float sampleRate);
	virtual void setBlockSize (VstInt32 blockSize);

	virtual void resume ();

	virtual bool getEffectName (char* name);
	virtual bool getVendorString (char* text);
	virtual bool getProductString (char* text);
	virtual VstInt32 getVendorVersion () { return 1000; }

protected:
	ARDistProgram programs[kNumPrograms];
	CPhaseMaker PhaseMaker[2];
	
	float fGain;
	float fDistortion;
	float fDenoise;
	float fThreshold;
	float fFeedback;
	float PreSample[2];
	double dFeedBack2;
};
