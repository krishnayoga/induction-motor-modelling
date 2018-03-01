/* Ida Bagus Krishna Yoga Utama 
1506716983 */

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME phase3
#include "simstruc.h" 3
#include <math.h> 

#define U(element) (*uPtrs[element])

static void mdlInitializeSizes(SimStruct *S){ 
if (!ssSetNumInputPorts(S, 1)) return; 
ssSetInputPortWidth(S, 0, 2); 
ssSetInputPortDirectFeedThrough(S, 0, 1); 
ssSetInputPortOverWritable(S, 0, 1); 
if (!ssSetNumOutputPorts(S, 1)) return; 
ssSetOutputPortWidth(S, 0, 3); 
ssSetNumSampleTimes(S, 1); 

ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE); } 

static void mdlInitializeSampleTimes(SimStruct *S) { 
ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME); 
ssSetOffsetTime(S, 0, 0.0); } 

static void mdlOutputs(SimStruct *S, int_T tid) { 
real_T *Y = ssGetOutputPortRealSignal(S,0); 
InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 

int_T i; 
real_T t = ssGetT(S); 

Y[0] = U(0)*cos(2*3.14*U(1)*t);
Y[1] = U(0)*cos((2*3.14*U(1)*t) + ((2*3.14)/3));
Y[2] = U(0)*cos((2*3.14*U(1)*t) - ((2*3.14)/3));

}
static void mdlTerminate(SimStruct *S) 
{ }  

#ifdef MATLAB_MEX_FILE 

#include "simulink.c" 
#else 
#include "cg_sfun.h" 
#endif 