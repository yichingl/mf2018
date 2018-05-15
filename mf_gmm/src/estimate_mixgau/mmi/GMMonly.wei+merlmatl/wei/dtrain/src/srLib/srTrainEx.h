///////////////////////////////////////////////////////////////
//                                                           //
//  Copyright (c) 2008                              //
//  All rights reserved.                                     //
//                                                           //
//  srTrainEx.h: interface for the CGMM class.                     //
//                                                           //
///////////////////////////////////////////////////////////////


#ifndef _SR_TRAIN_EX
#define _SR_TRAIN_EX


#define LOGZERO	(-1.0E20)   		/* ~log(0) */
#define LSMALL		(-0.5E20)   		/* log values < LSMALL are set to LOGZERO */
#define VAR_FLOOR  (0.01)			// Variance floor, make sure the variance to be large enough
#define MIX_FLOOR	(0.05)			// Mixture floor, make sure the mixture has enough samples
#define MAX_CUS_TAG		11		// Maximum length of Tags in Custom GMM definition file
#define MAX_HTK_TAG		11		// Maximum length of Tags in HTK format GMM definition file
#define TAG_CUS_NUM		9		// Maximum types of Tags in GMM definition file
#define TAG_HTK_NUM		16		// Maxinum types of Tags in HTK format GMM definition file
#define MAX_FILES	20000		// Maximum parameter files to be loaded when training

#define VAR_DEFAULT 2			// Variance maximum

#ifdef PI
#undef PI                /* PI is defined in Linux */
#endif
#define PI   3.14159265358979



//-----------------------------------------------------------------------
// Member Variables

typedef struct _CGMM 
{
	bool bMLoaded;		// whether the model is loaded
	bool bPLoaded;		// whether the param are loaded
	
	int iFrm;				// number of total frames
	float **fBuf;
	float threshold;			//threshold for silence

	char cName[128];			// Model name
	char cParm[128];			// Param Type
	bool bDiag;				// if the covariance matrix diagonal
	int iDim;				// Vector size
	int iMix;			// Mixtures of the GMM
	
	double dmini,dmaxi;		   	 // if iVectorSize is 1,and statistical the range
	double *dLdet;				// determinant of the covariance matrix [mixture]
	double *dWt;				// weight of each mixture[mixture]
	double *dLwt;				// log of the weight coeficient
	double *dVfl;				// variance floor
	double *dVec;				// a temp vector;
	double **dtMu;				// a temp matrix;
	double **dtDCov;			// a temp matrix
	double ***dtFCov;			// a temp matrix
	double **dMu;				// mean vector [mixture,dimension]
	double **dDCov;			// covariance (diagonal) [mixture,dimension]
	double **dDInv;				// inverse of convariance
	double ***dFCov;			// covariance (diagonal) [mixture,dimension]
	double ***dFInv;				// inverse of convariance
	
	double **dTrans;
	
} CGMM;


void SayBye(char *msg);
void LookOut(char *msg);

int LoadGMMSet(char *fname);

int Initialize(CGMM **pCGMM);
int Release(CGMM *pCGMM);

int WriteGMMSet(char *fname, int imod);
int CleanUpMdl();
int CleanUpPar();
int ReesMdl(int MaxLoop); //10
int InitMdl(int MaxLoop); //20

int LoadParaList(FILE *fp, int ifilenum);

int GetNextTag(FILE*fp,char*tag_buffer);
double LogAdd(double lvar1,double lvar2);

int ChkVarFloor(bool bUpdateWVD);
int UpdateWVD();
int CalVarFloor();
double LMixProb(float *buffer,int MixIndex);

//int VQCheckZero(int*Counter);
int VQNearest(float*Vector,double*distance);
int VQInit(bool bMaxMin);
double VQDistortion(float*vec1,double*vec2);

#endif // _SR_TRAIN_EX