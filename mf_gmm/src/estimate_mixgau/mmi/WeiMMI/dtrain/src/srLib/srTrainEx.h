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

#define VAR_FLOOR  (0.01)			// Variance floor, make sure the variance to be large enough
#define MIX_FLOOR	(0.05)			// Mixture floor, make sure the mixture has enough samples
#define MAX_CUS_TAG		11		// Maximum length of Tags in Custom GMM definition file
#define MAX_HTK_TAG		11		// Maximum length of Tags in HTK format GMM definition file
#define TAG_CUS_NUM		9		// Maximum types of Tags in GMM definition file
#define TAG_HTK_NUM		16		// Maxinum types of Tags in HTK format GMM definition file
#define MAX_FILES	20000		// Maximum parameter files to be loaded when training

#define VAR_DEFAULT 2			// Variance maximum
#define FILE_BLOCK 500
#ifdef PI
#undef PI                /* PI is defined in Linux */
#endif
#define PI   3.14159265358979

#define PROTO_MLE	11
#define PROTO_MMI	12
#define PROTO_MCE	13
#define MODEL_MLE	21
#define MODEL_MMI	22
#define MODEL_MCE	23



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
	
	double dLPrior;				// log prior probability
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

//-------------------------------------------------------------------
// Initialize the model
int Initialize(CGMM ***ppp, int imodnum);
int Release(CGMM **pp, int imodnum);

//--------------------------------------------------------------
// Load GMM from the definition files
// the definination can be either proto or trained models
int LoadGMMSet(char *fname, CGMM ***pCGMM);

//-------------------------------------------------------------------
// Load Feature into p->fBuf
// ifilenum : utterances to read from the current fp 
int LoadParaList(FILE *fp, int ifilenum, CGMM *p);
//-------------------------------------------------------------------
// show GMM configuration
void ShowGMMCfg(CGMM *p, int ii);

//------------------------------------------------------------------
// Write GMM definition
int WriteGMMSet(char *fname, int imod, CGMM *p);
int CleanUpMdl(CGMM *p);
int CleanUpPar(CGMM *p);


//--------------------------------------------------------------------
// Init the model for EM by VQ
//   MaxLoop :  Max loop of VQ iterations
int InitMdl(int MaxLoop, CGMM *p);
int ReesMdl(int MaxLoop, CGMM *p); //10
int ReesMdlMMI(int MaxLoop, CGMM **p, int imod);

int GetNextTag(FILE*fp,char*tag_buffer);
double LogAdd(double lvar1,double lvar2);

// if the diag elememts is replaced by the var floor, update logwt, inv, logdet
int ChkVarFloor(bool bUpdateWVD, CGMM *p);
//--------------------------------------------------------------
// Recalculate the Determinant for each mixture after the covariance matrixes are changed.
// plus update the log weight and the inverse covarance
int UpdateWVD(CGMM *p);
int CalVarFloor(CGMM **p, int imod, int itrained);
//----------------------------------------------------------------
// Return the log likelihood of the given vector to the given mixture 
double LMixProb(float *fparam, CGMM *p, int iMix);

//int VQCheckZero(int*Counter);
int VQNearest(float*Vector,double*distance);

//-----------------------------------------------------------------
// Init the mean vectors of VQ procedure
//   bMaxMin: 0 choose the very first points to initialize each mixture
//	      1 use points among extreme points to initialize each mixture
int VQInit(bool bMaxMin, CGMM **pp, int imod);
//-----------------------------------------------------------------
// Init the model for EM by VQ
//   MaxLoop :  Max loop of VQ iterations
int InitMdl(int MaxLoop, CGMM *pp, int imod);

//---------------------------------------------------------------
// Get the nearest mean vector to the given vector according to
// Euclidean distance.
//    fvector : Vector to calculate
//        err : minimal distance
int VQNearest(float *fvector, double *err, CGMM *p);

double VQDistortion(float *vec1, float *vec2, CGMM *p);

#endif // _SR_TRAIN_EX