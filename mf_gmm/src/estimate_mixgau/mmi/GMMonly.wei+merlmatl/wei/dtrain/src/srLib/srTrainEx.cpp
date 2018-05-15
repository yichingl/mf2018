///////////////////////////////////////////////////////////////
//                                                           //
//                                                           //
//  GMM.cpp: implementation of the CGMM class.               //
//                                                           //
///////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>
#include <float.h>
#include <stddef.h>
#include <stdbool.h>

#include "srShellEx.h"
#include "srMathEx.h"
#include "srTrainEx.h"


static char msg[1000];
static CGMM *p;
static double minLogExp;
static double dsqrt2pi; // sqrt (2 * PAI)
static double dlog2pi; // log(2*pai)


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void SayBye(char *msg)
{
	//fprintf(stderr,"\n *** Error around LINE %d at FILE \"%s\" !\n\n     Error Message: %s\n\n",__LINE__,__FILE__,msg);
	printf("Error Message: %s\n\n", msg);
	exit(-1);
}
void LookOut(char *msg)
{
	printf("Warining  %s !\n", msg);
}

//--------------------------------------------------------------
// Free model's memory before exit
int CleanUpMdl()
{
	if (p->bMLoaded) {
		free( p->dWt  );
		free( p->dLwt );
		free( p->dVfl );
		free( p->dVec);
		
		free2d( (char **) p->dMu );
		free2d( (char **) p->dDCov );
		free2d( (char **) p->dDInv );
		free2d( (char **) p->dtMu );
		free2d( (char **) p->dtDCov );
		
		free3d( (char ***)p->dFCov );
		free3d( (char ***)p->dFInv );
		free3d( (char ***)p->dtFCov );
		
		p->bMLoaded = false;
	}
	return 1;
}

//-------------------------------------------------------------------
// Free the parameter's memoroy
int CleanUpPar()
{
	int i;
	if (p->bPLoaded){
		free2d((char **)p->fBuf);
	
		p->bPLoaded = false;
	}
	return 1;
}

//-------------------------------------------------------------------
// Initialize the model
int Initialize(CGMM **pCGMM)
{
	*pCGMM = (CGMM *) malloc(sizeof(CGMM));
	p = *pCGMM;
	p->bMLoaded		= false;
	p->bPLoaded		= false;
	p->bDiag			= true;
	p->fBuf			= NULL;
	p->iDim			= 0;
	p->iMix			= 0;
	p->dLdet			= NULL;
	p->dWt			= NULL;
	p->dLwt			= NULL;
	p->dtMu			= NULL;
	p->dtDCov		= NULL;
	p->dtFCov			= NULL;
	p->dMu			= NULL;
	p->dDCov			= NULL;
	p->dDInv			= NULL;
	p->dFCov			= NULL;
	p->dFCov			= NULL;
	p->dTrans			= NULL;

	minLogExp = -log(-LOGZERO);
	dsqrt2pi = sqrt(2 * PI);
	dlog2pi = log(2 * PI);
	return 1;
}

int Release(CGMM *pCGMM)
{
	CleanUpMdl();
	CleanUpPar();
	return 1;
}

//--------------------------------------------------------------
// Load GMM from the definition files
// RETURN VALUE:
//	-1 : Error opening file
//	1 : Success
int LoadGMMSet(char *fname)
{
	char *cTags[TAG_CUS_NUM] = {
			"NAME",		// 0
			"VECSIZE",	// 1
			"DIAG",		// 2
			"NUMMIXES",	// 3
			"MIXTURE",	// 4
			"MEAN",		// 5
			"VARIANCE",	// 6
			"ENDGMM",	// 7
			"FULL" 		// 8
	};

	char ctag[MAX_CUS_TAG];
	int iflag, i, j;
	FILE *fp;
	bool bNAME		= false;
	bool bVECSIZE		= false;
	bool bCOVTYPE		= false;
	bool bNUMMIXES	= false;
	int iCurMix=0;
	int  iMEAN			= -1;
	int  iVAR			= -1;
	int iD				= sizeof(double);
	
	if ((fp = fopen(fname,"rt")) == NULL)
		return -1;

	while ( (iflag = GetNextTag(fp,ctag)) == 1)
	{
		for (i=0; i < TAG_CUS_NUM; i++) {
			if ( strcmp(ctag,cTags[i]) == 0 )
				break;
		}
		switch (i) {
			case 0:		// <NAME>
				if (bNAME) { sprintf(msg,"Duplicated TAG <%s> in %s",ctag,fname); SayBye(msg); }

				fscanf(fp,"%s", p->cName);	// read name
				fscanf(fp,"%s", p->cParm);	// read parameter type

				bNAME = true;
				break;
			case 1:		// <VECSIZE>
				if (bVECSIZE) { sprintf(msg,"Duplicated TAG <%s> in %s",ctag,fname); SayBye(msg); } 
				if (!bNAME)	{ sprintf(msg,"Missing TAG <NAME> in %s",fname); SayBye(msg); }

				fscanf(fp,"%d",&p->iDim);
				if (p->iDim<=0) SayBye("Vector size should be positive");

				bVECSIZE = true;
				break;
			case 2:		// <DIAG>
				if (bCOVTYPE) { sprintf(msg,"Duplicated TAG <%s> in %s",ctag,fname); SayBye(msg); } 
				if (!bVECSIZE) { sprintf(msg,"Missing TAG <VECSIZE> in %s",fname); SayBye(msg); }
				
				p->bDiag = true;
				bCOVTYPE = true;
				break;
			case 3:		// <NUMMIXES>
				if (bNUMMIXES) { sprintf(msg,"Duplicated TAG <%s> in %s",ctag,fname); SayBye(msg); } 
				if (!bCOVTYPE) { sprintf(msg,"Missing TAG <DIAG> in %s",fname); SayBye(msg); }

				fscanf(fp,"%d",&p->iMix);
				if (p->iMix<=0) SayBye("Mixture of Gaussian should be positive");
				
				
				p->dLdet 		= (double *) calloc( p->iMix, iD );
				p->dWt		= (double *) calloc( p->iMix, iD );
				p->dLwt		= (double *) calloc( p->iMix, iD );
				p->dVfl		= (double *) calloc( p->iDim, iD);
				p->dVec 		= (double *) calloc( p->iDim, iD );
				p->dMu 		= (double **) alloc2d(p->iMix, p->iDim, iD);
				p->dtMu 		= (double **) alloc2d(p->iMix, p->iDim, iD);
				if (p->bDiag) {
					p->dDCov 	= (double **) alloc2d( p->iMix, p->iDim, iD);
					p->dDInv 	= (double **) alloc2d( p->iMix, p->iDim, iD);
					p->dtDCov= (double **) alloc2d( p->iMix, p->iDim, iD);
				}
				else {
					p->dFCov 	= (double ***) alloc3d( p->iMix, p->iDim, p->iDim, iD);
					p->dFInv 	= (double ***) alloc3d( p->iMix, p->iDim, p->iDim, iD);
					p->dtFCov	= (double ***) alloc3d( p->iMix, p->iDim, p->iDim, iD);
				}

				bNUMMIXES = true;
				break;
			case 4:		// <MIXTURE> - read in the value one by one
				if (!bNUMMIXES) { sprintf(msg,"Missing TAG <NUMMIXES> in %s",fname); SayBye(msg); }

				fscanf(fp,"%d", &iCurMix);
				fscanf(fp,"%lf", p->dWt + iCurMix);
				p->dLwt[iCurMix] = log(p->dWt[iCurMix]);
				if ( iCurMix!=(iMEAN+1) || iCurMix!=(iVAR+1) )
					{ sprintf(msg,"Wrong <%s> index in %s",ctag,fname); SayBye(msg); } 
				break;
			case 5:		// <MEAN>
				iMEAN ++;
				if (iMEAN != iCurMix) { sprintf(msg,"Wrong <MIXTURE> index in %s",fname); SayBye(msg); } 

				for (i=0; i< p->iDim; i++)
					fscanf(fp,"%lf",&p->dMu[iCurMix][i]);
				break;
			case 6:		// <VARIANCE>
				iVAR ++;
				if (iVAR != iCurMix) { sprintf(msg,"Wrong <MIXTURE> index in %s",fname); SayBye(msg); } 

				p->dLdet[iCurMix] = .0;
				if (p->bDiag) {
					for (i = 0; i < p->iDim; i++) {
						fscanf(fp,"%lf",&p->dDCov[iCurMix][i]);
						p->dDInv[iCurMix][i] = 1.0 / p->dDCov[iCurMix][i];
						p->dLdet[iCurMix] += log(p->dDCov[iCurMix][i]);
					}
				} else {
					for (i = 0; i < p->iDim; i++) {
						for (j = i; j < p->iDim; j++) {
							fscanf(fp, "%lf", &p->dFCov[iCurMix][i][j]);
							p->dFCov[iCurMix][j][i] = p->dFCov[iCurMix][i][j];
						}
					}
					p->dLdet[iCurMix] = InvertCov(p->dFCov[iCurMix], p->dFInv[iCurMix], p->iDim);
				}
				break;
			case 7:		// <ENDGMM>
				break;
			case 8:		// <FULL>
				if (bCOVTYPE) { sprintf(msg,"Duplicated TAG <%s> in %s",ctag,fname); SayBye(msg); } 
				if (!bVECSIZE) { sprintf(msg,"Missing TAG <VECSIZE> in %s",fname); SayBye(msg); }
				
				p->bDiag = false;
				bCOVTYPE = true;
				break;
			default:
				sprintf(msg,"Unknown TAG \"%s\"",ctag);
				LookOut(msg);
		}
	}
	if (iflag == 0) {
		sprintf(msg,"Unsolved TAG in \"%s\"",fname);
		SayBye(msg);
	}
	fclose(fp);
	// check integraty
	if ( ! (bNAME					// found <NAME> tag
		&& bVECSIZE					// found <VECSIZE> tag
		&& bCOVTYPE					// found <DIAG> tag
		&& bNUMMIXES				// found <NUMMIXES> tag
		&& (iMEAN==p->iMix-1)	// means loaded
		&& (iVAR ==p->iMix-1)) )	// covariances loaded
	{
		sprintf(msg,"Format Error in %s",fname);
		SayBye(msg);
	}

	p->bMLoaded = true;

	return 1;
}


//-------------------------------------------------------------------
// Get Next Tag from 'fp', tag stored in 'ctag'
// Return Value:
//      0 : Unresolved TAG
//      1 : Success
//     -1 : End of file
int GetNextTag(FILE *fp, char *ctag)
{
	int i;
	char ch;

	while (!feof(fp))
	{
		ch=getc(fp);
		if (ch!='<')
			continue;
		// ch = '<'
		i=0;
		while ((ctag[i]=getc(fp))!='>' && i<MAX_CUS_TAG) i++;
		if (i>=MAX_CUS_TAG)
			return 0;
		ctag[i]=0;
		return 1;
	}
	return -1;
}


#define FILE_BLOCK 500
int LoadParaList(FILE *fp, int ifilenum)
{	
	int i,j,k,inum,idim,isampledim;
	FILE * fptemp;
	HTKheader strHeader;
	
	char filename[500];

	int ilen, itmp, iframenum, icur, isample; 	//maybe long long
	long lpos;
	float ** buf;

	idim = p->iDim;
	
	//time_t ltime;
	//time( &ltime );
//	printf( "Loading start:\t\t%s", ctime( &ltime ) );	
	// open the file list
	
	lpos = ftell(fp);
	iframenum = 0;
	for (i = 0; i < ifilenum; i++) {
		if ( fscanf(fp, "%s\n",  filename) == EOF ){
			break;
		}
		if ( (fptemp = fopen(filename, "rb")) == NULL ) {
			printf("fail to open file %s %d\n", filename, i ); return -1;
		}
		if ( (LONG_MAX - iframenum) < 1e3  ) {
			printf("iframenum - %l overflow!\n", iframenum ); return -1;
		}
		
		fread(&strHeader, sizeof(HTKheader), 1, fptemp);
		ilen  = fseek(fptemp, 0L, SEEK_END);
		ilen = ftell(fptemp);
		fclose(fptemp);
		
		itmp = (int) ( ilen - sizeof(HTKheader) ) / ( idim * sizeof(float) );
		if (strHeader.nSamples != itmp) {
			printf("feature %s len and header mismatch.\n", filename); return -1;}
		isampledim = (int) strHeader.sampleSize/sizeof(float);
		if (isampledim != idim) {
			printf("feature dim %d and model dim %d mismath.\n", isampledim, idim); return -1;}
		iframenum += itmp;
	}
	if (ifilenum >= MAX_FILES){
		SayBye("Maximum files in file list is set to 20000");
	}

	// allocate a big buffer for reading all the feature files
	buf = (float **) alloc2d(iframenum, idim, sizeof(float));
	
	fseek(fp, lpos, SEEK_SET);
	icur =0;
	for (i = 0; i < ifilenum; i++)
	{
		fscanf(fp, "%s",  filename);
		if ( (fptemp = fopen(filename, "rb")) == NULL ) {
			printf("fail to open file %s %d\n", filename, i ); return -1;
		}
		fread(&strHeader, sizeof(HTKheader), 1, fptemp);
		isample = strHeader.nSamples;
		for (j = 0; j < isample; j++)
		{
			fread(buf[icur], sizeof(float), idim, fptemp);
			icur++;
		}
		fclose(fptemp);
	}
	p->fBuf = buf;
	//time( &ltime );
	//printf( "Loading OK:\t\t%s", ctime( &ltime ) );
	
	
	p->iFrm = iframenum;
	if (ifilenum!=0){
		p->bPLoaded = true;
	}

	return ifilenum;
}


// if the diag elememts is replaced by the var floor, update logwt, inv, logdet
int ChkVarFloor(bool bUpdateWVD)
{
	bool bFloor;
	int i, j;
	
	bFloor = false;
	for (i = 0; i < p->iMix; i++)
	{
		for (j = 0; j < p->iDim; j++) {
			if (p->bDiag) {
				if (p->dDCov[i][j] < p->dVfl[j]) {
					p->dDCov[i][j] = p->dVfl[j];
					bFloor = true; 
				} 
			}
			else {
				if (p->dFCov[i][j][j] < p->dVfl[j]) {
					p->dFCov[i][j][j] = p->dVfl[j];
					bFloor = true; 
				}
			}
		}
	}
	
	if (bFloor && bUpdateWVD)
		UpdateWVD();
	return 1;
}

//--------------------------------------------------------------
// Recalculate the Determinant for each mixture after the covariance matrixes are changed.
// plus update the log weight and the inverse covarance
int UpdateWVD()
{
	int i,j;
	
	for (i = 0; i < p->iMix; i++) {
		if (p->dWt[i] <= 0)
			p->dLwt[i] = LOGZERO;
		else
			p->dLwt[i] = log(p->dWt[i]);
		
		if (p->bDiag) {
			p->dLdet[i] = .0;
			for (j = 0; j < p->iDim; j++) {
				p->dDInv[i][j] = 1.0 / p->dDCov[i][j];
				p->dLdet[i] += log(p->dDCov[i][j]);
			}
		}
		else {
			p->dLdet[i] = InvertCov(p->dFCov[i], p->dFInv[i], p->iDim);
		}
	}
	
	return 1;
}

//---------------------------------------------------------------------
// Euclidean Distance
// (alternative: pre-calculate the total variance of each feature, 
// then weighting each feature with the variance)
double VQDistortion(float *vec1, float *vec2)
{
	int i;
	double result=0;
	for (i = 0; i < p->iDim; i++)
		result += (vec1[i]-vec2[i]) * (vec1[i]-vec2[i]);
	return result;
}


//-----------------------------------------------------------------
// Init the mean vectors of VQ procedure
// Return Value:
//   1 :  Success
//   0 :  Nothing to init
int VQInit(bool bMaxMin)
{
	int i, j, k, iMax, iMin;
	double dMax, dMin, dCur, *pdInter, da;
	// check if the model and the param have been loaded
	if (!p->bMLoaded || !p->bPLoaded){
		LookOut("Model or Param not loaded yet");return 0;
	}
	if (bMaxMin && p->iMix > 1) {
		dMax = DBL_MIN;
		dMin = DBL_MAX;
		for (i=0; i < p->iFrm; i++)	{
			dCur = .0;
			for (j=0; j < p->iDim; j++) {
				dCur += p->fBuf[i][j];
			}
			if (dCur < dMin) {
				dMin = dCur;
				iMin = i;
			}
			else if (dCur > dMax){
				dMax = dCur;
				iMax = i;
			}
		}
		
		for (i=0; i < p->iMix; i++) {
			for (j=0; j < p->iDim; j++) {
				da = (p->fBuf[iMax][j] - p->fBuf[iMin][j]) / (p->iMix - 1);
				p->dMu[i][j] = p->fBuf[iMin][j] + da * i;
			}
		}
	}
	else {
		int *iCount = (int *) calloc(p->iMix, sizeof(int));
		for (i = 0; i < p->iMix; i++) {
			iCount[i] = 0;
			for (j = 0; j < p->iDim; j++) {
				p->dMu[i][j]  = 0.0;
				if (p->bDiag) {
					p->dDCov[i][j] = 0.0;
				}
				else {
					for (k = 0; k < p->iDim; k++) {
						p->dFCov[i][j][k] = 0.0;
					}
				}
			}
		}

		for (i=0; i < p->iFrm; i++)	{
			for (j=0; j < p->iDim; j++) {
				p->dMu[i % p->iMix][j]  += p->fBuf[i][j];
				if (p->bDiag) {
					p->dDCov[i % p->iMix][j] += p->fBuf[i][j] * p->fBuf[i][j];
				}
				else {
					for (k = 0; k < p->iDim; k++) {
						p->dFCov[i % p->iMix][j][k] += p->fBuf[i][j] * p->fBuf[i][k];
					}
				}
				
			}
			iCount[i % p->iMix]++ ;
		}

		for (i = 0; i < p->iMix; i++) {
			for (j = 0; j < p->iDim; j++) {
				p->dMu[i][j]  /= iCount[i];
				if (p->bDiag) {
					p->dDCov[i][j] /= iCount[i];
					p->dDCov[i][j] -= p->dMu[i][j] * p->dMu[i][j];	
				}
				else	{
					for (k = 0; k < p->iDim; k++) {
						p->dFCov[i][j][k] /= iCount[i];;
						p->dFCov[i][j][k] -= p->dMu[i][j] * p->dMu[i][k];
					}
				}
			}
		}

		CalVarFloor();
		ChkVarFloor(true);
		
		free( iCount );
	}
	return 1;
}

//---------------------------------------------------------------
// Get the nearest mean vector to the given vector according to
// Euclidean distance.
//    fvector : Vector to calculate
//        err : minimal distance
int VQNearest(float *fvector, double *err)
{
	int i, j, k, imin	= 0;
	double dtmp	= 0.0;	// buffer to accumulate distances
	
	if (!p->bMLoaded || !p->bPLoaded) {
		LookOut("Model or Param not loaded yet");return 0;
	}	
	
	*err = DBL_MAX;
	for (i = 0; i < p->iMix; i++) {
		dtmp = 0.0;
		if (p->bDiag) {
			for (j = 0; j < p->iDim; j++) {
				dtmp += (fvector[j] - p->dMu[i][j] ) * (fvector[j] - p->dMu[i][j] ) * p->dDInv[i][j];
			}
		}
		else {
			for (j= 0; j < p->iDim; j++) {
				for (k= 0; k < p->iDim; k++) {
					dtmp += (fvector[j] - p->dMu[i][j] ) * (fvector[k] - p->dMu[i][k] ) * p->dFInv[i][j][k];
				}
			}
		}
		if (dtmp < *err) {
			*err = dtmp;
			imin = i;
		}
	}
	return imin;
}

// Get the variance floor
int CalVarFloor()
{
	int i, j;
	
	for (j = 0; j < p->iDim; j++) // for all frames
	{
		p->dVec[j] = 0;
		p->dVfl[j] = 0;
		for (i = 0; i < p->iFrm; i++)
		{
			p->dVec[j] += p->fBuf[i][j] / p->iFrm;
			p->dVfl[j] += p->fBuf[i][j] * p->fBuf[i][j] / p->iFrm;
		}
		p->dVfl[j] -= p->dVec[j] * p->dVec[j];
		p->dVfl[j] *= VAR_FLOOR;	// * 0.01
	}
	return 1;
}


//--------------------------------------------------------------------
// Init the model for EM by VQ
//   MaxLoop :  Max loop of VQ iterations
// Return Value:
//   1 : Success
int InitMdl(int MaxLoop)
{
	bool bFloor = false;
	int i, j, k, imin, sum;
	int ct = 0;					// loop times
	
	double dPrevAvgDist, dAvgDist; // 
	double dTh = 1e-5;		// for convergence
	double dtmp;
	
	int *iCount;			// counter of the assigned frames for each mixture
	int iframenum = p->iFrm;
	int idim = p->iDim;
	int imix = p->iMix;
	int bd = p->bDiag;	

	if (!p->bMLoaded || !p->bPLoaded)
		SayBye("Model or Param not loaded");
	if (p->iFrm / p->iMix < 30)
		SayBye("Too few training data [AVG < 30]");
	
	iCount	= (int *) calloc(imix, sizeof(int));

	// the current covariance matrix is loaded from the proto
	// make sure the proto is correct
	if (!VQInit(1))
		SayBye("VQInit not succeed");

	printf("\n\n - Initializing ... \n");
	dAvgDist = 0;
	
	if (MaxLoop <= 0)
		return 0;
	
	do {	
		//-- clean up tmp vectors
		for (i = 0; i < p->iMix; i++) {
			iCount[i] = 0;
			for (j = 0; j < p->iDim; j++) {
				p->dtMu[i][j] = .0;
				if (bd) {
					p->dtDCov[i][j] = .0;	
				}
				else {
					for (k = 0; k < p->iDim; k++) {
						p->dtFCov[i][j][k] = .0;
					}
				}
			}			
		}
		
		dPrevAvgDist = dAvgDist;
		dAvgDist = 0.0;
		for (i = 0; i < p->iFrm; i++) { // for all frames
			imin = VQNearest(p->fBuf[i], &dtmp);
			dAvgDist += dtmp / p->iFrm;
			
			// sum the err
			for (j = 0; j < p->iDim; j++) {
				p->dtMu[imin][j] += p->fBuf[i][j];
				if (bd) {
					p->dtDCov[imin][j] += p->fBuf[i][j] * p->fBuf[i][j];
				} 
				else {
					for (k = 0; k < p->iDim; k++) {
						p->dtFCov[imin][j][k] += p->fBuf[i][j] * p->fBuf[i][k];
					}					
				}
			}
			iCount[imin]++;
		}

		// recalculate the core
		for (i = 0; i < p->iMix; i++) {
			if (iCount[i] != 0) {
				for (j = 0; j < p->iDim; j++) {
					p->dMu[i][j] = p->dtMu[i][j] / iCount[i];
					if (bd) {
						p->dDCov[i][j] = p->dtDCov[i][j] / iCount[i] - p->dMu[i][j] * p->dMu[i][j];
					}
					else {
						for (k = 0; k < p->iDim; k++) {
							p->dFCov[i][j][k] = p->dtFCov[i][j][k] / iCount[i] - p->dMu[i][j] * p->dMu[i][k];
						}					
					}
				}
			}
		}
		
		//-- check for the non-assigned core
		// VQCheckZero(iCount);
		
		// check varFloor and update WVD
		ChkVarFloor(true);		
		printf("   - Loop: %d ,  Residual: %.6e,  Decreased : %.4f\%\n", ct, dAvgDist, -100*(dAvgDist - dPrevAvgDist) / (1e-5 + dPrevAvgDist) );
		ct++;
	}while( (fabs((dPrevAvgDist-dAvgDist)/(1e-5+dPrevAvgDist)) > dTh) && (ct < MaxLoop) );

	printf("Init - OK!\n");

	// Update weight, WVD
	for (i = 0; i < p->iMix; i++){
		p->dWt[i] = (double) iCount[i] / p->iFrm;
		//printf("%.3f ", p->dWt[i]);
	}
	UpdateWVD();
	
	free( iCount );
	return 1;
}

//-----------------------------------------------------------------
// EM algorithm
int ReesMdl(int MaxLoop)
{
	int i, j, k, l, iloop;
	double dLogFrmProb;
	double *dSumGama,*dLogMixWt,*dGama;
	double dAvgDist,dPrevAvgDist;
	const double dTh = 1e-5;
	double *dMixWt;		// A copy of mixture weights

	bool bd = p->bDiag;
	int imix = p->iMix;
	int idim = p->iDim;
	int iframenum = p->iFrm;
	
	if (MaxLoop <= 0)
		return 0;
	
	dMixWt 		= (double *) calloc(imix, sizeof(double));
	dLogMixWt	= (double *) calloc(imix, sizeof(double));
	dGama		= (double *) calloc(imix, sizeof(double));
	dSumGama	= (double *) calloc(imix, sizeof(double));
	
	// main procedure >>
	dAvgDist	= .0;
	iloop		= 0;
	printf("\n\n - Re-estimating ... \n");
	do {
		//-----------------------------------------------
		// E-step
		//-----------------------------------------------
		// clean up temporary values
		for (i = 0; i < p->iMix; i++)
		{			
			dMixWt[i] = dSumGama[i] = .0;
			for (j = 0; j < p->iDim; j++) {
				p->dtMu[i][j] = .0;
				if (bd) {
					p->dtDCov[i][j] = .0;	
				}
				else {
					for (k = 0; k < p->iDim; k++) {
						p->dtFCov[i][j][k] = .0;	
					}
				}
			}
		}
		
		dPrevAvgDist = dAvgDist;
		dAvgDist = .0;		
		for (i = 0; i< p->iFrm; i++)
		{
			double maxdGama;
			dGama[0]   = LMixProb(p->fBuf[i], 0);
			dGama[0] += p->dLwt[0];
			dLogFrmProb = maxdGama = dGama[0];
			for (j = 1; j < p->iMix; j++)
			{
				dGama[j]   = LMixProb(p->fBuf[i], j);
				dGama[j] += p->dLwt[j];
				dLogFrmProb = LogAdd(dGama[j], dLogFrmProb);
			}
			dAvgDist += dLogFrmProb;	// sum the dAvgDist
			
			for (j = 0; j < p->iMix; j++)
			{
				dGama[j] -= dLogFrmProb;
				dGama[j]  = exp( dGama[j] );

				// for updating weights
				dMixWt[j] += dGama[j];
				dSumGama[j]  += dGama[j];
				// for updating mean & covariance
				for (k = 0; k < p->iDim; k++)
				{					 
					p->dtMu[j][k] += dGama[j] * p->fBuf[i][k];
					if (bd) {
						p->dtDCov[j][k] += dGama[j] * p->fBuf[i][k] * p->fBuf[i][k];
					}
					else {
						for (l = 0; l < p->iDim; l++) {
							p->dtFCov[j][k][l] +=  dGama[j] * p->fBuf[i][k] * p->fBuf[i][l];
						}
					}
				}
			}
		}
		dAvgDist /= p->iFrm;
		
		//-----------------------------------------------
		// M-step
		//-----------------------------------------------
		// update weight, mean & covariance
		for (i = 0; i < p->iMix; i++)
		{
			p->dWt[i] = dMixWt[i] / p->iFrm;
			for (j = 0; j < p->iDim; j++)
			{
				p->dMu[i][j] = p->dtMu[i][j] / dSumGama[i];
				if (bd) {
					p->dDCov[i][j]  = p->dtDCov[i][j] / dSumGama[i];
					p->dDCov[i][j] -= p->dMu[i][j] * p->dMu[i][j];
				}
				else {
					for (k = 0; k < p->iDim; k++) {
						p->dFCov[i][j][k]  = p->dtFCov[i][j][k] / dSumGama[i];
						p->dFCov[i][j][k] -= p->dMu[i][j] * p->dMu[i][k];
					}
				}
			}
		}
		// update determinate
		ChkVarFloor(true);

		iloop ++;
		printf("   - Loop: %d ,  Avg Log Likelihood: %.6e,  Increment : %.4f\%\n", iloop, dAvgDist, 100*(dAvgDist - dPrevAvgDist) / fabs(dPrevAvgDist + 1e-5));
	} while( fabs((dAvgDist-dPrevAvgDist)/(dPrevAvgDist+1e-5))> dTh && iloop < MaxLoop);
	
	if (iloop >= MaxLoop)
		printf(" - Break at loop %d!\n",iloop);
	else
		printf(" - Converged at loop %d!\n",iloop);

	free( dMixWt );
	free( dGama );
	free( dSumGama );
	free( dLogMixWt );

	return 1;
}

//----------------------------------------------------------------
// Return the log likelihood of the given vector to the given
// mixture (the score should be multiplied with the mixture weight).
//    fparam :  pointer to the feature vector in float
//  iMix :  Mixture index
// Return Value:
//    The log-likelihood in double
double LMixProb(float *fparam, int iMix)
{
	int i, j;
	double dsum,dtmp;
	const double dfloor   = 0;

	if (!p->bMLoaded)
		return .0;

	dsum = .0;
	if (p->bDiag) {
		for (i = 0; i < p->iDim; i++) {
			dsum += (fparam[i] - p->dMu[iMix][i]) * (fparam[i] - p->dMu[iMix][i]) * p->dDInv[iMix][i];
		}
	}
	else {
		for (i = 0; i < p->iDim; i++) {
			for (j = 0; j < p->iDim; j++) {
				dsum = (fparam[i] - p->dMu[iMix][i]) * (fparam[j] - p->dMu[iMix][j]) * p->dFInv[iMix][i][j];
			}
		}	
	}
	dsum  = - p->iDim * dlog2pi - p->dLdet[iMix] - dsum;
	dsum *= 0.5;
	return (dsum<LOGZERO? LOGZERO:dsum);
}


//----------------------------------------------------------------
// Routine for adding two log-values in a linear scale, return the log
// result in double
double LogAdd(double logVar1,double logVar2)
{
	double temp,diff,z;
	
	if (logVar1<logVar2) {
		temp = logVar1; logVar1 = logVar2; logVar2 = temp;
	}
	diff = logVar2-logVar1;
	if (diff<minLogExp) 
		return  (logVar1<LSMALL)?LOGZERO:logVar1;
	else {
		z = exp(diff);
		return logVar1+log(1.0+z);
	}
}


//------------------------------------------------------------------
// Write GMM definition to *fname in ASC
// Return Value :
//    1 :  Success
//    0 :  No data to write
//   -1 :  error
int WriteGMMSet(char *fname, int imod)
{
	int i, j, k;
	FILE *fp;
	
	if (imod != 0) {
	// write model nums
		if( (fp = fopen(fname, "wt")) == NULL) {
			printf("cannot open %s for write.\n", fname); return -1;
		}
		fprintf(fp, "<MODNUM> %d\n\n", imod);
		fclose(fp);
	}
	else {
	// write a model into the file
		if (!p->bMLoaded)
			return 0;
		
		if ((fp=fopen(fname,"at+"))==NULL) {
			printf("cannot open %s for write.\n", fname);return -1;}
		
		//fprintf(fp, "<BEGINGMM>\n");
			fprintf(fp,"<NAME> %s %s\n", p->cName, p->cParm);
			fprintf(fp,"<VECSIZE> %d", p->iDim);
			if (p->bDiag)
				fprintf(fp," <DIAG>\n");
			else
				fprintf(fp," <FULL>\n");
			fprintf(fp, "<NUMMIXES> %d\n", p->iMix);
			for (i = 0; i < p->iMix; i++)
			{
				fprintf(fp, "<MIXTURE> %d %E\n", i, p->dWt[i]);
				fprintf(fp, " <MEAN> \n");
				for (j = 0; j < p->iDim; j++)
					fprintf(fp, "  %E", p->dMu[i][j]);
				fprintf(fp, "\n");
				fprintf(fp, " <VARIANCE> \n");
				if (p->bDiag) {
					for (j = 0; j < p->iDim; j++)
						fprintf(fp, "  %E", p->dDCov[i][j]);
					fprintf(fp,"\n");
				}
				else {
					for (j = 0; j < p->iDim; j++) {
						for (k = 0; k < p->iDim; k++)
							fprintf(fp, "  %E", p->dFCov[i][j][k]);
						fprintf(fp,"\n");
					}
					
				}
			}
			fprintf(fp, "<ENDGMM>\n\n");
			fclose(fp);
	}
	
	return 1;
}

//------------------------------------------------------------------
// Check if any mean vector (core) are not assigned any vectors, if
// yes, delete this one and divide the most populated one into two
// cores.
// Return Value:
//   1 : Success
//   0 : Model or Param not loaded yet
//  -1 : Not enough training data
/*
int VQCheckZero(int *counter)
{
	if (!p->bMLoaded || !p->bPLoaded){
		LookOut("Model or Param not loaded yet");return 0;}

	int i, j, k, opti, optc;
	static unsigned int nSeed = 0;
	double dBias = .0;
	double dPro = .0;							// proportion of smallest mixture in samples  
	float minimum; 							// threshold for smallest mixture
	double **dL = (double **)alloc2d(p->iDim, p->iDim, sizeof(double)); 
	
	minimum = p->iFrm / p->iMix * MIX_FLOOR;	// should be more than 2% of the average.
	for (i = 0; i < p->iMix; i++) {
		if (counter[i] <= minimum) {
			// find the largest mixture: opti
			optc = -1;
			for (j = 0; j < p->iMix; j++) {
				if (counter[j] > optc) {
					optc = counter[j];
					opti = j;
}
}

			
			// split the largest mixture, then update the smallest mixture:
			// small = a * small + (1-a) * large
			dPro = counter[i] / ( counter[opti] + counter[i] );
			if (p->bDiag) { 
				for (j = 0; j < p->iDim; j++) {
					// uj' = uj +- r*stdj					
					srand(nSeed++);								// do not worry about overflow
					dBias = (double)(rand()%10000)/10000. - .5;			// dFrac -0.5 ~ +0.5
					dBias *= sqrt(p->dDCov[opti][j]);					// std = sqrt of var !
					p->dMu[i][j]	= dPro * p->dMu[i][j] + (1 - dPro) * (p->dMu[opti][j] - dBias);
					p->dMu[opti][j]	+= dBias;
					p->dDCov[i][j]	= p->dDCov[opti][j];
}
}
			else {
				// u = u +- Lr
				for (j = 0; j < p->iDim; j++) {
					srand(nSeed++);
					p->dVec[k] = (double)(rand()%10000)/10000. - .5;
}
				Cholesky(p->dFCov[opti],  dL, p->iDim);
				for (k = 0; k < p->iDim; k++) {
					p->dMu[opti][k]	= p->dMu[opti][j] + dBias;
					// to be finished ..
}
}

			//update weights
			counter[i] += (int)(counter[opti]*0.5);
			counter[opti] -= (int)(counter[opti]*0.5);
}
}
	
	free( (char **) dL);
	return 1;
}
*/

/*
	FILE * fdb = fopen("debuginfo.txt", "wt");
	for (int i = 0; i < imix; i++)
{
		fprintf(fdb, "mix %d\n", i);
		for (int j = 0; j < idim; j++)
			fprintf(fdb, "%.3f\t", p->dMu[i][j]);
		fprintf(fdb, "\n");
		for (int j = 0; j < idim; j++)
			fprintf(fdb, "%.3f\t", p->dDCov[i][j]);
		fprintf(fdb, "\n\n");
}
	fclose(fdb);
*/
