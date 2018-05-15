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
//static CGMM *p;
static double *dgVfl;
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
int CleanUpMdl(CGMM *p)
{
	if (p->bMLoaded) {
		free( p->dWt  );
		free( p->dLwt );
		free( p->dVec);
		
		free2d( (char **) p->dMu );
		free2d( (char **) p->dtMu );
		if (p->bDiag) {
			free2d( (char **) p->dDCov );
			free2d( (char **) p->dDInv );
			free2d( (char **) p->dtDCov );
		}
		else {		
			free3d( (char ***) p->dFCov );
			free3d( (char ***) p->dFInv );
			free3d( (char ***) p->dtFCov );
		}
		
		p->bMLoaded = false;
	}
	return 1;
}

//-------------------------------------------------------------------
// Free the parameter's memoroy
int CleanUpPar(CGMM *p)
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
int Initialize(CGMM ***ppp, int imodnum)
{
	CGMM *p;
	int m;
	
	*ppp = (CGMM **) calloc(imodnum, sizeof(CGMM *));	//GMM num - imodnum
	for (m = 0; m < imodnum; m++) {
		// Initialization
	 	p = (CGMM *) malloc(sizeof(CGMM));
	 	(*ppp)[m] = p;
		p->bMLoaded		= false;
		p->bPLoaded		= false;
		p->bDiag		= true;
		p->fBuf			= NULL;
		p->iDim			= 0;
		p->iMix			= 0;
		p->dLdet		= NULL;
		p->dWt			= NULL;
		p->dLwt			= NULL;
		p->dtMu			= NULL;
		p->dtDCov		= NULL;
		p->dtFCov		= NULL;
		p->dMu			= NULL;
		p->dDCov		= NULL;
		p->dDInv		= NULL;
		p->dFCov		= NULL;
		p->dFCov		= NULL;
		p->dTrans		= NULL;
		p->dLPrior		= (double)log(1.0/imodnum); //(double)log(1.0/imodnum); //equal prior 1.0/imodnum
	}
	
	dgVfl = (double *) calloc(p->iDim, sizeof(double));
	minLogExp = -log(-LOGZERO);
	dsqrt2pi = sqrt(2 * PI);
	dlog2pi = log(2 * PI);
	return 1;
}

int Release(CGMM **pp, int imodnum)
{
	CGMM *p;
	int m;
	for (m = 0; m  < imodnum; m++) {
		p = pp[m];
		CleanUpMdl(p);
		CleanUpPar(p);	
	}
	//free( dgVfl );
	return 1;
}

//--------------------------------------------------------------
// Load GMM from the definition files
// the definination can be either proto or trained models
int LoadGMMSet(char *fname, CGMM ***pCGMM)
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
	
	CGMM *p 	= NULL;
	char ctag[MAX_CUS_TAG];
	int iflag, i, j, m;
	FILE *fp;
	bool bNAME, bVECSIZE, bCOVTYPE, bNUMMIXES, bNext;
	int iCurMix, iMEAN, iVAR; 
	int iD		= sizeof(double);
	int iModNum	= 0;
	
	if ((fp = fopen(fname,"rt")) == NULL)
		return -1;
	fscanf(fp, "<MODNUM> %d\n", &iModNum);
	
	for (m = 0; m < iModNum; m++) {
		p = (*pCGMM)[m];
		bNAME		= false;
		bVECSIZE	= false;
		bCOVTYPE	= false;
		bNUMMIXES	= false;
		bNext		= false;
		iCurMix		= 0;
		iMEAN		= -1;
		iVAR		= -1;
		while ( (iflag = GetNextTag(fp,ctag)) == 1) {
			
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
					
					
					p->dLdet 	= (double *) calloc( p->iMix, iD );
					p->dWt		= (double *) calloc( p->iMix, iD );
					p->dLwt		= (double *) calloc( p->iMix, iD );
					p->dVec 	= (double *) calloc( p->iDim, iD );
					p->dMu 		= (double **) alloc2d(p->iMix, p->iDim, iD);
					p->dtMu 	= (double **) alloc2d(p->iMix, p->iDim, iD);
					if (p->bDiag) {
						p->dDCov = (double **) alloc2d( p->iMix, p->iDim, iD);
						p->dDInv = (double **) alloc2d( p->iMix, p->iDim, iD);
						p->dtDCov= (double **) alloc2d( p->iMix, p->iDim, iD);
					}
					else {
						p->dFCov = (double ***) alloc3d( p->iMix, p->iDim, p->iDim, iD);
						p->dFInv = (double ***) alloc3d( p->iMix, p->iDim, p->iDim, iD);
						p->dtFCov= (double ***) alloc3d( p->iMix, p->iDim, p->iDim, iD);
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
	
					for (i = 0; i < p->iDim; i++)
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
							for (j = 0; j < p->iDim; j++) {
								fscanf(fp, "%lf", &p->dFCov[iCurMix][i][j]);
							}
						}
						p->dLdet[iCurMix] = InvertCov(p->dFCov[iCurMix], p->dFInv[iCurMix], p->iDim);
					}
					break;
				case 7:		// <ENDGMM>
					bNext = true;
					break;
				case 8:		// <FULL>
					if (bCOVTYPE) { sprintf(msg,"Duplicated TAG <%s> in %s",ctag,fname); SayBye(msg); } 
					if (!bVECSIZE) { sprintf(msg,"Missing TAG <VECSIZE> in %s",fname); SayBye(msg); }
					
					p->bDiag = false;
					bCOVTYPE = true;
					break;
				default:
					break;
					sprintf(msg,"Unknown TAG \"%s\"",ctag);
					LookOut(msg);
			}
			
			// start looking for the next gmm def
			if (bNext) {
				// check integraty
				if ( ! (bNAME				// found <NAME> tag
					&& bVECSIZE		// found <VECSIZE> tag
					&& bCOVTYPE			// found <DIAG> tag
					&& bNUMMIXES			// found <NUMMIXES> tag
					&& (iMEAN==p->iMix-1)		// means loaded
					&& (iVAR ==p->iMix-1)) )	// covariances loaded
				{
					sprintf(msg,"Format Error in %s",fname);
					SayBye(msg);
				}

				p->bMLoaded = true;
				break;
			}
		}
	}
	return 1;
}


//-------------------------------------------------------------------
// show GMM configuration
void ShowGMMCfg(CGMM *p, int ii)
{
	printf("\nGMM Training for <%s>\n", p->cName);
	printf(" Mixture		:  %d\n", p->iMix);
	if (p->bDiag )
		printf(" Cov Matrix Type	:  Diagonal\n");
	else
		printf(" Cov Matrix Type	:  Full\n");
	printf(" Vector Size		:  %d\n", p->iDim);
	printf(" Number of Frames	:  %d\n", p->iFrm);
		
	if (p->bDiag )
		sprintf(p->cName, "g.n%d.d%d.m%d.d.tr%d", (ii + 1), p->iDim, p->iMix, p->iFrm);
	else
		sprintf(p->cName, "g.n%d.d%d.m%d.d.tr%d", (ii + 1), p->iDim, p->iMix, p->iFrm);
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



//-------------------------------------------------------------------
// Load Feature into p->fBuf
// ifilenum : utterances to read from the current fp 
int LoadParaList(FILE *fp, int ifilenum, CGMM *p)
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
	//printf( "Loading start:\t\t%s", ctime( &ltime ) );	
	
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

	return 1;
}


// if the diag elememts is replaced by the var floor, update logwt, inv, logdet
int ChkVarFloor(bool bUpdateWVD, CGMM *p)
{
	bool bFloor;
	int i, j;
	
	bFloor = false;
	for (i = 0; i < p->iMix; i++)
	{
		for (j = 0; j < p->iDim; j++) {
			if (p->bDiag) {
				if (p->dDCov[i][j] < dgVfl[j]) {
					p->dDCov[i][j] = dgVfl[j];
					bFloor = true; 
				} 
			}
			else {
				if (p->dFCov[i][j][j] < dgVfl[j]) {
					p->dFCov[i][j][j] = dgVfl[j];
					bFloor = true; 
				}
			}
		}
	}
	
	if (bFloor)
		printf("\tcovariance matrix is floored\n");
	if (bFloor || bUpdateWVD)
		UpdateWVD(p);
	return 1;
}

//--------------------------------------------------------------
// Recalculate the Determinant for each mixture after the covariance matrixes are changed.
// plus update the log weight and the inverse covarance
int UpdateWVD(CGMM *p)
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
double VQDistortion(float *vec1, float *vec2, CGMM *p)
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
//   bMaxMin: 0 choose the very first points to initialize each mixture
//	      1 use points among extreme points to initialize each mixture
int VQInit(bool bMaxMin, CGMM **pp, int imod)
{
	int m, i, j, k, iMax, iMin;
	double dMax, dMin, dCur, *pdInter, da;
	CGMM *p;
	// check if the model and the param have been loaded
	for (m = 0; m < imod; m++) {
		p = pp[m];
		if (!p->bMLoaded || !p->bPLoaded){
			LookOut("Model or Param not loaded yet");return 0;
		}
	}
	if (bMaxMin && p->iMix > 1) {
		// init each model
		for (m = 0; m < imod; m++) {
			dMax = DBL_MIN;
			dMin = DBL_MAX;
			p = pp[m];
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
		
	}
	else {
		/*
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
		*/
		for (m = 0; m < imod; m++) {
			p = pp[m];
			for (i = 0; i < p->iMix; i++) {
				for (j = 0; j < p->iDim; j++) {
					p->dMu[i][j]  = p->fBuf[i][j];
				}
			}
		}		
		//free( iCount );
	}
	return 1;
}

//---------------------------------------------------------------
// Get the nearest mean vector to the given vector according to
// Euclidean distance.
//    fvector : Vector to calculate
//        err : minimal distance
int VQNearest(float *fvector, double *err, CGMM *p)
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

// Get the global variance floor
int CalVarFloor(CGMM **pp, int imod, int itrained)
{
	int i, j, m, iAllFrm;
	CGMM *p = pp[0];
	double *dVec = (double *) calloc(p->iDim, sizeof(double));
	for (j = 0; j < p->iDim; j++) {
		dVec[j] = .0;
		dgVfl[j] = .0;
	}
	iAllFrm = 0;
	for (m = 0; m < imod; m++) {
		iAllFrm += pp[m]->iFrm;
	}
	
	for (j = 0; j < p->iDim; j++) {
		for (m = 0; m < imod; m++) {
			p = pp[m];
			for (i = 0; i < p->iFrm; i++) {
				dVec[j] += p->fBuf[i][j] / iAllFrm;
				dgVfl[j] += p->fBuf[i][j] * p->fBuf[i][j] / iAllFrm;
			}
		}
		dgVfl[j] -= dVec[j] * dVec[j];
		
		//set the global variance
		if (itrained == 0) {
			// start from the proto
			if (p->bDiag) {
				for (i = 0; i < p->iMix; i++)
					p->dDCov[i][j] = dgVfl[j];  
			}
			else {
				for (i = 0; i < p->iMix; i++)
					p->dFCov[i][j][j] = dgVfl[j];
			}
		}
		dgVfl[j] *= VAR_FLOOR;	// * 0.01
	}
	
	free(dVec);
	return 1;
}


//--------------------------------------------------------------------
// Init the model for EM by VQ
//   MaxLoop :  Max loop of VQ iterations
int InitMdl(int MaxLoop, CGMM *p)
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
	if (p->iFrm / p->iMix < 10)
		SayBye("Too few training data [AVG < 10]");
	
	iCount	= (int *) calloc(imix, sizeof(int));

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
			imin = VQNearest(p->fBuf[i], &dtmp, p);
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
				}
				for (j = 0; j < p->iDim; j++) {
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
		printf("   - Loop: %d ,  Residual: %.6e,  Decreased : %.4f\%\n", ct, dAvgDist, -100*(dAvgDist - dPrevAvgDist) / (1e-5 + dPrevAvgDist) );
		ChkVarFloor(true, p);		
		ct++;
	}while( (fabs((dPrevAvgDist-dAvgDist)/(1e-5+dPrevAvgDist)) > dTh) && (ct < MaxLoop) );

	printf("Init - OK!\n");

	// Update weight, WVD
	for (i = 0; i < p->iMix; i++){
		p->dWt[i] = (double) iCount[i] / p->iFrm;
		//printf("%.3f ", p->dWt[i]);
	}
	UpdateWVD(p);
	
	free( iCount );
	return 1;
}

//-----------------------------------------------------------------
// MLE re-estimate
int ReesMdl(int MaxLoop, CGMM *p)
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
			dLogFrmProb = LOGZERO;		// log scale probability of a given frame
			for (j = 0; j < p->iMix; j++)
			{
				dGama[j]   = LMixProb(p->fBuf[i], p, j);
				dGama[j] += p->dLwt[j];
				dLogFrmProb = LogAdd(dGama[j], dLogFrmProb);
			}
			dAvgDist += dLogFrmProb/p->iFrm;	// sum the dAvgDist
			
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
		
		//-----------------------------------------------
		// M-step
		//-----------------------------------------------
		// update weight, mean & covariance
		for (i = 0; i < p->iMix; i++)
		{
			p->dWt[i] = dMixWt[i] / p->iFrm;
			for (j = 0; j < p->iDim; j++) {
				p->dMu[i][j] = p->dtMu[i][j] / dSumGama[i];
			}
			for (j = 0; j < p->iDim; j++) {
				if (bd) {
					p->dDCov[i][j]  = p->dtDCov[i][j] / dSumGama[i] - p->dMu[i][j] * p->dMu[i][j];
				}
				else {
					for (k = 0; k < p->iDim; k++) {
						p->dFCov[i][j][k]  = p->dtFCov[i][j][k] / dSumGama[i] - p->dMu[i][j] * p->dMu[i][k];
					}
				}
			}
		}
		// update determinate
		ChkVarFloor(true, p);

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


int ReesMdlMMI(int MaxLoop, CGMM **pp, int imod)
{
	int m, n, i, j, k, l, it, iloop, iAllFrm;
	double dAvgDist,dPrevAvgDist, dDSm, dCSm, dt1, dt2;
	double dLProb, dNewGama, dPostProb;
	double *dLJoinProb, *dLPostProb, *dNextMixWt, *dCurrMixWt, *dKMix;
	double **dGama, **dCondGama, **dJoinGama;
	double **dNumGama, **dDenGama, **dSumGama;	
	const double dTh = 1e-4;
	const double dDSmFac = -2.0;
	const double dCSmFac = 2.0;
	const int iWtIter = 100;
	const int iDan = 1;
	CGMM *q, *p;
	
	if (MaxLoop <= 0)
		return 0;

	p = pp[0];
	bool bd  = p->bDiag;
	int imix = p->iMix;
	int idim = p->iDim;
	
	dCurrMixWt	= (double *) calloc(imix, sizeof(double));
	dNextMixWt	= (double *) calloc(imix, sizeof(double));
	dKMix		= (double *) calloc(imix, sizeof(double));
	dLJoinProb	= (double *) calloc(imod, sizeof(double));
	dLPostProb	= (double *) calloc(imod, sizeof(double));
	dSumGama	= (double **) alloc2d(imod, imix, sizeof(double));
	dGama		= (double **) alloc2d(imod, imix, sizeof(double));	
	dNumGama	= (double **) alloc2d(imod, imix, sizeof(double));
	dDenGama	= (double **) alloc2d(imod, imix, sizeof(double));
	dCondGama	= (double **) alloc2d(imod, imix, sizeof(double));
	dJoinGama	= (double **) alloc2d(imod, imix, sizeof(double));
	
	iAllFrm = 0;
	for (m = 0; m < imod; m++) {
		iAllFrm += pp[m]->iFrm;
	}
	// main procedure >>
	dAvgDist	= .0;
	iloop		= 0;
	printf("\n\n - Re-estimating ... \n");
	do {
		//-----------------------------------------------
		// E-step
		//-----------------------------------------------
		// clean up temporary values
		for (m = 0; m < imod; m++) {
			p = pp[m];
			for (i = 0; i < p->iMix; i++) {	
				dSumGama[m][i] = .0; 
				dNumGama[m][i] = .0;
				dDenGama[m][i] = .0;
				for (j = 0; j < p->iDim; j++) {
					p->dtMu[i][j] = .0;
					if (bd) {
						p->dtDCov[i][j] = .0;	
					} else {
						for (k = 0; k < p->iDim; k++) {
							p->dtFCov[i][j][k] = .0;	
						}
					}
				}
			}
		}
		
		dPrevAvgDist = dAvgDist;
		dAvgDist = .0;
		for (m = 0; m < imod; m++) {
			p = pp[m];
			for (i = 0; i < p->iFrm; i++) {
				//every frame, need clear
				for (n = 0; n < imod; n++) {
					dLJoinProb[n] = LOGZERO;		// log scale probability of a given frame
					for (j = 0; j < p->iMix; j++) {
						dGama[n][j] = .0;
						dCondGama[n][j] = dJoinGama[n][j] = .0;
					}
				}		

				dLProb = LOGZERO;
				for (n = 0; n < imod; n++) {
					q = pp[n];
					// logP(Ot|Mt') = LOGSUM [ logP(Ot|ck, Mt') + logP(ck|Mt') ]
					for (j = 0; j < q->iMix; j++) {
						dGama[n][j] = LMixProb(p->fBuf[i], q, j) + q->dLwt[j];
						dLJoinProb[n] = LogAdd(dGama[n][j], dLJoinProb[n]);
					}
					// logGamma(ck|Mt') = logp(ck|Mt',Ot) = logP(Ot,Ck|Mt') - logP(Ot|Mt')
					for (j = 0; j < q->iMix; j++) {
						dCondGama[n][j] = dGama[n][j] - dLJoinProb[n];
					}
					// logP(Ot, Mt') = logP(Ot|Mt') + logP(Mt');
					dLJoinProb[n] += pp[n]->dLPrior;
					// logP(Ot) = LOGSUM logP(Ot, Mt') 
					dLProb = LogAdd(dLJoinProb[n], dLProb);
				}

				//only after going through all the models, the calculation of dLProb is over.
				for (n = 0; n < imod; n++) {
					//logP(Mt'|Ot) = logP(Ot,Mt') - logP(Ot)
					dLPostProb[n] = dLJoinProb[n] - dLProb;
				}
				//Dist = SUM_t logP(Mt|Ot)
				dAvgDist += dLPostProb[m] / iAllFrm;
				//really postprob of current frame
				dPostProb = exp(dLProb);
				
				for (n = 0; n < imod; n++) {
					q = pp[n];
					for (j = 0; j < p->iMix; j++) {
						//logGamma(ck, Mt') = logP(ck,Mt'|Ot) = logP(ck|Mt',Ot) + logP(Mt'|Ot)
						//                 = logGamma(ck|Mt') + logP(Mt'|Ot)
						dJoinGama[n][j] = dCondGama[n][j] + dLPostProb[n];
						dt1 = exp(dCondGama[n][j]);
						dt2 = exp(dJoinGama[n][j]);
						dNewGama =  dt1 - dt2;
						dSumGama[n][j] += dNewGama;
						dNumGama[n][j] += dt1;
						dDenGama[n][j] += dt2;
						// for updating mean & covariance
						for (k = 0; k < p->iDim; k++)
						{					 
							q->dtMu[j][k] += dNewGama * q->fBuf[i][k];
							if (bd) {
								q->dtDCov[j][k] += dNewGama * q->fBuf[i][k] * q->fBuf[i][k];
							}
							else {
								for (l = 0; l < p->iDim; l++) {
									q->dtFCov[j][k][l] +=  dNewGama * q->fBuf[i][k] * q->fBuf[i][l];
								}
							}
						}
					}
				}
			}
		}
			
		//-----------------------------------------------
		// M-step
		//-----------------------------------------------
		// update weight, mean & covariance
		for (m = 0; m < imod; m++) {
			p = pp[m];
			for (i = 0; i < p->iMix; i++) {
				// find the D - smoothing factor
				dDSm = dDSmFac * dDenGama[m][i];
				
				// update mean and cov
				for (j = 0; j < p->iDim; j++) {
					if (bd) {
						dt1 = p->dtDCov[i][j] + dDSm*(p->dDCov[i][j] + p->dMu[i][j]*p->dMu[i][j]);
						p->dDCov[i][j] = dt1 / (dSumGama[m][i] + dDSm);
					} else {
						for (k = 0; k < p->iDim; k++) {
							dt1 = p->dtFCov[i][j][k] + dDSm*(p->dFCov[i][j][k] + p->dMu[i][j]*p->dMu[i][k]);
							p->dFCov[i][j][k]  = dt1 / (dSumGama[m][i] + dDSm);
						}
					}
				}
				for (j = 0; j < p->iDim; j++) {
					p->dMu[i][j] = (p->dtMu[i][j] + dDSm*p->dMu[i][j]) / (dSumGama[m][i] + dDSm);
				}
				for (j = 0; j < p->iDim; j++) {
					if (bd) {
						p->dDCov[i][j] -= p->dMu[i][j] * p->dMu[i][j];
					} else {
						for (k = 0; k < p->iDim; k++) {
							p->dFCov[i][j][k] -= p->dMu[i][j] * p->dMu[i][k];
						}
					}
				}
			}
		}
			
		// update weight
		if (iDan == 1) {
			// find max_k gamma^den /P(c'k)
			for (m = 0; m < imod; m++) {
				dt1 = DBL_MIN;
				p = pp[m];
				for (i = 0; i < p->iMix; i++) {
					if (dDenGama[m][i]/p->dWt[i] > dt1) {
						dt1 = dDenGama[m][i]/p->dWt[i];
					}
					dCurrMixWt[i] = p->dWt[i];
				}
				for (i = 0; i < p->iMix; i++) {
					dKMix[i] = dt1 - dDenGama[m][i]/p->dWt[i];
				}
				
				for (it = 0; it < iWtIter; it++) {
					dt1 = .0;
					for (i = 0; i < p->iMix; i++) {
						dNextMixWt[i] = dKMix[i]*dCurrMixWt[i] + dNumGama[m][i];
						dt1 += dNextMixWt[i];
					}
					for (i = 0; i < p->iMix; i++) {
						dNextMixWt[i] /= dt1;
						dCurrMixWt[i] = dNextMixWt[i]; 
					}	
				}
				
				for (i = 0; i < p->iMix; i++) {
					p->dWt[i] = dCurrMixWt[i];
				}
			}
		} else {
			for (m = 0; m < imod; m++) {
				p = pp[m];
				dt1 = dt2 = .0;
				dCSm = DBL_MIN;				
				for (i = 0; i < p->iMix; i++) {
					dt1 += dNumGama[m][i];
					dt2 += dDenGama[m][i];
				}
				for (i = 0; i < p->iMix; i++) {
					dNumGama[m][i] /= dt1;
					dDenGama[m][i] /= dt2;
					// max dDenGama[m][i]
					if (dCSm < dCSmFac * dDenGama[m][i] ) {
						dCSm = dCSmFac * dDenGama[m][i];
					}
				}
				dt1 = .0;
				for (i = 0; i < p->iMix; i++) {
					p->dWt[i] = p->dWt[i] * (dNumGama[m][i] - dDenGama[m][i] + dCSm);
					dt1 += p->dWt[i];
				}
				for (i = 0; i < p->iMix; i++) {
					p->dWt[i] /= dt1;
				}
			}
		}
			
		// update determinate
		ChkVarFloor(true, p);

		iloop ++;
		printf("   - Loop: %d ,  Avg Log Likelihood: %.6e,  Increment : %.4f\%\n", iloop, dAvgDist, 100*(dAvgDist - dPrevAvgDist) / fabs(dPrevAvgDist + 1e-5));
	} while( fabs((dAvgDist-dPrevAvgDist)/(dPrevAvgDist+1e-5))> dTh && iloop < MaxLoop);
	
	if (iloop >= MaxLoop)
		printf(" - Break at loop %d!\n",iloop);
	else
		printf(" - Converged at loop %d!\n",iloop);

	free( dKMix );
	free( dCurrMixWt );
	free( dNextMixWt );
	free( dLJoinProb );
	free( dLPostProb );
	free2d( (char **)dGama );
	free2d( (char **)dSumGama );
	free2d( (char **)dNumGama );
	free2d( (char **)dDenGama );
	free2d( (char **)dJoinGama );
	free2d( (char **)dCondGama );
	
	return 1;
}

//----------------------------------------------------------------
// Return the log likelihood of the given vector to the given
// mixture (the score should be multiplied with the mixture weight).
//    fparam :  pointer to the feature vector in float
//  iMix :  Mixture index
// Return Value:
//    The log-likelihood in double
double LMixProb(float *fparam, CGMM *p, int iMix)
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
				dsum += (fparam[i] - p->dMu[iMix][i]) * (fparam[j] - p->dMu[iMix][j]) * p->dFInv[iMix][i][j];
			}
		}	
	}
	dsum  = - p->iDim * dlog2pi - p->dLdet[iMix] - dsum;
	dsum *= 0.5;
	return (dsum<LSMALL? LOGZERO:dsum);
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
// imod == 0, write GMM definition in p to *fname
// imod != 0, write num of GMM definition to *fname
int WriteGMMSet(char *fname, int imod, CGMM *p)
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
