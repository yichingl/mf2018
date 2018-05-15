/***************************************************************************
 *   Copyright (C) 2008 by Wei Chu @ UCLA                                         *
 *                                                                                                      *
 *   Function :  GMM Training                                                              *
 *   Created  :  06/18/2008                                                                  *
 *   Modified :                                                                                     *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>

#include "pconf.h"
#include "srTrainEx.h"
#include "gmmtr.h"


char msg[1024];

void PrintUsage();



extern FILE *fdb;

int main(int argc, char *argv[])
{
	int ii,	ifnum,iIter1 = 5,iIter2 = 10, imodnum;

	char * PathProto	= "seed.pro";
	char * PathOutput	= "model.def";
	char * trList		= "list";
	CGMM *p;

	char fn[512];
	char strCmd[512];
	char PathTrainList[512];
	FILE *fp, fpmod;
	static Confit params[] = {
		
		"Proto Model", "Model for GMM training to start with",
		"-i", cmuSTRING, (caddr_t)&PathProto,

		"Train List", "list of the training file list: each line - GMM proto path \t GMM training file list",
		"-l", cmuSTRING, (caddr_t)&trList,

		"Output Model", "Output Dir of the GMM",
		"-o", cmuSTRING, (caddr_t)&PathOutput,
		
		"Iter1", "Max iterations when initializing the model",
		"-i1", cmuINT, (caddr_t)&iIter1,

		"Iter2", "Max iterations when re-estimating the model",
		"-i2", cmuINT, (caddr_t)&iIter2,

		0, 0, 0, NOTYPE, 0
	};

	if (argc < 5)
	{
		pusage (argv[0], params);
		exit(0);
	}

	if ( ppconf( argc, argv, params, true ) ) 
	{
		pusage(argv[0], params);
		exit(1);
	}

	if( (fp = fopen(trList, "rt")) == NULL) {
		printf("cannot open %s for read.\n", trList); return -1;
	}
	fscanf(fp, "%d\n", &imodnum);
	if (imodnum <= 0) {
		printf("modnum should be larger than 1\n"); return -1;
	}
	WriteGMMSet(PathOutput, imodnum);
	for (ii = 0; ii < imodnum; ii++)
	{
		// Initialization
		Initialize(&p);
	
		// Load GMM set
		if (LoadGMMSet(PathProto) == -1) {
			sprintf(msg,"Can't open \"%s\" for reading", PathProto);	SayBye(msg);}

		// Load Features
		fscanf(fp, "%d\n", &ifnum);	
		ifnum = LoadParaList(fp, ifnum);
		if ( ifnum == -1) {
			sprintf(msg,"Can't open \"%s\" for reading", PathTrainList);	SayBye(msg);
		} else if (ifnum == 0) {
			SayBye("No training files loaded");
		}
	
		printf("\nGMM Training for <%s>\n", p->cName);
		printf(" Mixture		:  %d\n", p->iMix);
		if (p->bDiag )
			printf(" Cov Matrix Type	:  Diagonal\n");
		else
			printf(" Cov Matrix Type	:  Full\n");
		printf(" Vector Size		:  %d\n", p->iDim);
		printf(" Number of Files	:  %d\n", ifnum);
		printf(" Number of Frames	:  %d\n", p->iFrm);
		
		if (p->bDiag )
			sprintf(p->cName, "g.n%d.d%d.m%d.d.tr%d", (ii + 1), p->iDim, p->iMix, p->iFrm);
		else
			sprintf(p->cName, "g.n%d.d%d.m%d.d.tr%d", (ii + 1), p->iDim, p->iMix, p->iFrm);

		//K-Means
		InitMdl(iIter1);
		// re-estimate model EM
		ReesMdl(iIter2);
	
		// Write the GMM set
		WriteGMMSet(PathOutput, 0);

			
		if (Release(p) == 0) {
			printf("Release model and feature array failed\n"); 
		};
	}
		
	printf("Training Finished. :-)\n\n");
	
	return 0;
}


void PrintUsage(void)
{
	printf("                                                                   \n");
	printf("lGmmTr: GMM-Training Utility.                                      \n");
	printf("                                                                   \n");
	printf(" [USAGE]: lGmmTr <ProtoFile> <OutFile> <Tr_File_List> [Ite1] [Ite2]\n");
	printf("                                                                   \n");
	printf("     <ProtoFile> :  Proto definition of the GMM                    \n");
	printf("       <OutFile> :  Output file name of the GMM                    \n");
	printf("  <Tr_File_List> :  HTK format training file list                  \n");
	printf("          [Ite1] :  Max iterations when initializing the model     \n");
	printf("          [Ite2] :  Max iterations when re-estimating the model    \n");
	printf("                                                                   \n");
	printf("                                                                   \n");
}