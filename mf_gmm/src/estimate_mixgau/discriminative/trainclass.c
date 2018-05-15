/*---------------------------------------------------------------------------
* Program to train a Gaussian mixture distribution for one or more
* classes discriminatively.
* Initialization is done with ML trained densities
*----------------------------------------------------------------------------*/
#include "discriminative.h"

int main(int argc, char **argv)
{
    int nclass, *nmodes, dim, ndclass, *nseqs, **seqlen;
    float ***mean, ***var, **c, ***newmean, ***newvar, **newc, ****data;
    float *classpriors;
    char *distctllist, *priorsfile=NULL, *datactllist, *outdistctllist;

    if (argc != 4 && argc != 5)
        QUIT(("USAGE: %s datactllist outdistctllist distctllist [priorsfile]\n",
              argv[0]))
    datactllist = argv[1]; outdistctllist = argv[2]; distctllist = argv[3];
    priorsfile = argc == 5 ? argv[4] : NULL;

    // Read all distributions for all the classes
    readalldistributions(distctllist,&nclass,&nmodes,&dim,&c,&mean,&var);

    // Read priors for all classes (set to uniform if reqd)
    readpriors(priorsfile,nclass,&classpriors);

    // Read data for all classes
    readallclassdata(datactllist, dim, &data, &ndclass, &nseqs, &seqlen);
    if (ndclass != nclass) 
        QUIT(("%d distributions, data for only %d classes\n",nclass, ndclass))

    // Reestimate all models
    reestimateallmodels(data,nclass,nseqs,seqlen,dim,mean,var,c,nmodes,
                        &newmean,&newvar,&newc,classpriors);
		
    // Write out models
    storealldistributions(outdistctllist,nclass,nmodes,dim,newc,newmean,newvar);

    return 0;
}
