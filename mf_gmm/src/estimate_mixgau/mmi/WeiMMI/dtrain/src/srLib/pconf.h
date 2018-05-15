/* PCONF.H
 *------------------------------------------------------------*
 * Copyright 1988, Fil Alleva and Carnegie Mellon University.
 *------------------------------------------------------------*
 * HISTORY
 */

#ifndef _PCONF_
#define _PCONF_

#define TRUE 1
#define FALSE 0

typedef enum {NOTYPE,
   cmuBYTE, cmuSHORT, cmuINT, cmuLONG,
   U_BYTE, U_SHORT, U_INT, U_LONG,
   cmuFLOAT, cmuDOUBLE,
   cmuBOOL, cmuCHAR, cmuSTRING
} art;

typedef char* caddr_t;
typedef unsigned char u_char;
typedef unsigned short u_short;
typedef unsigned int u_int;
typedef unsigned long u_long;


typedef union _ptr {
   char *CharP;
   char *ByteP;
   u_char  *UByteP;
   short   *ShortP;
   u_short *UShortP;
   int  *IntP;
   u_int   *UIntP;
   long *LongP;
   u_long  *ULongP;
   float   *FloatP;
   double  *DoubleP;
   int      *BoolP;
   char **StringP;
} ptr_t;

typedef struct _config {
   char  *LongName;     /* Long Name */
   char  *Doc;       /* Documentation string */
   char  *swtch;        /* Switch Name */
   art artype;      /* Argument Type */
   caddr_t  var;        /* Pointer to the variable */
} confit;

typedef struct _InternalConfig {
   char  *LongName;     /* Long Name */
   char  *Doc;       /* Documentation string */
   char  *swtch;        /* Switch Name */
   art artype;      /* Argument Type */
   ptr_t var;        /* Pointer to the variable */
} Confit;

int required(char **variable, Confit *vparam, int n);
void pusage (char *prog, Confit *cp);
int pconf (
             int    argc,    // arg count
             char  *argv[], // arguments
             Confit *confip   // argument configuratiuon
             );
int ppconf (
              int    argc,    // arg count
              char  *argv[], // arguments
              Confit *confip,  // argument configuratiuon
              char   last    // flag indicating that all remain arguments should be parsed
              );


#endif _PCONF_
/* pconf.h */


