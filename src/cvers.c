/*
 *    FILE:   cvers.c
 *    AUTHOR: Brian E. Mitchell
 *    DATE:   September 1, 1993
 *
 *      Fortran callable routines to interface with versions.c
 */
 
#include <stdio.h>
#include <string.h>

#include "version.h"
#include "ctof.h"

void concat( char rs[]    , int len, 
              char prefix[], char post[] ,
	      int Lprefix,   int Lpost        )
{
   char *result, *cprefix, *cpost;
   
   cprefix = stringToC(prefix, Lprefix);
   cpost   = stringToC(post,   Lpost  );

   result = (char *) malloc(strlen(cprefix)+strlen(cpost)+2) ;
   result = strcpy(result,cprefix) ;
   result = strcat(result,cpost) ;
   
   stringToF77(result,rs,len);

   free(cprefix);
   free(cpost);
   free(result);
}

void makefilename( char rs[]    , int len, 
                    char path[]  , char prefix[], int *version ,
		    int Lpath    , int Lprefix                  )
{
   char *result, *cpath, *cprefix ;
   
   cpath   = stringToC(path  , Lpath  );
   cprefix = stringToC(prefix, Lprefix);

   result = MakeFileName(cpath,cprefix,*version) ;
   
   stringToF77(result,rs,len);
   free(cpath);  
   free(cprefix);
   free(result);
}


int nextversion(char  path[], char prefix[], 
                  int  Lpath,   int Lprefix   ) 
{
   char *cpath,*cprefix;
   int  val;
   
   cpath   = stringToC(path  , Lpath  );
   cprefix = stringToC(prefix, Lprefix);
   val = nextVersion(cpath,cprefix);
   free(cpath);
   free(cprefix);

   return val;
}

int minversion(char  path[], char prefix[], 
                  int  Lpath,  int Lprefix   ) 
{
   char *cpath,*cprefix;
   int  val;
   
   cpath   = stringToC(path  , Lpath  );
   cprefix = stringToC(prefix, Lprefix);
   val = minVersion(cpath,cprefix);
   free(cpath);
   free(cprefix);

   return val;
}

int maxversion(char  path[], char prefix[], 
                  int  Lpath,   int Lprefix   ) 
{
   char *cpath,*cprefix;
   int  val;
   
   cpath   = stringToC(path  , Lpath  );
   cprefix = stringToC(prefix, Lprefix);
   val = maxVersion(cpath,cprefix);
   free(cpath);
   free(cprefix);

   return val;
}
