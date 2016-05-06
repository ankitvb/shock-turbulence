/*
 *    version number support
 *    Brian E. Mitchell
 *    August 31, 1993
 *
 *       The following code provides a number of routines
 *   to support the use of version numbers of files.  The general
 *   form of file names under this system is - prefix.n - where
 *   n represents an integer.
 *
 *      The publicly provided routines are:
 *
 *     int nextVersion(char *path, char *prefix)
 *             Returns a version number one larger than any current file.
 *             To specify the current directory, set path = ".".
 *
 *     int minVersion(char *path, char *prefix)
 *             Returns the lowest version number currently existing.
 *
 *     int maxVersion(char *path, char *prefix)
 *             Returns the highest version number currently existing.
 *
 *     char *MakeFileName(char *path, char *prefix, int version)
 *            Returns a character string containing the full file name.  
 *
 *    NOTES:
 *       It is assumed that the operating system is AT&T Version 5 or
 *     compatible.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <sys/types.h>
#include <dirent.h>

#include "version.h"

#define PRIVATE static
#define PUBLIC


/**************************************************************************/
PRIVATE int getVersion(char *fn,char *prefix)
   /*
    *  Return the version number associated with a given file name
    *  associated with the given prefix.  If the version number is 
    *  not legal, -1 is returned.
    */
{
   long val;
   char *p;
   
   fn += strlen(prefix)+1;
   
   if ( strspn(fn,"0123456789") == strlen(fn) ) {
       errno = 0;
       val = strtol(fn,&p,10);
       if (errno || (p==fn) || p[0] != '\0') val = -1;
   } else {
       val = -1;
   }
   
   return val;
}


PRIVATE void searchDir(char *path,char *prefix, int *min,int *max)
   /*
    *  Search the directory specified by path to determine the 
    *  the lowest and highest version numbers of a file.
    */
{
   DIR *dirp;
   struct dirent *dp;
   char *fn;
   int  ver;
   
   *min = 0;
   *max = 0;

   dirp = opendir(path);
   if (dirp == NULL) {
      fprintf(stderr,"(%s) %s: Directory not found\n",__FILE__,path);
      exit (-1);
   }
   
   while ((dp = readdir(dirp)) != NULL) 
   {
       fn = dp->d_name;
       if (      !strncmp(prefix,fn,strlen(prefix)) 
	      && (fn[strlen(prefix)] == '.')      
	      && ((ver = getVersion(fn,prefix)) != -1)       ) 
       {
          if  (ver > *max)               *max = ver;
	  if ((ver < *min) || (*min==0)) *min = ver;
       }
   }
   (void) closedir(dirp);
   return;
}
/**************************************************************************/


/**************************************************************************/
PUBLIC char *MakeFileName(char *path, char *prefix, int version)
   /*
    *  Return a string of the form path/prefix.version
    */
{
    char *s;
    char n[10];
    
    (void) sprintf(n,"%d",version); 
    s = (char *) malloc(strlen(path)+strlen(prefix)+strlen(n)+3);
    (void) sprintf(s,"%s/%s.%d",path,prefix,version);

    return s;
}

PUBLIC int nextVersion(char *path, char *prefix)
{
    int min,max;
    
    searchDir(path,prefix,&min,&max);
    return max+1;
}

PUBLIC int minVersion(char *path, char *prefix)
{
    int min,max;
    
    searchDir(path,prefix,&min,&max);
    return min;
}

PUBLIC int maxVersion(char *path, char *prefix)
{
    int min,max;
    
    searchDir(path,prefix,&min,&max);
    return max;
}
/**************************************************************************/



#ifdef DEBUG
int main (int argc, char *argv[])
{
  int min,max,next;
  
  
  if (argc != 2) {
     fprintf(stderr,"Usage: %s prefix\n",argv[0]);
     exit(-1);
  }
  
  searchDir(".",argv[1],&min,&max);
  printf ("%d %d\n",min,max);
  
  min  = minVersion (".",argv[1]);
  max  = maxVersion (".",argv[1]);
  next = nextVersion(".",argv[1]);
  
  printf ("%s  %s  %s\n", 
      MakeFileName(".",argv[1],min) ,
      MakeFileName(".",argv[1],max) ,
      MakeFileName(".",argv[1],next)   ) ;
      
  return 1;
}

#endif
