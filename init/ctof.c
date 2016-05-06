/*
 *   FILE: ctof.c
 *   AUTHOR: Brian E. Mitchell
 *   DATE:   August 31, 1993
 *
 *   Some routines to aid the conversion of strings 
 *   between FORTRAN and C.  
 *
 *
 *   sgi:
 *     The FORTRAN string format consists of two data elements, an array
 *     of characters padded to the right with spaces, i.e. ' ', and a integer
 *     denoting the length of the array.  There is no null character.
 * 
 *     The provided routines are:
 *          char *stringToC(char p[],int len)
 *          void stringToF77 (char *s, char p[],int len)
 *
 *
 *   cray:
 *     Fortran strings are right padded with spaces and the length of the
 *     string is encoded into the string pointer itself.  The cray provides
 *     some useful utility functions and definitions in the file fortran.h.
 *
 *     The provided routines are:
 *          char *strinToC(_fcd p) 
 *          void stringToF77(char *s, _fcd p)
 *
 */

#include <string.h>
#include "ctof.h"

#define PUBLIC
#define PRIVATE static

PRIVATE void strip(char *s)
   /*
    *    Remove the trailing blanks from the string s.
    */
{
   int i;
   for ( i = strlen(s)-1 ; (i >= 0) && (s[i] == ' ') ; s[i--]='\0') ;
}

PUBLIC char *stringToC(char p[],int len)
   /*
    *   Convert the Fortran string p of length len to C
    *   format.
    */
{
   char *s;
   
   s = (char *) malloc( len*sizeof(char) + 1);
   strncpy(s,p,len);
   s[len] = '\0' ;
   strip(s);
   return s;
}

PUBLIC void stringToF77(char *s, char p[],int len)
   /*
    *   Convert the C format string s to the Fortran format
    *   string p.
    */
{
   int i;
   
   for (i=0         ; i < strlen(s) ; i++) p[i] = s[i] ;
   for (i=strlen(s) ; i < len       ; i++) p[i] = ' '  ;
}
