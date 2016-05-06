/************************************************************/
/* Subroutine: Change_Dir                                   */
/*                                                          */
/* Function: Create and change to subdirectory              */
/*           corresponding to current time step             */
/*                                                          */
/************************************************************/

# include <stdio.h>
# include <sys/stat.h>
# include <unistd.h>
# include <string.h>
# include "dirs.h"

# define STRINGSIZE 128 

void make_dir_(int* index, char* vis_path, int len)
{
	char dirname[STRINGSIZE], curr_dir[STRINGSIZE], idxstr[STRINGSIZE];
    int dirlen, dir;

    strncpy(dirname, vis_path, len);
    sprintf(idxstr, "/%d%d%d", *index/100, (*index%100)/10, *index%10);
    strncat(dirname, idxstr, strlen(idxstr));
    getcwd(curr_dir, STRINGSIZE);
    strncat(curr_dir, dirname, strlen(dirname));
    dir = mkdir(dirname, S_IRWXU|S_IRWXG);

	return;
}
