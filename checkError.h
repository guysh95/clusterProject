/*
* Created by guy on 02/09/2020.
*/

#ifndef CHECKP_CHECKERROR_H
#define CHECKP_CHECKERROR_H

void checkMalloc(void* ptr);
void checkFileAccess(int returnVal, int val);
void checkFilePtr(void *filePtr);

#endif /*CHECKP_CHECKERROR_H*/
