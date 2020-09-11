#ifndef CHECKP_CHECKERROR_H
#define CHECKP_CHECKERROR_H
/*
 * this module check for allocations and file access errors,
 * via three defined functions listed below.
 * if an error is discovered the program will abort.
 */

/*
 * to be called after a memory block is allocated to ptr.
 * if ptr == NULL the program will abort.
 */
void checkMalloc(void* ptr);
/*
 * to be called after a file has been accessed for fread\fwrite
 * if returnVal != val then the wrong number of values
 * were read from the file or written to it
 * in this case the program will abort.
 */
void checkFileAccess(int returnVal, int val);
/*
 * to be called after a file has been opened and allocated to filePtr
 * if filePtr == NULL the program will abort
 */
void checkFilePtr(void *filePtr);

#endif /*CHECKP_CHECKERROR_H*/
