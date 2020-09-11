#include<stdlib.h>
#include <stdio.h>
void checkMalloc(void* ptr);
void checkFilePtr(void *filePtr);
void checkFileAccess(int returnVal, int val);

/*
 * to be called after a memory block is allocated to ptr.
 * if ptr == NULL the program will abort.
 */
void checkMalloc(void* ptr){
    if(ptr == NULL){
        printf("malloc malfunction\n");
        exit(1);
    }
}

/*
 * to be called after a file has been opened and allocated to filePtr
 * if filePtr == NULL the program will abort
 */
void checkFilePtr(void *filePtr){
    if(filePtr == NULL){
        printf("could not open input or output file\n");
        exit(1);
    }
}

/*
 * to be called after a file has been accessed for fread\fwrite
 * if returnVal != val then the wrong number of values
 * were read from the file or written to it
 * in this case the program will abort.
 */
void checkFileAccess(int returnVal, int val){
    if(returnVal != val){
        printf("return value of read or write is %d instead of %d\n", returnVal, val);
        exit(1);
    }
}
