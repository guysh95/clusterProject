/*
* Created by guy on 02/09/2020.
*/
#include<stdlib.h>
#include <stdio.h>
void checkMalloc(void* ptr);
void checkFilePtr(void *filePtr);
void checkFileAccess(int returnVal, int val);

void checkMalloc(void* ptr){
    if(ptr == NULL){
        printf("malloc malfunction\n");
        exit(1);
    }
    return;
}

void checkFilePtr(void *filePtr){
    if(filePtr == NULL){
        printf("could not open input or output file\n");
        exit(1);
    }
}

void checkFileAccess(int returnVal, int val){
    if(returnVal != val){
        printf("return value of read or write is %d instead of %d\n", returnVal, val);
        exit(1);
    }
    return;
}
