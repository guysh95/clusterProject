/*
 * spmat.c
 *
 *  Created on: May 13, 2020
 *      Author: Guy Shnaider
 */
#include "spmat.h"
#include<stdlib.h>
#include <stddef.h>
#include"checkError.h"
#include "g.h"

/*arrays*/
spmat* spmat_allocate(int n, int nnz);
void add_row_array_ver_project (struct _spmat *A, const int *row, int rowLen, int index);
void free_array(struct _spmat *A);
void mult_array_verProject(const struct _spmat *A, const void *v, double *result, int vType);
int normRowHelp_array(const struct _spmat *A, group *G, int row, int index);
void preNormRow_array(const struct _spmat *A, group *G, int row);
int spmatRowSum_array(const struct _spmat *A, group *G, int index);
double modularityAlterPartA_array(const struct _spmat *A, int *S, int index);


typedef struct {
    int *colind;
    int *rowptr;
    int *currentCol; /*additional pointer for colind*/
    int *currentRow;
    int count; /* keeping track of how many values were inserted to colind*/
    int nnz; /* counter for number of vals*/
}compMat;


/*
 * spmat_allocate allocates a new empty spmat in the size nXn
 */
spmat* spmat_allocate(int n, int nnz){
    compMat *matImp;
    spmat *mat;
    /*allocating room for the implementation*/
    matImp = (compMat*)malloc(sizeof(compMat));
    checkMalloc(matImp);
    matImp->colind = (int*)malloc(nnz * sizeof(int));
    checkMalloc(matImp->colind);
    matImp->rowptr = (int*)malloc((n+1) * sizeof(int));
    checkMalloc(matImp->colind);
    matImp->currentCol = matImp->colind;
    matImp->currentRow = matImp->rowptr;
    matImp->nnz = nnz;
    matImp->count = 0;
    /*creating spmat to point at implementation */
    mat = (spmat*)malloc(sizeof(spmat));
    checkMalloc(mat);
    mat->n = n;
    mat->private = matImp;
    mat->add_row = &add_row_array_ver_project;
    mat->free = &free_array;
    mat->multSpmat = &mult_array_verProject;
    mat->spmatRowSum = &spmatRowSum_array;
    mat->spmat_modularityAlter = &modularityAlterPartA_array;
    return mat;
}
/*
 * add_row_list adds row to the smpat. note that spmat is an adjacency matrix, and it is being implemented with the structure linked_list,
 * so if vertex i is neighbor to vertex j, in the row i there will be a linked list with an ELEMENT containing (data: 1, column: j).
 * input: A - the spmat we mean to add a row to; i - the row we wish to add (i is the vertex that we want to add to the spmat it's neighbors);
 * inputRow - an array containing the neighbors of vertex i (thus, the row we want to add to the spmat); verk - the size of inputRow (the
 * degree of i).
 */

void add_row_array_ver_project(struct _spmat *A, const int *row, int rowLen, int index){
    const int* p;
    int n = A->n;
    compMat *mat = A->private;
    *(mat->currentRow) = mat->count;
    mat->currentRow++;
    /* values treatment + colind treatment */
    for(p = row ; p < row + rowLen ; p++)
        *(mat->currentCol++) = *p;
    mat->count += rowLen;
    if((index + 1) == n)
        *(mat->currentRow) = mat->count;
}

/*
 * free_array frees the space spmat A takes in the memory once A is no longer in use.
 */
void free_array(struct _spmat *A){
    compMat *mat;
    mat = A->private;
    free(mat->rowptr);
    free(mat->colind);
    free(mat);
    free(A);
}

/*
 * mult_array multiplies vector v by spmat A. the the result vector is result.
 * input: A - the matrix; v- the vector; result - pointer for the result vector; vType - v is of type void to allow a multiplication of A by
 * either int or double (if vType == 0 - v is a double array, if vType == 1 - v is an int array).
 */
void mult_array_verProject(const struct _spmat *A, const void *v, double *result, int vType){
    compMat *mat = A->private;
    int nnz = mat->nnz, i, curRowLen, row = 0,
            *columnPointer = mat->colind, *rowPtr = mat->rowptr, *vInt;
    double *vDouble;
    if(vType == 0){
        vDouble = (double*)v;
        while(*rowPtr != nnz){
            curRowLen = *(rowPtr+1) - *(rowPtr);
            for(i = 0; i < curRowLen; i++)
                result[row] += vDouble[*columnPointer++];
            row++, rowPtr++;
        }
    }
    else{
        vInt = (int*)v;
        while(*rowPtr != nnz){
            curRowLen = *(rowPtr+1) - *(rowPtr);
            for(i = 0; i < curRowLen; i++)
                result[row] += vInt[*columnPointer++];
            row++, rowPtr++;
        }
    }
}

/*
 * spmatRowSum_array calculates the sum of the row index in spmat A
 */
int spmatRowSum_array(const struct _spmat *A, group *G, int index){
    compMat *mat;
    int sum = 0, *members = G->members, *col, *goal, *endMem, rowLen;
    mat = A->private;
    rowLen = mat->rowptr[G->members[index] + 1] - mat->rowptr[G->members[index]];
    col = mat->colind + mat->rowptr[G->members[index]];
    if(rowLen == 0)
        return 0;
    goal = col + rowLen;
    endMem = members + G->size;
    while(members < endMem && col < goal){
        if(*members == *col){
            sum += 1;
            members++, col++;
        }
        else{
            if(*members > *col)
                col++;
            else
                members++;
        }
    }
    return sum;
}

void preNormRow_array(const struct _spmat *A, group *G, int row){
    compMat *mat = A->private;
    mat->currentCol = mat->colind + mat->rowptr[G->members[row]];
}

int normRowHelp_array(const struct _spmat *A, group *G, int row, int index){
    compMat *mat = A->private;
    int *endOfRow = mat->colind + mat->rowptr[G->members[row] + 1];
    while(mat->currentCol < endOfRow && *mat->currentCol < G->members[index])
        mat->currentCol++;
    if(mat->currentCol == endOfRow)
        return 0;
    if(*mat->currentCol == G->members[index]){
        mat->currentCol++;
        return 1;
    }
    /*if(*mat->currentCol > G->members[index])*/
    return 0;
}


/*
 * calculates the modularity of spmat A according to vector S, after the muliplication S[index] = S[index]*(-1)
 */
double modularityAlterPartA_array(const struct _spmat *A, int *S, int index){
    compMat *mat = A->private;
    int *col = mat->colind + mat->rowptr[index], *ptr;
    int rowLen = mat->rowptr[index+1] - mat->rowptr[index];
    int sum = 0;
    for(ptr = col; ptr < col + rowLen; ptr++)
        sum += S[*ptr];
    return sum;
}

