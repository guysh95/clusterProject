#include "spmat.h"
#include<stdlib.h>
#include"checkError.h"
#include "g.h"

/*arrays*/
spmat* spmat_allocate(int n, int nnz);
void add_row_array(struct _spmat *A, const int *row, int rowLen, int index);
void free_array(struct _spmat *A);
void mult_array_verProject(const struct _spmat *A, const void *v, double *result, int vType);
int normRowHelp_array(const struct _spmat *A, group *G, int row, int index);
void preNormRow_array(const struct _spmat *A, group *G, int row);
int spmatRowSum_array(const struct _spmat *A, group *G, int index);
double modularityAlterPartA_array(const struct _spmat *A, int *S, int index);


typedef struct {
    /*
     * because this struct is a representation of an adjacency matrix, all of its
     * values are 1 so they are not stored by themselves.
     */
    int *colind; /* this array holds the collumn indices of the matrix values by row to row order */
    int *rowptr; /* this array holds the indices in colind where each row values are stored */
    int *currentCol; /*additional pointer for colind */
    int *currentRow; /*additional pointer for rowptr */
    int count; /* keeping track of how many values were inserted to colind */
    int nnz; /* total number of values to be expected in compMat */
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
    mat->add_row = &add_row_array;
    mat->free = &free_array;
    mat->multSpmat = &mult_array_verProject;
    mat->spmatRowSum = &spmatRowSum_array;
    mat->spmat_modularityAlter = &modularityAlterPartA_array;
    return mat;
}
/*
 * add_row_array adds row to the smpat. note that spmat is an adjacency matrix,
 * and it is being implemented with the structure compMat,
 * so if vertex j is neighbor to vertex i, in the row i that begins in
 * colind[rowptr[i]] there will an array contatining all of i's neighbours including j.
 * input: A - the spmat we mean to add a row to;
 * index - the row we wish to add (index is the vertex i that we want to add to the spmat it's neighbors);
 * row - an array containing the neighbors of vertex i (thus, the row we want to add to the spmat);
 * rowLen - the size of row, as in the number of non zero values in index row.
 */

void add_row_array(struct _spmat *A, const int *row, int rowLen, int index){
    const int* p;
    int n = A->n;
    compMat *mat = A->private;
    /*
     * rowptr treatment: current mat->count is the place where the
     * values of current index row begin
     */
    *(mat->currentRow) = mat->count;
    mat->currentRow++;
    /* colind treatment */
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
 * input: A - the matrix; v- the vector;
 * result - pointer for the result vector;
 * vType - v is of type void to allow a multiplication of A by
 * either int or double (if vType == 0 - v is a double array,
 * if vType == 1 - v is an int array).
 */
void mult_array_verProject(const struct _spmat *A, const void *v, double *result, int vType){
    compMat *mat = A->private;
    int nnz = mat->nnz, i, curRowLen, row = 0,
            *columnPointer = mat->colind, *rowPtr = mat->rowptr, *vInt;
    double *vDouble;
    if(vType == 0){ /* V is double */
        vDouble = (double*)v;
        while(*rowPtr != nnz){
            curRowLen = *(rowPtr+1) - *(rowPtr);
            for(i = 0; i < curRowLen; i++)
                result[row] += vDouble[*columnPointer++];
            row++, rowPtr++;
        }
    }
    else{ /* V is int */
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

/*
 * This function sets a pointer to the current colind values needed
 * by the norm computation of B[g] before normRowHelp_array is being called.
 */
void preNormRow_array(const struct _spmat *A, group *G, int row){
    compMat *mat = A->private;
    mat->currentCol = mat->colind + mat->rowptr[G->members[row]];
}

/*
 * Assuming preNormRow_array has already been called.
 * this function accepts an index of a member in G and
 * checks if that member is a neighbour of the vertex with the index row.
 * if it is the return value will be 1, otherwise 0.
 */
int normRowHelp_array(const struct _spmat *A, group *G, int row, int index){
    compMat *mat = A->private;
    int *endOfRow = mat->colind + mat->rowptr[G->members[row] + 1];
    while(mat->currentCol < endOfRow && *mat->currentCol < G->members[index])
        mat->currentCol++;

    if(mat->currentCol == endOfRow)
        return 0; /* checked all of this row values without finding the current G's member */

    if(*mat->currentCol == G->members[index]){
        mat->currentCol++;
        return 1; /* the current G's member is a neighbour of this row vertex */
    }

    /* if(*mat->currentCol > G->members[index]) */
    return 0;
}


/*
 * calculates the change in the modularity matrix B[g] that depends on A values
 * according to vector S, after the multiplication S[index] = S[index]*(-1)
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

