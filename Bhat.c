/*
* Created by guy on 09/08/2020.
*/
#include "Bhat.h"
#include<stdlib.h>
#include <stddef.h>
#include <math.h>
#include<time.h>
#include"checkError.h"


void allocate_bmat(spmat *A, int *degs, int M, bmat *B);
double matrixShifting(bmat *B, group *G);
double rowNorm(bmat *B, group *g, int row);
void matrixUnshift(bmat *B);
double rowSum(bmat *B, group *g, int row);
double vecMult(void *A, void *B, int dim, int aType, int bType);
int vecSub(double *b0, double *b1, int dim);
void vecNorm(double *b1, int dim);
void vecCopy (double *b0, double *b1, int dim);
void findEigenPair(bmat *B, group *G, double *eigenVec, double *eigenVal, double *randVec);
void mult(bmat *B, group *G, const double* v, double* result);
double computeModularity(bmat *B, group *G, int *S, double *temp);
void free_bmat(bmat *B);
double modularityAlteration(bmat *B, group *G, int *s, int index);


void allocate_bmat(spmat *A, int *degs, int M, bmat *B){
    B->degs = degs;
    B->A = A;
    B->dim = A->n;
    B->M = M;
    B->norm = 0.0;
}

void free_bmat(bmat *B){
    B->A->free(B->A);
    free(B->degs);
    free(B);
}


double vecMult(void *A, void *B, int dim, int aType, int bType){
    int i, *bInt, *aInt;
    double sum = 0, *bDouble, *aDouble;
    if(aType == 1){
        bInt = (int *)B;
        aInt = (int *)A;
        for(i = 0; i < dim; i++)
            sum += (*(aInt++)) * (*(bInt++));
    }
    else{
        aDouble = (double *)A;
        if(bType == 1){
            bInt = (int *)B;
            for(i = 0; i < dim; i++)
                sum += (*(aDouble++)) * (*(bInt++));
        }
        else{
            bDouble = (double *)B;
            for(i = 0; i < dim; i++)
                sum += (*(aDouble++)) * (*(bDouble++));
        }
    }
    return sum;
}

/*
 * compute F vector (sums of B[g] rows)
 * assuming g is not the original G = V
 */
void computeF(bmat *B, group *g){
    double *ptr = g->F;
    int i;
    for(i = 0; i < g->size; i++)
        *ptr++ = rowSum(B, g, i);
}

/*
 * returns regular sum of a row in B, for F vector
 * assuming row is according to B[g]
 */
double rowSum(bmat *B, group *g, int row){
    int *k, *kIndex, k_row, size = g->size;
    double sumA = 0, sumB = 0;
    k = B->degs;
    k_row = k[g->members[row]];
    sumA += B->A->spmatRowSum(B->A, g, row);
    for (kIndex = g->members; kIndex < g->members + size; kIndex++)
        sumB += (k_row)*(k[*kIndex]);
    sumB = (double)sumB / B->M;
    return sumA - sumB;
}

/*
 * computes row Norm for matrix shifting using spmat's array implementation
 */
double rowNorm(bmat *B, group *G, int row){
    int i, dim = G->size;
    int *currentK, k_row = B->degs[G->members[row]];
    double sum = 0.0;
    currentK = B->degs;
    preNormRow_array(B->A, G, row);
    for(i = 0; i < dim; i++){
        if(i != row)
            sum += fabs(normRowHelp_array(B->A, G, row, i)
                        - ((double)k_row*(currentK[G->members[i]]))/B->M);
        else
            sum += fabs(- ((double)k_row*(currentK[G->members[i]]))/B->M - G->F[i]);
    }
    return sum;
}
/* result has to be set to zeroes */
void mult(bmat *B, group *G, const double* v, double* result){
    int i, j, *K, n = B->dim,
    *members = G->members, *endMem = members + G->size;
    spmat *A;
    double kPart, *F = G->F;
    A = B->A;
    K = B->degs;
    A->multSpmat(A, v, result, 0);
    for(i = 0; i < n; i++){
        if(members == endMem || i != *members){
            result[i] = 0;
            continue;
        }
        members++;
        kPart = 0;
        for(j = 0; j < n; j++)
            kPart += K[i]*K[j]*v[j];
        kPart = ((double)kPart) / B->M;
        result[i] += -kPart + (B->norm -*(F++))*v[i];
    }
}

double matrixShifting(bmat *B, group *G){
    double max = -HUGE_VAL, norm;
    int n = G->size, i;

    for (i = 0; i < n; i++){
        norm = rowNorm(B, G, i);
        if (norm > max)
            max = norm;
    }
    B->norm = max;
    return max;
}

void matrixUnshift(bmat *B){
    B->norm = 0.0;
}

/* functions for power iteration */

/* returns 1 only if all the ints in the first dim places of both vectors are different. else - returns 0 */
int vecSub(double *b0, double *b1, int dim){
    double sigma = 0.0001, *ptr, *goal = b0 + dim;

    for (ptr = b0; ptr < goal; ptr++, b1++){
        if (fabs(*b1 - *ptr) > sigma)
            return 0;
    }
    return 1;

}
/* normalizes a vector (divides it by it's norm) */
void vecNorm(double *b1, int dim){
    double* ptr, *goal = b1 + dim;
    double nrm = 0;
    for (ptr = b1; ptr < goal; ptr++)
        nrm +=  *ptr * *ptr;
    nrm = sqrt(nrm);
    for (ptr = b1; ptr < goal; ptr++)
        *ptr = *ptr/nrm;
}

/* copies vector b0 to b1 and sets b1 to be a vector of zeroes */
void vecCopy (double *b0, double *b1, int dim){
    double *goal = b0 + dim;
    for ( ; b0 < goal; b0++, b1++){
        *b0 = *b1;
        *b1 = 0.0;
    }
}

/*
 * power iteration and calculating eigen value
 * same as hw2 but eigenVec = b1 (assuming it's already been initialized to all zeroes)
 * and rand = b0, is a random vector we create
 */
void findEigenPair(bmat *B, group *G, double *eigenVec, double *eigenVal, double *randVec){
    double *randPointer;
    double numerator, denominator;
    int isClose, dim, *members = G->members, count = 0, i;
    computeF(B, G);
    dim = B->dim;
    srand(time(NULL));
    randPointer = randVec;
    for(i = 0; i < dim; i++){
        if(count != G->size && i == *members){
            *randPointer = rand();
            members++, count++;
        }
        else
            *randPointer = 0.0;
        randPointer++;
    }

    /*
     * matrix shifting before power iteration in order to
     * find leading eigen vector instead of dominant eigen vector
     */
    matrixShifting(B, G);
    /* power iteration begins here */
    isClose = 0;
    while(isClose == 0){
        mult(B, G, randVec, eigenVec);
        vecNorm(eigenVec, dim);
        isClose = vecSub(randVec, eigenVec, dim);
        if(isClose == 0)
            vecCopy(randVec, eigenVec, dim);
    }
    matrixUnshift(B);
    /*
     * calculating eigen value, using randVec to hold the result of Bhat x eigenVec
     * hence reseting randVec to contain zeroes
     */
    for(randPointer = randVec; randPointer < randVec + dim; randPointer++)
        *randPointer = 0.0;
    mult(B, G, eigenVec, randVec);
    numerator = vecMult(eigenVec, randVec, dim, 0, 0);
    denominator = vecMult(eigenVec, eigenVec, dim, 0, 0);
    *eigenVal = numerator / denominator;
}
void build_temp(double* temp, int dim){
    double *goal = temp + dim;
    for( ; temp < goal; temp++)
        *temp = 0.0;
}

/* assuming S is vector containing only +1/-1 */
double computeModularity(bmat *B, group *G, int *S, double *temp){
    double partA, partB = 0, partF = 0, fixed,
    *F = G->F, *brute, *goal;
    spmat *A = B->A;
    int *K, dim = B->dim, *goalS = S + dim;
    build_temp(temp, dim);
    A->multSpmat(A, S, temp, 1);
    partA = vecMult(temp, S, dim, 0, 1);
    K = B->degs;
    fixed = vecMult(K, S, dim, 1, 1);
    goal = F + G->size;
    for(brute = F; brute < goal; brute++)
        partF += *brute;
    for( ; S < goalS; S++, K++)
        partB += (*S) * (*K) * fixed;
    partB /= B->M;
    return partA - partB - partF;
}

double modularityAlteration(bmat *B, group *G, int *s, int index){
    int *K = B->degs, *ptr, gIndex = G->members[index];
    double partA, partB = 0, alteration;

    partA = B->A->spmat_modularityAlter(B->A, s, gIndex);
    for(ptr = G->members; ptr < G->members + G->size; ptr++){
        if(*ptr != gIndex)
            partB += K[*ptr]*s[*ptr];
    }
    partB = (((double)K[gIndex]) / B->M) * partB;
    alteration = 4 * s[gIndex] * (partA - partB);
    return alteration;
}
