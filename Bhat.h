/*
* Created by guy on 09/08/2020.
*/
#ifndef CHECKP_BHAT_H
#define CHECKP_BHAT_H
#include "spmat.h"
#include<stdio.h>
#include"g.h"

typedef struct _bmat{
    spmat *A; /* sparse matrix representing Adjacency matrix */
    int dim; /* number of vertices in graph */
    int *degs; /* array containing the vertices degrees  - k*/
    int M; /* sum of the graph degrees */
    double norm;
} bmat;


void allocate_bmat(spmat *A, int *degs, int M, bmat *B);
void findEigenPair(bmat *B, group *G, double *eigenVec, double *eigenVal, double *randVec);/* get leading eigen vec and val via matrix shifting and power iteration */
double computeModularity(bmat *B, group *G, int *S, double *temp);
double modularityAlteration(bmat *B, group *G, int *s, int index);
void free_bmat(bmat *B);

#endif
