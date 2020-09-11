#ifndef CHECKP_BHAT_H
#define CHECKP_BHAT_H
#include "spmat.h"
#include<stdio.h>
#include"g.h"
/*
 *  bmat is a representation of the modularity matrix of the given network.
 *  We hold the adjacency matrix in "A", a sparse matrix data structure,
 *  the network's vertices degrees are kept in a vector "degs" and their sum is held in "M".
 *  In this representation we compute the needed values in place for each modularity computation,
 *  while these computation are depended on the current subgroup of vertices, those of which we
 *  want to check for their ideal partition.
 */
typedef struct _bmat{
    spmat *A; /* sparse matrix representing the adjacency matrix */
    int dim; /* number of vertices in graph */
    int *degs; /* array containing the vertices degrees */
    int M; /* sum of the vertices degrees */
    double norm; /* used when shifting the matrix to hold the current B[g] norm */
} bmat;

/*
 * the function accepts pre allocated A, degs and computed M,
 * and assign them to a pre allocated pointer B
 */
void allocate_bmat(spmat *A, int *degs, int M, bmat *B);

/* get leading eigen vec and val of B[g] via matrix shifting and power iteration */
void findEigenPair(bmat *B, group *G, double *eigenVec, double *eigenVal, double *randVec);

/* computes modularity of B[g] given a partition S */
double computeModularity(bmat *B, group *G, int *S, double *temp);

/*
 *  compute the alteration to the modularity of B[g] and a partition S
 *  when a single value of S was changed.
 */
double modularityAlteration(bmat *B, group *G, int *s, int index);

/*
 * free B and all the memory it is pointing to.
 */
void free_bmat(bmat *B);

#endif
