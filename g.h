#ifndef CHECKP_G_H
#define CHECKP_G_H
#include<stdio.h>

/*
 * This module define "group", a struct that holds a subgroup of
 * the network vertices, represented by their indices.
 * This struct is built for helping the in place computation of B[g].
 */
typedef struct newGroup{
    int *members; /* the group indices in an increasing order */
    double *F; /* a vector holding the sum of each row in B[g] */
    int size; /* a variable holding the group's size */
} group;

/*
 * Allocate the initial group, containing all of the network indices
 * e.g. size == B->dim
 */
group* allocate_first(int size);
/*
 * Free the memory allocated to the group g
 */
void free_group(group *g);
/*
 * This function split the group G, according to the partition S,
 * to two new groups and points gPos and gNeg to those new groups.
 * Assuming gPos and gNeg are pre allocated.
 * After execution G memory is deallocated.
 */
void split(group *G, int *S, group *gPos, group *gNeg);
/*
 * this function writes the info about group G
 * to the output file of the program
 */
void writeGroupToFile(group *G, FILE *output);


#endif /*CHECKP_G_H*/
