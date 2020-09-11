#include "g.h"
#include<stdlib.h>
#include<stdio.h>
#include"checkError.h"

group* allocate_first(int size);
void split(group *g, int *S, group *gPos, group *gNeg);
void writeGroupToFile(group *G, FILE *output);
void free_group(group *g);

/*
 * Allocate the initial group, containing all of the network indices,
 * e.g. size == B->dim.
 */
group* allocate_first(int size){
    int i, *members;
    group *g = malloc(sizeof(group));
    checkMalloc(g);
    g->members = malloc(sizeof(int) * size);
    checkMalloc(g->members);
    g->size = size;
    members = g->members;
    for(i = 0; i < size; i++)
        *members++ = i;
    g->F = calloc(size, sizeof(double));
    checkMalloc(g->F);
    return g;
}

/*
 * This function split the group G, according to the partition S,
 * to two new groups and points gPos and gNeg to those new groups.
 * Assuming gPos and gNeg are pre allocated.
 * After execution G memory is deallocated.
 */
void split(group *G, int *S, group *gPos, group *gNeg){
    int numPos = 0, numNeg = 0;
    int *memPos, *memNeg, *tempPos, *tempNeg, *tempG;
    for(tempG = G->members; tempG < G->members + G->size; tempG++)
        (S[*tempG] == 1) ? numPos++ : numNeg++;
    memPos = malloc(sizeof(int)*numPos);
    checkMalloc(memPos);
    memNeg = malloc(sizeof(int)*numNeg);
    checkMalloc(memNeg);
    tempPos = memPos;
    tempNeg = memNeg;
    for(tempG = G->members; tempG < G->members + G->size; tempG++){
        if(S[*tempG] == 1)
            *tempPos++ = *tempG;
        else
            *tempNeg++ = *tempG;
    }
    gPos->size = numPos;
    gPos->members = memPos;
    gPos->F = malloc(sizeof(double)*numPos);
    checkMalloc(gPos->F);
    gNeg->size = numNeg;
    gNeg->members = memNeg;
    gNeg->F = malloc(sizeof(double)*numNeg);
    checkMalloc(gNeg->F);
    free_group(G);
}

/*
 * this function writes the info about group G
 * to the output file of the program.
 */
void writeGroupToFile(group *G, FILE *output){
    checkFileAccess((int)fwrite(&G->size, sizeof(int), 1, output), 1);
    checkFileAccess((int)fwrite(G->members, sizeof(int), G->size, output), G->size);
}

/*
 * Free the memory allocated to the group g.
 */
void free_group(group *g){
    free(g->members);
    free(g->F);
    free(g);
}
