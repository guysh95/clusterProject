#ifndef UPGRADECLUSTER_QUEUE_H
#define UPGRADECLUSTER_QUEUE_H
#include "g.h"

/*
 * This module defines a queue data structure implemented as a linked list.
 * The queue is used in algorithm 3 to hold the current groups that were divided by the program,
 * and maintain a work flow of FIFO to mend the bigger sized groups first,
 * and also to hold the groups which their partition wouldn't increase the network's modularity.
 */

/*a node in queue*/
typedef struct item{
    group *data; /* the current group this node holds */
    struct item *next; /* a pointer to the next node in the queue */
} item;

/* queue struct definition */
typedef struct queue{
    item *head; /* the first node in the queue */
    item *last; /* the last node in the queue */
    int size; /* the number of nodes in queue */
}queue;

/* allocating an empty queue and returning a pointer to it's location */
queue* allocate_queue();

/* This function accepts a queue Q and a group G and inserts G the end of the queue */
void insert_queue(queue *Q, group *G);

/* This function pops the first group in Q and returns a pointer to that group. */
group* pop_queue(queue *Q);

/* this pointer returns TRUE(== 1) if Q doesn't have any nodes left and FALSE(== 0) otherwise. */
int isEmpty(queue *Q);

/* This function free's all of Q and it's nodes allocated memory. */
void free_queue(queue *Q);

#endif /*UPGRADECLUSTER_QUEUE_H*/
