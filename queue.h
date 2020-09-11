#ifndef UPGRADECLUSTER_QUEUE_H
#define UPGRADECLUSTER_QUEUE_H
#include "g.h"


/*an item in queue*/
typedef struct item{
    group *data;
    struct item *next;
} item;

/*the queue for P and O from alg 3*/
typedef struct queue{
    item *head;
    item *last;
    int size;
}queue;


queue* allocate_queue();
void insert_queue(queue *Q, group *G);
group* pop_queue(queue *Q);
int isEmpty(queue *Q);
void free_queue(queue *Q);

#endif /*UPGRADECLUSTER_QUEUE_H*/
