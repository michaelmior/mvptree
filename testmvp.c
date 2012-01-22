/*

    MVPTree c library 
    Copyright (C) 2008-2009 by D. Grant Starkweather
    All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    D Grant Starkweather - starkd88@gmail.com

*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "mvptree.h"

#define MVP_BRANCHFACTOR 2
#define MVP_PATHLENGTH   5
#define MVP_LEAFCAP     25

void testmvpfree(void *ptr){
    free(ptr);
}

static unsigned long long nbcalcs = 0;

int next_poisson_int(int lambda){
/* generate poisson random deviate with average lambda */
    double L = exp((double)-lambda);
    double u, p = 1.0;
    int k = 0;
    do {
	k++;
	u = (float)rand()/((float)RAND_MAX-1);
	p = p*u;
    } while (p > L);
    return k-1;
}


float point_l1_distance(MVPDP *pointA, MVPDP *pointB){
    /* L1 distance */
    if (!pointA || !pointB) return -1.0f;
    unsigned int i;
    unsigned int sum = 0;
    uint8_t *data1 = (uint8_t*)pointA->data;
    uint8_t *data2 = (uint8_t*)pointB->data;
    for (i = 0; i < pointA->datalen; i++){
	int d1 = (int)data1[i];
	int d2 = (int)data2[i];
	int diff = abs(d1 - d2);
	sum += diff;
    }
    nbcalcs++;

    return (float)sum/(float)pointA->datalen;
}
float point_l2_distance(MVPDP *pointA, MVPDP *pointB){
    /* L2 distance */
    if (!pointA || !pointB) {
	return -1.0f;
    } else if (pointA->datalen != pointB->datalen){
	fprintf(stdout," dataA %d, dataB %d\n", pointA->datalen, pointB->datalen);
	return -2.0f;
    }
    int i;
    int sum = 0;
    uint8_t *data1 = (uint8_t*)pointA->data;
    uint8_t *data2 = (uint8_t*)pointB->data;
    for (i=0;i<pointA->datalen;i++){
	int d1 = (int)data1[i];
	int d2 = (int)data2[i];
	int diff = d1 - d2;
	sum += diff*diff;
    }
    nbcalcs++;

    return sqrt((float)sum)/(float)pointA->datalen;
}

MVPDP* mvpdp_dup(MVPDP *dp){
    MVPDP *newpnt = dp_alloc(dp->type);
    if (!newpnt) return NULL;
    newpnt->datalen = dp->datalen;
    newpnt->data = malloc(dp->datalen*dp->type);
    memcpy(newpnt->data, dp->data, dp->datalen*dp->type);
    newpnt->id = strdup(dp->id);
    return newpnt;
}

MVPDP* generate_point(const unsigned int dp_length){
    /* generate one datapoint with BYTEARRAY type data of length dp_length */
    static unsigned long long uid = 0;
    char scratch[32];

    MVPDP *newpnt = dp_alloc(BYTEARRAY);
    if (newpnt == NULL) return NULL;

    newpnt->datalen = dp_length;
    newpnt->data = malloc(dp_length*sizeof(uint8_t));
    if (newpnt->data == NULL) {
	free (newpnt);
	return NULL;
    }
    int8_t *row = (int8_t*)newpnt->data;

    int i;
    for (i=0;i<dp_length;i++){
	row[i] = (int8_t)rand();
    }
    snprintf(scratch, 32, "point%llu", ++uid);
    newpnt->id = strdup(scratch);
    if (newpnt->id == NULL){
	free(newpnt);
	free(newpnt->data);
	return NULL;
    }

    return newpnt;
}

/* generate random data points for a tree */
MVPDP** generate_uniform_points(const unsigned int nbpoints, const unsigned int dp_length){
    /* generate nbpoints uniformly distributed points all of dp_length BYTEARRAY */
    MVPDP **pointlist = (MVPDP**)malloc(nbpoints*sizeof(MVPDP*));
    if (pointlist == NULL){
	return NULL;
    }
    unsigned int i;
    for (i=0;i < nbpoints;i++){
	pointlist[i] = generate_point(dp_length);
	if (!pointlist[i]) return NULL;
    }

    return pointlist;
}

MVPDP** generate_cluster(unsigned int nbpoints, unsigned int dplength, int var, CmpFunc dist){
    /* generate nbpoints cluster with data of type BYTEARRAY and dplength long */

    MVPDP **points = (MVPDP**)malloc(nbpoints*sizeof(MVPDP*));
    if (!points)return NULL;

    /* generate first point around which all other points shall cluster */
    points[0] = generate_point(dplength);
    if (points[0] == NULL){
	free(points);
	return NULL;
    }

    /* generate the other points */
    int i, j, val;
    uint8_t *orig_row = (uint8_t*)points[0]->data;
    float max_distance = 0.0f, d;
    for (i=1;i<nbpoints;i++){
	points[i] = generate_point(dplength);
	if (!points[i]){
	    free(points);
	    return NULL;
	}
	
	uint8_t *new_row  = (uint8_t*)points[i]->data;
	for (j=0;j<dplength;j++){
	    int diff = next_poisson_int(var);
	    float toss = (float)rand()/(float)RAND_MAX-1;
	    val = (toss > 0.5) ? diff : -diff;
	    val += (int)orig_row[j];
	    if (val > 255) val = 255;
	    else if (val < 0) val = 0;
	    new_row[j] = (uint8_t)val;
	}
    
	d = dist(points[0], points[i]);
	if (d > max_distance) {
	    max_distance = d;
	}
    }

    fprintf(stdout,"cluster - maximum distance from main point, %f\n", max_distance);
    return points;
}

int main(int argc, char **argv){

    const unsigned int nbpoints = 100;
    const float radius = 5.0f;
    const unsigned int nbcluster1 = nbpoints/10;
    const unsigned int knearest = nbpoints;
    const unsigned int dplength = 10;
    const int var = 10;
    const char *testfile = "testfile.mvp";
    CmpFunc distance_func = point_l2_distance;
    srand(98293928);

    /* generate points */
    MVPDP  **pointlist = generate_uniform_points(nbpoints, dplength);
    assert(pointlist);

    /* generate cluster */
    MVPDP **cluster1 = generate_cluster(nbcluster1, dplength, var, distance_func);
    assert(cluster1);

    /* save center of cluster for later retrieval */
    MVPDP *cluster_center = mvpdp_dup(cluster1[0]);

    /* create tree */
    MVPTree *tree = mvptree_alloc(NULL, distance_func, MVP_BRANCHFACTOR, MVP_PATHLENGTH, MVP_LEAFCAP);
    assert(tree);

    printf("Add list of %d points to tree.\n", nbpoints);
    /* add points to tree */
    MVPError err = mvptree_add(tree, pointlist, nbpoints);
    assert(err == MVP_SUCCESS);


    printf("Add cluster of %d points to tree.\n", nbcluster1);
    /* add cluster to tree */
    err = mvptree_add(tree, cluster1, nbcluster1);
    assert(err == MVP_SUCCESS);
 
    /* free the point lists but not the points themselves */
    /* the tree takes ownership of the points!!! */
    free(pointlist);
    free(cluster1);


    printf("Write tree out to file, %s.\n", testfile);
    /* write out the tree */
    err = mvptree_write(tree, testfile, 00755);
    assert(err == MVP_SUCCESS);


    mvptree_clear(tree, free);
    free(tree);


    printf("Read tree from %s.\n", testfile);
    /* read it back in again */
    tree = mvptree_read(testfile, distance_func, MVP_BRANCHFACTOR, MVP_PATHLENGTH, MVP_LEAFCAP, &err);
    assert(tree);


    /* add some points to a tree one at a time */
    MVPDP *savedpoint = NULL;
    int count = 0, total = nbpoints/10;
    do {
	MVPDP *pnt = generate_point(dplength);
	if (count == 0) savedpoint = mvpdp_dup(pnt);

	printf("Add point, %s to tree.\n", pnt->id);
	err = mvptree_add(tree, &pnt, 1);
	assert(err == MVP_SUCCESS);
    } while (++count < total);


    printf("Retrieve point, %s.\n", cluster_center->id);

     /* retrieve the cluster center */
    nbcalcs = 0;
    unsigned int nbresults;
    MVPDP **results = mvptree_retrieve(tree, cluster_center , knearest, radius, &nbresults, &err);
    assert(results);
    assert(err == MVP_SUCCESS);

    unsigned int i;
    for (i = 0;i < nbresults;i++){
	fprintf(stdout,"  FOUND --> (%d) %s\n", i, results[i]->id);
    }
    free(results);


    printf("\nRetrieve point, %s.\n", savedpoint->id);

    /* retrieve the first saved point */
    nbcalcs = 0;
    results = mvptree_retrieve(tree, savedpoint, knearest, radius, &nbresults, &err);
    assert(results);
    assert(err == MVP_SUCCESS);

    for (i=0;i<nbresults;i++){
	fprintf(stdout,"  FOUND --> (%d) %s\n", i, results[i]->id);
    }

    free(results);
    mvptree_clear(tree, free);
    free(tree);

    dp_free(cluster_center, free);
    dp_free(savedpoint, free);

    printf("Done.\n");
    return 0;
}
