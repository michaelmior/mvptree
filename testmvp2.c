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
    double sum = 0;
    int32_t *data1 = (int32_t*)pointA->data;
    int32_t *data2 = (int32_t*)pointB->data;
    for (i = 0; i < pointA->datalen; i++){
	double diff = (double)data1[i] - (double)data2[i];
	sum = sum + abs(diff);
    }
    nbcalcs++;
    return (float)sum/(float)pointA->datalen;
}

float point_l2_distance(MVPDP *pointA, MVPDP *pointB){
    /* L2 distance */
    if (!pointA || !pointB || pointA->datalen != pointB->datalen) return -1.0f;
    int i;
    double sum = 0.0;
    int32_t *data1 = (int32_t*)pointA->data;
    int32_t *data2 = (int32_t*)pointB->data;
    for (i=0;i<pointA->datalen;i++){
	double diff =  (double)data1[i] - (double)data2[i];
	sum = sum + pow(diff,2.0f) ;
    }
    nbcalcs++;
    float result = (float)(sqrt(sum)/(float)pointA->datalen);
    return result;
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
    /* generate one datapoint with UINT32ARRAY type data of length dp_length */
    static unsigned long long uid = 0;
    char scratch[32];

    MVPDP *newpnt = dp_alloc(UINT32ARRAY);
    if (newpnt == NULL) return NULL;

    newpnt->datalen = dp_length;
    newpnt->data = malloc(dp_length*sizeof(int32_t));
    if (newpnt->data == NULL) {
	free (newpnt);
	return NULL;
    }
    int32_t *row = (int32_t*)newpnt->data;

    int i;
    for (i=0;i<dp_length;i++){
	row[i] = rand()%1024;
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
    /* generate nbpoints uniformly distributed points all of dp_length UINT32ARRA */
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
    /* generate nbpoints cluster with data of type UINT32ARRAY and dplength long */

    MVPDP **points = (MVPDP**)malloc(nbpoints*sizeof(MVPDP*));
    if (!points)return NULL;

    /* generate first point around which all other points shall cluster */
    points[0] = generate_point(dplength);
    if (points[0] == NULL){
	free(points);
	return NULL;
    }

    /* generate the other points */
    int i, j;
    int32_t *orig_row = (int32_t*)points[0]->data;
    float max_distance = 0.0f, d;
    for (i=1;i<nbpoints;i++){
	points[i] = generate_point(dplength);
	if (!points[i]){
	    free(points);
	    return NULL;
	}
	
	int32_t *new_row  = (int32_t*)points[i]->data;
	for (j=0;j<dplength;j++){
	    int diff = next_poisson_int(var);
	    float toss = (float)rand()/((float)RAND_MAX-1);
	    new_row[j] = (toss > 0.5) ? orig_row[j] + diff : orig_row[j] - diff;
	}

	d = dist(points[0], points[i]);
	if (d > max_distance) {
	    max_distance = d;
	}
    }

    return points;
}

int main(int argc, char **argv){
    if (argc < 3){
	fprintf(stdout,"not enough input args\n");
	fprintf(stdout,"prog <nbpoints> <radius>\n");
	return 0;
    }

    const unsigned int nbpoints = atoi(argv[1]);
    const float radius = atof(argv[2]);
    const unsigned int nbcluster1 = nbpoints/10;
    const unsigned int knearest = nbpoints/2;
    const unsigned int dplength = 10;
    const int var = 10;
    const char *testfile = "testfile.mvp";
    CmpFunc distance_func = point_l2_distance;

    srand(98293928);

    /* generate points */
    MVPDP  **pointlist = generate_uniform_points(nbpoints, dplength);
    assert(pointlist);


    MVPDP **cluster1 = generate_cluster(nbcluster1, dplength, var, distance_func);
    assert(cluster1);


    MVPTree *tree = mvptree_alloc(NULL, distance_func,  MVP_BRANCHFACTOR, MVP_PATHLENGTH, MVP_LEAFCAP);
    assert(tree);

    printf("Add %d points to tree.\n", nbpoints);

    /* add points to tree */
    MVPError err = mvptree_add(tree, pointlist, nbpoints);
    assert(err == MVP_SUCCESS);

    printf("Add cluster of %d points to tree.\n", nbcluster1);

    /* add cluster to tree */
    err = mvptree_add(tree, cluster1, nbcluster1);
    assert(err == MVP_SUCCESS);
 

    /* save the cluster center for later retrieval */
    MVPDP *cluster_center = mvpdp_dup(cluster1[0]);

    /* free the list of points but not the points themselves */
    /* the tree takes ownership of the points  */
    free(pointlist);
    free(cluster1);


    printf("Write the tree to file, %s.\n", testfile);

    /* write the tree out */
    err = mvptree_write(tree, testfile, 00755);
    assert(err == MVP_SUCCESS);

    /* destroy the tree */
    mvptree_clear(tree, free);
    free(tree);


    printf("Read the tree from file, %s.\n", testfile);

    /* read tree back in */
    tree = mvptree_read(testfile, distance_func, MVP_BRANCHFACTOR, MVP_PATHLENGTH, MVP_LEAFCAP, &err);
    assert(tree);

    MVPDP *savedpoint = NULL;
    int count = 0, total = 10;
    do {
	MVPDP *pnt = generate_point(dplength);
	assert(pnt);
	if (count==0) savedpoint = mvpdp_dup(pnt);
	printf("Add point %s.\n", pnt->id);
	err = mvptree_add(tree, &pnt, 1);
	assert(err == MVP_SUCCESS);

    } while (++count < total);


    printf("Retrieve cluster center point, %s.\n", cluster_center->id);

    /* retrieve cluster */ 
    nbcalcs = 0;
    unsigned int nbresults = 0;
    MVPDP **results = mvptree_retrieve(tree, cluster_center, knearest, radius, &nbresults, &err);
    assert(results);
    assert(err == MVP_SUCCESS);

    /* print out results */
    unsigned int i;
    for (i = 0;i < nbresults;i++){
	printf("  FOUND --> (%d) %s\n", i, results[i]->id);
    }

    free(results);

    printf("Retrieved point, %s.\n", savedpoint->id);
    nbcalcs = 0;
    results = mvptree_retrieve(tree, savedpoint, knearest, radius, &nbresults, &err);
    assert(results);
    assert(err == MVP_SUCCESS);
    
    for (i=0;i<nbresults;i++){
	printf("  FOUND --> (%d) %s\n", i, results[i]->id);
    }

    free(results);

    err = mvptree_write(tree, testfile, 00755);
    assert(err == MVP_SUCCESS);


    /* cleanup */
    dp_free(cluster_center, free);
    dp_free(savedpoint, free);
    mvptree_clear(tree, free);
    free(tree);
    printf("Done.\n");
    return 0;
}
