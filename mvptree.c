/*
 *  MVPTree c library
 *  Copyright (C) 2008-2009 by D. Grant Starkweather.
 *  All rights reserved.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  D Grant Starkweather - starkd88@gmail.com
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include "mvptree.h"

#define HEADER_SIZE 32
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE

const char *tag = "phashmvp2010";
const int version = 0x01000000;

const char *error_msgs[] = {
    "no error",
    "bad argument",
    "no distance function found",
    "mem alloc error",
    "no leaf node created",
    "no internal node created",
    "no path array alloc'd",
    "could not select vantage points",
    "could not calculate range from an sv1",
    "could not calculate range from an sv2",
    "points too compact",
    "could not sort points",
    "could not open file",
    "could not close file",
    "mmap error",
    "unmap eror",
    "no write",
    "could not extend file",
    "could not remap file",
    "datatypes in conflict",
    "no. retrieved exceeds k",
    "empty tree",
    "distance value either NaN or less than zero",
    "could not open file",
    "unrecognized node" };


const char* mvp_errstr(MVPError err) {
    return error_msgs[(int) err];
}

MVPDP* dp_alloc(MVPDataType type) {
    MVPDP *newdp = (MVPDP*) malloc(sizeof(MVPDP));
    newdp->id = NULL;
    newdp->data = NULL;
    newdp->datalen = 0;
    newdp->type = type;
    newdp->path = NULL;
    return newdp;
}

void dp_free(MVPDP *dp, MVPFreeFunc free_func) {
    if (dp) {
        if (dp->path) { free(dp->path); }
        if (free_func) {
            if (dp->id)    free_func(dp->id);
            if (dp->data) free_func(dp->data);
        }
        free(dp);
    }
}

MVPTree* mvptree_alloc(MVPTree *tree, CmpFunc distance, unsigned int bf, unsigned int p,
                                                           unsigned int k) {
    if (distance == NULL) {
        return NULL;
    }

    MVPTree *retTree;
    if (tree == NULL) {
        retTree = (MVPTree*) malloc(sizeof(MVPTree));
        if (retTree  == NULL) { return NULL; }
    } else {
        retTree = tree;
    }

    retTree->branchfactor = bf;
    retTree->pathlength   = p;
    retTree->leafcap      = k;
    retTree->dist         = distance;
    retTree->datatype     = 0;
    retTree->node         = NULL;
    retTree->fd           = 0;
    retTree->k            = 0;
    retTree->size         = 0;
    retTree->pos          = 0;
    retTree->buf          = NULL;
    retTree->pgsize       = sysconf(_SC_PAGESIZE);

    return retTree;
}

/* custom isnan function */
static int is_nan(float x) {
    float var = x;
    return (var != var) ? 1 : 0;
}

static Node* create_leaf(unsigned int leafcap) {
    Node *node = (Node*) malloc(sizeof(Node));
    node->leaf.sv1 = NULL;
    node->leaf.sv2 = NULL;
    node->leaf.points = (MVPDP**) calloc(leafcap, sizeof(MVPDP*));
    node->leaf.d1 = (float*) calloc(leafcap,sizeof(float));
    node->leaf.d2 = (float*) calloc(leafcap,sizeof(float));
    node->leaf.nbpoints = 0;
    node->leaf.type = LEAF_NODE;

    return node;
}

static Node* create_internal(unsigned int bf) {
    Node *node = (Node*) malloc(sizeof(Node));
    node->internal.sv1 = NULL;
    node->internal.sv2 = NULL;
    node->internal.M1 = (float*) calloc(bf-1, sizeof(float));
    node->internal.M2 = (float*) calloc(bf, sizeof(float));
    node->internal.child_nodes = calloc(bf*bf, sizeof(Node*));
    node->internal.type = INTERNAL_NODE;

    return node;
}

static void free_node(Node *node) {
    if (node) {
        if (node->leaf.type == LEAF_NODE) {
            free(node->leaf.points);
            free(node->leaf.d1);
            free(node->leaf.d2);
        } else if (node->internal.type == INTERNAL_NODE) {
            free(node->internal.M1);
            free(node->internal.M2);
            free(node->internal.child_nodes);
        }
        free(node);
    }
}

static void _mvptree_clear(MVPTree *tree, Node *node, MVPFreeFunc free_func, int lvl) {
    if (!node) { return; }

    int i;
    if (node->internal.type == INTERNAL_NODE) {
        int fanout = tree->branchfactor*tree->branchfactor;
        for (i = 0; i < fanout; i++) {
            _mvptree_clear(tree, node->internal.child_nodes[i], free_func, lvl+1);
        }
        dp_free(node->internal.sv1, free_func);
        dp_free(node->internal.sv2, free_func);
    } else {
        dp_free(node->leaf.sv1, free_func);
        dp_free(node->leaf.sv2, free_func);

        for (i=0; i< node->leaf.nbpoints; i++) {
            dp_free(node->leaf.points[i], free_func);
        }
    }
    free_node(node);
}

void mvptree_clear(MVPTree *tree, MVPFreeFunc free_func) {
    if (!tree || !tree->node) { return; }
    _mvptree_clear(tree, tree->node, free_func, 0);
}

/*
   Select the two points at maximum distance from each other using the dist metric.
   Return the positions in list of points in sv1_pos and sv2_pos
*/

static int select_vantage_points(MVPDP **points, unsigned int nb,int *sv1_pos, int *sv2_pos,
                                                      CmpFunc dist) {
    if (!points || !sv1_pos || !sv2_pos || !dist || nb == 0) { return -1; }

    *sv1_pos = (nb >= 1) ? 0 : -1;
    *sv2_pos = -1;

    float max_dist = 0.0f, d;
    int i, j;
    for (i = 0; i < nb; i++) {
        for (j = i+1; j < nb; j++) {
            d = dist(points[i], points[j]);
            if (is_nan(d) || d < 0.0f) {
                return -2;
            }

            if (d > max_dist) {
                max_dist = d;
                *sv1_pos = i;
                *sv2_pos = j;
            }
        }
    }

    return 0;
}


static int find_splits(MVPDP **points, unsigned int nb, MVPDP *vp, MVPTree *tree,
                                                     float *M, unsigned int lengthM) {
    if (!points || nb == 0 || !M || lengthM == 0) return -1;

    CmpFunc distfunc = tree->dist;
    float *dist = (float*) malloc(nb*sizeof(float));

    int i, j;
    for (i = 0; i < nb; i++) {
        dist[i] = distfunc(points[i], vp);
        if (is_nan(dist[i]) || dist[i] < 0.0f) {
            free(dist);
            return -2;
        }
    }

    int min_pos;
    for (i = 0;i < nb-1;i++) {
        min_pos = i;
        for (j = i+1; j < nb; j++) {
            if (dist[j] < dist[min_pos]) {
                min_pos = j;
            }
        }

        if (min_pos != i) {
            float tmp = dist[min_pos];
            dist[min_pos] = dist[i];
            dist[i] = tmp;
        }
    }

    for (i = 0; i < lengthM; i++) {
        int index = (i+1)*nb/(lengthM+1);
        if (index <= 0) { index = 0; }
        if (index >= nb) { index = nb-1; }
        M[i] = dist[index];
    }

    free(dist);
    return 0;
}

/* Sort points into bins by distance(points[i], dp) for each i in list, skipping
 * points[sv1_pos] and points[sv2_pos]. Use pivot[LengthM1] array as pivot points
 * to determine which bins. */
static MVPDP*** sort_points(MVPDP **points, unsigned int nbpoints, int sv1_pos, int sv2_pos,
                         MVPDP *vp, MVPTree *tree, int **counts, float *pivots) {
    if (!points || !vp || !tree || !counts || !pivots || nbpoints == 0) { return NULL; }

    CmpFunc distfunc = tree->dist;
    int bf = tree->branchfactor;
    int lengthM1 = bf-1;

    MVPDP*** bins = (MVPDP***) malloc(bf*sizeof(MVPDP**));
    if (!bins) return NULL;

    *counts = (int*) calloc(bf, sizeof(int));
    if (!counts) {
        free(bins);
        return NULL;
    }

    int i, k;
    for (i=0; i < bf; i++) {
        bins[i] = (MVPDP**) malloc(nbpoints*sizeof(MVPDP*));
        if (!bins[i]) { return NULL; }
    }

    for (i=0;i<nbpoints;i++) {
        if (i == sv1_pos || i == sv2_pos) { continue; }

        float d = distfunc(vp, points[i]);
        if (is_nan(d) || d < 0.0f) {
            free(counts);
            free(bins);
            return NULL;
        }

        for (k = 0; k < lengthM1; k++) {
            if (d <= pivots[k]) {
                bins[k][(*counts)[k]] = points[i];
                (*counts)[k]++;
                break;
            }
        }

        if (d > pivots[lengthM1-1]) {
            bins[lengthM1][(*counts)[lengthM1]] = points[i];
            (*counts)[lengthM1]++;
        }
    }

    return bins;
}

/* Calculate distances for all points from given vantage point, vp, and
   assign that distance into each points path using the lvl parameter. */
static int find_distance_range_for_vp(MVPDP **points, unsigned int nbpoints, MVPDP *vp,
                                      MVPTree *tree, int lvl) {
    if (!points || nbpoints == 0 || !vp || !tree || !tree->dist) {
        return -1;
    }

    CmpFunc func = tree->dist;
    int i, error = 0;
    for (i = 0; i < nbpoints; i++) {
        float d = func(vp, points[i]);
        if (is_nan(d) || d < 0.0f) {
            return -2;
        }
        if (lvl < tree->pathlength) {
            points[i]->path[lvl] = d;
        }
    }

    return error;
}

static Node* _mvptree_add(MVPTree *tree, Node *node, MVPDP **points, unsigned int nbpoints,
                                                               MVPError *error, int lvl) {
    Node *new_node = node;
    if (nbpoints == 0) { return new_node; }
    if (!tree || lvl < 0 || !points) {
        *error = MVP_ARGERR;
        return NULL;
    }

    CmpFunc dist_fnc = tree->dist;
    int bf = tree->branchfactor, lengthM1 = bf-1;

    if (new_node == NULL) { /* create new node */
        int sv1_pos, sv2_pos;
        if (nbpoints <= tree->leafcap + 2) {
            /* create leaf node */
            new_node = create_leaf(tree->leafcap);
            if (!new_node) {
                *error = MVP_NOLEAF;
                return NULL;
            }

            if (select_vantage_points(points, nbpoints, &sv1_pos, &sv2_pos, tree->dist) < 0) {
                *error = MVP_VPNOSELECT;
                free_node(new_node);
                return NULL;
            }

            new_node->leaf.sv1 = (sv1_pos >= 0) ? points[sv1_pos] : NULL;
            new_node->leaf.sv2 = (sv2_pos >= 0) ? points[sv2_pos] : NULL;

            if (find_distance_range_for_vp(points, nbpoints, new_node->leaf.sv1,tree,lvl) < 0) {
                *error = MVP_NOSV1RANGE;
                free_node(new_node);
                return NULL;
            }

            if (new_node->leaf.sv2) {
                if (find_distance_range_for_vp(points,nbpoints,new_node->leaf.sv2,tree,lvl+1) < 0) {
                    *error = MVP_NOSV2RANGE;
                    free_node(new_node);
                    return NULL;
                }
            }

            /* add remaining points to leaf */
            int i, count = 0;
            for (i=0; i < nbpoints; i++) {
                if (i == sv1_pos || i == sv2_pos) { continue; }
                new_node->leaf.d1[count] = dist_fnc(points[i], new_node->leaf.sv1);
                new_node->leaf.d2[count] = dist_fnc(points[i], new_node->leaf.sv2);
                new_node->leaf.points[count++] = points[i];
            }
            new_node->leaf.nbpoints = count;
        } else { /* create internal node */
            new_node = create_internal(tree->branchfactor);
            if (!new_node) {
                *error = MVP_NOINTERNAL;
                return NULL;
            }

            if (select_vantage_points(points, nbpoints, &sv1_pos, &sv2_pos, tree->dist) < 0) {
                *error = MVP_VPNOSELECT;
                free_node(new_node);
                return NULL;
            }

            new_node->internal.sv1 = points[sv1_pos];
            new_node->internal.sv2 = points[sv2_pos];

            if (find_distance_range_for_vp(points,nbpoints,new_node->internal.sv1,tree,lvl) < 0) {
                *error = MVP_NOSV1RANGE;
                free_node(new_node);
                return NULL;
            }

            if (find_splits(points, nbpoints, new_node->internal.sv1, tree,
                    new_node->internal.M1,lengthM1) < 0) {
                *error = MVP_NOSPLITS;
                free_node(new_node);
                return NULL;
            }

            int i, j;
            int *binlengths = NULL;
            MVPDP ***bins = sort_points(points, nbpoints, sv1_pos, sv2_pos,
                                     new_node->internal.sv1, tree, &binlengths, new_node->internal.M1);
            if (!bins) {
                *error = MVP_NOSORT;
                free_node(new_node);
                return NULL;
            }

            for (i=0; i < tree->branchfactor; i++) {
                /* for each bin */
                if (find_distance_range_for_vp(bins[i], binlengths[i], new_node->internal.sv2,
                                   tree, lvl+1) < 0) {
                    *error = MVP_NOSV2RANGE;
                    free_node(new_node);
                    for (j=0; j < tree->branchfactor; j++) { free(bins[j]); }
                    free(bins);
                    return NULL;
                }

                if (find_splits(bins[i], binlengths[i], new_node->internal.sv2, tree,
                        new_node->internal.M2 + i*lengthM1,lengthM1) < 0) {
                    *error = MVP_NOSPLITS;
                    free_node(new_node);
                    for (j=0; j < tree->branchfactor; j++) { free(bins[j]); }
                    free(bins);
                    return NULL;
                }

                int *bin2lengths = NULL;
                MVPDP ***bins2 = sort_points(bins[i], binlengths[i], -1, -1, new_node->internal.sv2,
                                              tree, &bin2lengths, new_node->internal.M2 + i*lengthM1);

                if (!bins2) {
                    *error = MVP_NOSORT;
                    for (j=0; j < tree->branchfactor; j++) { free(bins[j]); }
                    free(bins);
                    free_node(new_node);
                    return NULL;
                }

                int j;
                for (j=0; j < tree->branchfactor; j++) {
                    /* for each row of 2nd tier bins
                     * index into child node = i*branchfactor + j */
                    Node *child = _mvptree_add(tree, NULL, bins2[j], bin2lengths[j], error, lvl+2);
                    new_node->internal.child_nodes[i*tree->branchfactor+j] = child;
                }

                free(bin2lengths);
                for (j = 0; j < tree->branchfactor; j++) { free(bins2[j]); }
                free(bins2);
            }

            free(binlengths);
            for (i=0; i < tree->branchfactor; i++) { free(bins[i]); }
            free(bins);
        }
    } else { /* node already exists */
        if (new_node->leaf.type == LEAF_NODE) {
            if (new_node->leaf.nbpoints + nbpoints <= tree->leafcap) {
                /* add points into leaf - plenty of room */
                if (find_distance_range_for_vp(points,nbpoints,new_node->leaf.sv1,tree,lvl) < 0) {
                    *error = MVP_NOSV1RANGE;
                    return new_node;
                }

                int pos = 0;
                if (new_node->leaf.sv2 == NULL) {
                    new_node->leaf.sv2 = points[0];
                    pos = 1;
                }

                if (find_distance_range_for_vp(points,nbpoints,new_node->leaf.sv2,tree,lvl+1) < 0) {
                    *error = MVP_NOSV2RANGE;
                    return new_node;
                }

                int count = new_node->leaf.nbpoints;
                for (; pos < nbpoints;pos++) {
                    new_node->leaf.d1[count] = tree->dist(points[pos], new_node->leaf.sv1);
                    new_node->leaf.d2[count] = tree->dist(points[pos], new_node->leaf.sv2);
                    new_node->leaf.points[count++] = points[pos];
                }
                new_node->leaf.nbpoints = count;
            } else {
                /* not enough room in current leaf - create new node */
                unsigned int new_nb = new_node->leaf.nbpoints + nbpoints;
                if (new_node->leaf.sv1) { new_nb++; }
                if (new_node->leaf.sv2) { new_nb++; }

                MVPDP **tmp_pts = (MVPDP**) malloc(new_nb*sizeof(MVPDP*));
                if (!tmp_pts) {
                    *error = MVP_MEMALLOC;
                    return new_node;
                }

                int i, index = 0;
                if (new_node->leaf.sv1) { tmp_pts[index++] = new_node->leaf.sv1; }
                if (new_node->leaf.sv2) { tmp_pts[index++] = new_node->leaf.sv2; }
                for (i=0; i < new_node->leaf.nbpoints; i++) {
                    tmp_pts[index++] = new_node->leaf.points[i];
                }
                for (i=0; i < nbpoints; i++) {
                    tmp_pts[index++] = points[i];
                }
                Node *old_node = new_node;
                free_node(old_node);
                new_node = _mvptree_add(tree, NULL, tmp_pts, new_nb, error, lvl);

                free(tmp_pts);
            }
        } else { /* node is internal - must recurse on subnodes */
            if (find_distance_range_for_vp(points, nbpoints, new_node->internal.sv1,tree,lvl) < 0) {
                *error = MVP_NOSV1RANGE;
                return new_node;
            }

            int *binlengths = NULL;
            MVPDP ***bins = sort_points(points, nbpoints, -1, -1, new_node->internal.sv1,
                                                      tree, &binlengths, new_node->internal.M1);
            int i;
            if (!bins) {
                *error = MVP_NOSORT;
                return new_node;
            }

            for (i=0; i < tree->branchfactor; i++) {
                /* for each bin */
                if (binlengths[i] <= 0) {
                    continue;
                }
                int j;
                if (find_distance_range_for_vp(bins[i], binlengths[i],
                                                new_node->internal.sv2, tree, lvl+1) < 0) {
                    *error = MVP_NOSV2RANGE;
                    for (j=0; j < tree->branchfactor; j++) { free(bins[j]); }
                    free(bins);
                    return new_node;
                }

                int *bin2lengths = NULL;
                MVPDP ***bins2 = sort_points(bins[i], binlengths[i], -1, -1,new_node->internal.sv2,\
                           tree, &bin2lengths, new_node->internal.M2 + i*lengthM1);

                if (!bins2) {
                    *error = MVP_NOSORT;
                    for (j=0; j < tree->branchfactor; j++) { free(bins[j]); }
                    free(bins);
                    return new_node;
                }
                for (j=0; j < tree->branchfactor; j++) {
                    /* for each row of 2nd tier bins */
                    /* index into child node = i*branchfactor + j      */
                    Node *child;
                    child = _mvptree_add(tree,
                                          new_node->internal.child_nodes[i*tree->branchfactor + j],
                                          bins2[j],  bin2lengths[j],error, lvl+2);
                    new_node->internal.child_nodes[i*tree->branchfactor+j] = child;
                    if (*error != MVP_SUCCESS) { break; }
                }

                free(bin2lengths);
                for (j=0; j <tree->branchfactor;j++) { free(bins2[j]); }
                free(bins2);
            }

            free(binlengths);
            for (i=0;i<tree->branchfactor;i++) { free(bins[i]); }
            free(bins);
        }
    }
    return new_node;
}

MVPError mvptree_add(MVPTree *tree, MVPDP **points, unsigned int nbpoints) {
    MVPError err = MVP_SUCCESS;
    if (nbpoints == 0) { return err; }

    if (tree && points) {
        if (tree->datatype == 0) {
            tree->datatype = points[0]->type;
        }
        if (tree->datatype != points[0]->type) {
            return MVP_TYPEMISMATCH;
        }

        unsigned int i;
        for (i=0; i < nbpoints; i++) {
            points[i]->path = (float*) malloc(tree->pathlength*sizeof(float));
            if (points[i]->path == NULL) {
                return MVP_PATHALLOC;
            }
            memset(points[i]->path, 0, tree->pathlength*sizeof(float));
        }
        tree->node = _mvptree_add(tree, tree->node, points, nbpoints, &err, 0);
    } else {
        err = MVP_ARGERR;
    }
    return err;
}


static
MVPError _mvptree_retrieve(MVPTree *tree,Node *node,MVPDP *target, float radius, MVPDP** results,
                                                          unsigned int *nbresults, int lvl) {
    MVPError err = MVP_SUCCESS;
    int bf = tree->branchfactor;
    int lengthM1 = bf - 1;
    float d1, d2;
    if (node == NULL) return err;

    CmpFunc distance = tree->dist;
    unsigned int i, j;

    if (node->leaf.type == LEAF_NODE) {
        d1 = distance(target, node->leaf.sv1);
        if (is_nan(d1) || d1 < 0.0f) {
            return MVP_BADDISTVAL;
        }

        if (lvl < tree->pathlength) { target->path[lvl] = d1; }
        if (d1 <= radius) {
            results[(*nbresults)++] = node->leaf.sv1;
            if (*nbresults >= tree->k) { return MVP_KNEARESTCAP; }
        }
        if (tree->node->leaf.sv2) {
            d2 = distance(target, node->leaf.sv2);

            if (is_nan(d2) || d2 < 0.0f) {
                return MVP_BADDISTVAL;
            }
            if (d2 <= radius) {
                results[(*nbresults)++] = node->leaf.sv2;
                if (*nbresults >= tree->k) { return MVP_KNEARESTCAP; }
            }
            if (lvl+1 < tree->pathlength) { target->path[lvl+1] = d2; }

            for (i=0; i < node->leaf.nbpoints; i++) {
                /* check all points
                float d = distance(target,node->leaf.points[i]);
                fprintf(stdout, "pnt%d distance(Q,%s)=%f\n",i,node->leaf.points[i]->id, d);
                if (d <= radius) {
                    results[(*nbresults)++] = node->leaf.points[i];
                    if (*nbresults >= tree->k) {
                        return MVP_KNEARESTCAP;
                    }
                }
                */

                /* filter points before checking */
                if (d1 - radius <= node->leaf.d1[i] && d1 + radius >= node->leaf.d1[i]) {
                    if (d2 - radius <= node->leaf.d2[i] && d2 + radius >= node->leaf.d2[i]) {
                        int endpath = (lvl+1 < tree->pathlength) ? lvl+1 : tree->pathlength;
                        int skip = 0;
                        for (j=0; j < endpath; j++) {
                            if (target->path[j] - radius <= node->leaf.points[i]->path[j] &&
                                target->path[j] + radius >= node->leaf.points[i]->path[j]) {
                                continue;
                            } else {
                                skip = 1;
                                break;
                            }
                        }

                        if (!skip) {
                            float d = distance(target, node->leaf.points[i]);
                            if (is_nan(d) || d < 0.0) {
                                return MVP_BADDISTVAL;
                            }
                            if (d <= radius) {
                                results[(*nbresults)++] = node->leaf.points[i];
                                if (*nbresults >= tree->k) {
                                    return MVP_KNEARESTCAP;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else if (node->internal.type == INTERNAL_NODE) {
        d1 = distance(target, node->internal.sv1);
        if (is_nan(d1) || d1 < 0.0f) {
            return MVP_BADDISTVAL;
        }
        if (d1 <= radius) {
            results[(*nbresults)++] = node->internal.sv1;
            if (*nbresults >= tree->k) return MVP_KNEARESTCAP;
        }
        if (lvl < tree->pathlength) target->path[lvl] = d1;
        d2 = distance(target, node->internal.sv2);
        if (is_nan(d2) || d2 < 0.0f) {
            return MVP_BADDISTVAL;
        }
        if (d2 <= radius) {
            results[(*nbresults)++] = node->internal.sv2;
            if (*nbresults >= tree->k) { return MVP_KNEARESTCAP; }
        }
        if (lvl+1 < tree->pathlength) { target->path[lvl+1] = d2; }

        /* check <= each 1st level bins */
        for (i=0; i < lengthM1; i++) {
            if (d1 - radius <= node->internal.M1[i]) {
                /* check <= each 2nd level bins */
                for (j=0; j < lengthM1; j++) {
                    if (d2 - radius <= node->internal.M2[i*lengthM1+j]) {

                        err = _mvptree_retrieve(tree,node->internal.child_nodes[i*bf+j],target,
                                                                        radius, results, nbresults, lvl+2);

                        if (err != MVP_SUCCESS) { return err; }
                    }
                }
                /* check >= last 2nd level bin  */
                if (d2 + radius >= node->internal.M2[i*lengthM1+lengthM1-1]) {

                    err = _mvptree_retrieve(tree,node->internal.child_nodes[i*bf+lengthM1],
                                                target, radius, results, nbresults, lvl+2);
                    if (err != MVP_SUCCESS) { return err; }
                }
            }
        }

        /* check >= last 1st level bin */
        if (d1 + radius >= node->internal.M1[lengthM1-1]) {
            /* check <= each 2nd level bins */
            for (j=0; j < lengthM1; j++) {
                if (d2 - radius <= node->internal.M2[lengthM1*lengthM1+j]) {
                    err = _mvptree_retrieve(tree,node->internal.child_nodes[bf*lengthM1+j],
                                                    target, radius, results, nbresults, lvl+2);
                    if (err != MVP_SUCCESS) { return err; }
                }
            }
            /* check >= last 2nd level bin  */

            if (d2 + radius >= node->internal.M2[lengthM1*lengthM1+lengthM1-1]) {

                err = _mvptree_retrieve(tree,node->internal.child_nodes[bf*lengthM1+lengthM1],
                                                 target, radius, results, nbresults, lvl+2);
                if (err != MVP_SUCCESS) { return err; }
            }
        }
    } else {
        err = MVP_UNRECOGNIZED;
    }

    return err;
}

MVPDP** mvptree_retrieve(MVPTree *tree, MVPDP *target, unsigned int knearest, float radius,
                                                 unsigned int *nbresults,MVPError *error) {
    if (!tree || !target || !nbresults || knearest == 0 || radius < 0) {
        *error = MVP_ARGERR;
        return NULL;
    }

    if (!tree->dist) {
        *error = MVP_NODISTANCEFUNC;
        return NULL;
    }

    *nbresults = 0;
    *error = MVP_SUCCESS;

    if (!tree->node) {
        *error = MVP_EMPTYTREE;
        return NULL;
    }

    MVPDP **results = (MVPDP**) malloc(knearest*sizeof(MVPDP*));
    if (!results) {
        *error = MVP_MEMALLOC;
        return NULL;
    }

    target->path = (float*) malloc(tree->pathlength*sizeof(float));
    if (target->path == NULL) {
        *error = MVP_MEMALLOC;
        free(results);
        return NULL;
    }
    tree->k = knearest;

    *error = _mvptree_retrieve(tree, tree->node, target, radius, results, nbresults, 0);

    free(target->path);
    target->path = NULL;

    return results;
}

static off_t write_datapoint(MVPDP *dp, MVPTree *tree) {
    off_t start = tree->pos;
    off_t pos = tree->pos;
    uint8_t active = 0;
    uint32_t bytelength = 0;
    char *buf = tree->buf;

    if (dp == NULL) {
        memcpy(&buf[pos++], &active, 1);
        memcpy(&buf[pos], &bytelength, sizeof(uint32_t));
        pos += sizeof(uint32_t);
        tree->pos = pos;
        return start;
    }
    active = 1;
    uint8_t idlen = strlen(dp->id);
    uint32_t datalength = dp->datalen;
    uint8_t type = dp->type;
    bytelength = sizeof(uint8_t) + idlen + sizeof(uint32_t) +
                            datalength*type + (tree->pathlength)*sizeof(float);

    memcpy(&buf[pos++], &active    , 1);
    memcpy(&buf[pos]  , &bytelength, sizeof(uint32_t));
    pos += sizeof(uint32_t);
    memcpy(&buf[pos++], &idlen     , 1);
    memcpy(&buf[pos]  , dp->id     , idlen);
    pos += idlen;
    memcpy(&buf[pos]  , &datalength, sizeof(uint32_t));
    pos += sizeof(uint32_t);
    memcpy(&buf[pos]  , dp->data   , datalength*type);
    pos += datalength*type;
    memcpy(&buf[pos]  , dp->path   , (tree->pathlength)*sizeof(float));
    pos += (tree->pathlength)*sizeof(float);

    tree->pos = pos;
    return start;
}

static int extend_mvpfile(MVPTree *tree) {
    if (munmap(tree->buf, tree->size) < 0) {
        return -1;
    }
    if (ftruncate(tree->fd, tree->size + tree->pgsize) < 0) {
        return -2;
    }

    tree->size += tree->pgsize;
    char *buf = (char*) mmap(NULL, tree->size, PROT_READ|PROT_WRITE, MAP_SHARED, tree->fd, 0);
    if (buf == NULL) {
        return -3;
    }
    tree->buf = buf;

    return 0;
}

static off_t _mvptree_write(MVPTree *tree,Node *node,MVPError *error,int lvl) {
    off_t start_pos = tree->pos;
    if (node == NULL) { return 0; }

    uint8_t node_type = (uint8_t) node->leaf.type;
    if (node->leaf.type == LEAF_NODE) {
        uint32_t nbpoints = node->leaf.nbpoints;
        if (tree->pos >= tree->size - tree->pgsize/2) {
            if (extend_mvpfile(tree) < 0) {
                *error = MVP_FILETRUNCATE;
                return start_pos;
            }
        }

        /* save node */
        memcpy(&tree->buf[tree->pos++], &node_type, 1);
        write_datapoint(node->leaf.sv1, tree);
        write_datapoint(node->leaf.sv2, tree);
        memcpy(&tree->buf[tree->pos], &nbpoints, sizeof(uint32_t));
        tree->pos += sizeof(uint32_t);

        /* write points */
        int i;
        off_t saved_pos = tree->pos;
        tree->pos += (tree->leafcap)*(2*sizeof(float)+sizeof(off_t));
        for (i=0; i < nbpoints; i++) {
            if (tree->pos >= tree->size - tree->pgsize/2) {
                if (extend_mvpfile(tree) < 0) {
                    *error = MVP_FILETRUNCATE;
                    break;
                }
            }

            memcpy(&tree->buf[saved_pos], &(node->leaf.d1[i]), sizeof(float));
            saved_pos += sizeof(float);
            memcpy(&tree->buf[saved_pos], &(node->leaf.d2[i]), sizeof(float));
            saved_pos += sizeof(float);

            off_t offset = write_datapoint(node->leaf.points[i], tree);
            memcpy(&tree->buf[saved_pos], &offset, sizeof(off_t));
            saved_pos += sizeof(off_t);
        }
    } else if (node->internal.type == INTERNAL_NODE) {
        const uint8_t fileno = 0;
        int bf = tree->branchfactor;
        int lengthM1 = bf - 1;
        int lengthM2 = (bf - 1)*bf;
        int fanout   = bf*bf;

        memcpy(&tree->buf[tree->pos++], &node_type, 1);
        write_datapoint(node->internal.sv1, tree);
        write_datapoint(node->internal.sv2, tree);
        memcpy(&tree->buf[tree->pos], node->internal.M1, lengthM1*sizeof(float));
        tree->pos += lengthM1*sizeof(float);
        memcpy(&tree->buf[tree->pos], node->internal.M2, lengthM2*sizeof(float));
        tree->pos += lengthM2*sizeof(float);

        off_t saved_pos = tree->pos;
        tree->pos += fanout*(sizeof(uint8_t) + sizeof(off_t));
        int i;
        for (i=0; i < fanout; i++) {
            if (tree->pos >= tree->size - tree->pgsize/2) {
                if (extend_mvpfile(tree) < 0) {
                    *error = MVP_FILETRUNCATE;
                    break;
                }
            }
            off_t offset = _mvptree_write(tree, node->internal.child_nodes[i], error, lvl+2);
            memcpy(&tree->buf[saved_pos++], &fileno, 1);
            memcpy(&tree->buf[saved_pos]  , &offset, sizeof(off_t));
            saved_pos += sizeof(off_t);
        }
    } else {
        *error = MVP_UNRECOGNIZED;
    }

    /*    msync(tree->buf, tree->size, MS_ASYNC);   */

    return start_pos;
}

MVPError mvptree_write(MVPTree *tree, const char *filename, int mode) {
    if (!tree || !tree->dist || !tree->node || !filename) {
        return MVP_ARGERR;
    }

    tree->fd = open(filename, O_CREAT|O_RDWR|O_TRUNC, mode);
    if (tree->fd < 0) {
        return MVP_FILEOPEN;
    }
    tree->pgsize = sysconf(_SC_PAGESIZE);
    if (ftruncate(tree->fd, tree->pgsize) < 0) {
        close(tree->fd);
        return MVP_FILETRUNCATE;
    }
    tree->size = tree->pgsize;
    char *buf  = (char*) mmap(NULL, tree->size, PROT_READ|PROT_WRITE, MAP_SHARED, tree->fd, 0);
    if (buf == NULL) {
        close(tree->fd);
        return MVP_MEMMAP;
    }
    off_t pos = 0;

    uint8_t bf = tree->branchfactor;
    uint8_t pl = tree->pathlength;
    uint8_t lc = tree->leafcap;
    uint8_t ht = (uint8_t) tree->node->internal.sv1->type;

    /* write header */
    memcpy(&buf[pos], tag, strlen(tag)+1);
    pos += strlen(tag)+1;
    memcpy(&buf[pos], &version, sizeof(version));
    pos += sizeof(version);
    memcpy(&buf[pos++], &bf, 1);
    memcpy(&buf[pos++], &pl, 1);
    memcpy(&buf[pos++], &lc, 1);
    memcpy(&buf[pos++], &ht, 1);

    tree->buf = buf;
    pos = HEADER_SIZE;
    tree->pos = pos;

    /* write nodes */
    MVPError error = MVP_SUCCESS;
    _mvptree_write(tree, tree->node, &error, 0);

    /* cleanup */
    if (msync(tree->buf, tree->size, MS_SYNC) < 0) {
        error = MVP_MEMMAP;
    }

    if (munmap(tree->buf, tree->size) < 0) {
        error = MVP_MUNMAP;
    }
    tree->buf = NULL;

    if (close(tree->fd) < 0) {
        error = MVP_FILECLOSE;
    }

    return error;
}

static MVPDP* read_datapoint(MVPTree *tree) {
    uint8_t active, idlen;
    uint32_t bytelength, datalength;

    memcpy(&active, &tree->buf[tree->pos], sizeof(active));
    tree->pos += sizeof(active);
    memcpy(&bytelength, &tree->buf[tree->pos], sizeof(bytelength));
    tree->pos += sizeof(bytelength);

    if (active == 0 && bytelength == 0) return NULL;

    MVPDP *dp = dp_alloc(tree->datatype);
    if (!dp) { return NULL; }

    dp->path = (float*) malloc(tree->pathlength*sizeof(float));
    if (!dp->path) { return NULL; }

    memcpy(&idlen, &tree->buf[tree->pos], sizeof(uint8_t));
    tree->pos += sizeof(uint8_t);
    dp->id = malloc(idlen+1);
    memcpy(dp->id, &tree->buf[tree->pos], idlen);
    tree->pos += idlen;
    dp->id[idlen] = '\0';
    memcpy(&datalength, &tree->buf[tree->pos], sizeof(uint32_t));
    tree->pos += sizeof(uint32_t);

    dp->datalen = datalength;
    dp->data = malloc(datalength*tree->datatype);
    memcpy(dp->data, &tree->buf[tree->pos], datalength*tree->datatype);
    tree->pos += datalength*tree->datatype;
    memcpy(dp->path, &tree->buf[tree->pos], tree->pathlength*sizeof(float));
    tree->pos += tree->pathlength*sizeof(float);

    return dp;
}

static Node* _mvptree_read_node(MVPTree *tree, MVPError *error, int lvl) {
    uint8_t node_type;
    Node *node = NULL;
    memcpy(&node_type, &tree->buf[tree->pos++], sizeof(uint8_t));

    if (node_type == LEAF_NODE) {
        uint32_t nbpoints;
        node = create_leaf(tree->leafcap);
        if (!node) {
            *error = MVP_NOLEAF;
            return node;
        }
        node->leaf.sv1 = read_datapoint(tree);
        node->leaf.sv2 = read_datapoint(tree);

        memcpy(&nbpoints,&tree->buf[tree->pos], sizeof(uint32_t));
        tree->pos += sizeof(uint32_t);
        node->leaf.nbpoints = nbpoints;

        off_t saved_pos = tree->pos;
        int i;
        off_t offset;
        for (i = 0; i < nbpoints; i++) {
            memcpy(&(node->leaf.d1[i]), &tree->buf[saved_pos], sizeof(float));
            saved_pos += sizeof(float);
            memcpy(&(node->leaf.d2[i]), &tree->buf[saved_pos], sizeof(float));
            saved_pos += sizeof(float);
            memcpy(&offset, &tree->buf[saved_pos], sizeof(off_t));
            saved_pos += sizeof(off_t);

            tree->pos = offset;
            node->leaf.points[i] = read_datapoint(tree);
        }
    } else if (node_type == INTERNAL_NODE) {
        int bf = tree->branchfactor;
        int lengthM1 = bf - 1;
        int lengthM2 = (bf - 1)*bf;
        int fanout   = bf*bf;
        uint8_t fileno;

        node = create_internal(bf);
        if (node == NULL) {
            *error = MVP_NOINTERNAL;
            return node;
        }
        node->internal.sv1 = read_datapoint(tree);
        node->internal.sv2 = read_datapoint(tree);

        memcpy(node->internal.M1, &tree->buf[tree->pos], lengthM1*sizeof(float));
        tree->pos += lengthM1*sizeof(float);
        memcpy(node->internal.M2, &tree->buf[tree->pos], lengthM2*sizeof(float));
        tree->pos += lengthM2*sizeof(float);

        int i;
        off_t offset, saved_pos = tree->pos;
        for (i = 0;i < fanout; i++) {
            memcpy(&fileno, &tree->buf[saved_pos], sizeof(fileno));
            saved_pos += sizeof(fileno);
            memcpy(&offset, &tree->buf[saved_pos], sizeof(offset));
            saved_pos += sizeof(offset);

            tree->pos = offset;
            node->internal.child_nodes[i] = _mvptree_read_node(tree, error, lvl+2);
            if (*error != MVP_SUCCESS) { break; }
        }
    } else {
        *error = MVP_UNRECOGNIZED;
    }
    return node;
}

MVPTree* mvptree_read(const char *filename, CmpFunc fnc, int branchfactor, int pathlength,
                               int leafcapacity,MVPError *error) {
    if (!error) { return NULL; }
    *error = MVP_SUCCESS;
    if (!filename || !fnc) {
        *error = MVP_ARGERR;
        return NULL;
    }
    MVPTree *tree = NULL;

    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        /* file not found return empty tree */
        tree = mvptree_alloc(NULL, fnc, branchfactor, pathlength, leafcapacity);
        *error = MVP_FILENOTFOUND;
        return tree;
    }

    struct stat file_info;
    if (fstat(fd, &file_info) < 0) {
        *error = MVP_FILEOPEN;
        return NULL;
    }
    off_t size = file_info.st_size;

    char *buf = (char*) mmap(NULL, size, PROT_READ, MAP_SHARED, fd, 0);
    if (buf == NULL) {
        *error = MVP_MEMMAP;
        close(fd);
        return NULL;
    }
    off_t pos = 0;
    char line[16];
    int v;
    uint8_t bf, pl, lc, ht;

    memcpy(line, &buf[pos], strlen(tag)+1);
    pos += strlen(tag)+1;
    memcpy(&v, &buf[pos], sizeof(int));
    pos += sizeof(int);
    memcpy(&bf, &buf[pos++], 1);
    memcpy(&pl, &buf[pos++], 1);
    memcpy(&lc, &buf[pos++], 1);
    memcpy(&ht, &buf[pos++], 1);

    tree = mvptree_alloc(NULL, fnc, bf, pl, lc);
    if (!tree) {
        *error = MVP_MEMALLOC;
        return NULL;
    }
    tree->pgsize = sysconf(_SC_PAGESIZE);
    tree->size = size;
    tree->buf = buf;
    tree->pos = HEADER_SIZE;
    tree->fd = fd;
    tree->datatype = (MVPDataType) ht;
    tree->dist = fnc;
    tree->node = _mvptree_read_node(tree, error, 0);

    if (munmap(buf, size) < 0) {
        *error = MVP_MUNMAP;
    }
    tree->buf = NULL;

    if (close(fd) < 0) {
        *error = MVP_FILECLOSE;
    }
    tree->buf = NULL;
    tree->pos = 0;
    tree->fd  = 0;

    return tree;
}

static MVPError _mvptree_print(FILE *stream, MVPTree *tree, Node *node, int lvl) {
    MVPError error = MVP_SUCCESS;
    Node *next_node = node;
    int bf = tree->branchfactor, lengthM1 = bf-1, lengthM2 = bf, fanout = bf*bf;

    if (next_node) {
        if (next_node->leaf.type == LEAF_NODE) {
            fprintf(stream, "LEAF%d  (%d points)\n", lvl, next_node->leaf.nbpoints);
            if (next_node->leaf.sv1) {
                fprintf(stream, "    sv1: %s\n", next_node->leaf.sv1->id);
            }
            if (next_node->leaf.sv2) {
                fprintf(stream, "    sv2: %s\n", next_node->leaf.sv2->id);
            }

            int i;
            for (i = 0; i < next_node->leaf.nbpoints; i++) {
                fprintf(stream, "        point[%d]: %s\n", i, next_node->leaf.points[i]->id);
            }
        } else if (next_node->internal.type == INTERNAL_NODE) {
            fprintf(stream, "INTERNAL%d\n", lvl);
            fprintf(stream, "  sv1: %s\n", next_node->internal.sv1->id);
            fprintf(stream, "  sv2: %s\n", next_node->internal.sv2->id);

            int i;
            for (i=0; i < lengthM1; i++) {
                fprintf(stream, "  M1[%d] = %.4f;", i, next_node->internal.M1[i]);
            }
            for (i=0; i < lengthM2; i++) {
                fprintf(stream, "  M2[%d] = %.4f;", i, next_node->internal.M2[i]);
            }
            fprintf(stream, "\n");

            for (i=0; i < fanout; i++) {
                error = _mvptree_print(stream, tree, node->internal.child_nodes[i], lvl+2);
                if (error != MVP_SUCCESS) { break; }
            }
        } else {
            error = MVP_UNRECOGNIZED;
        }
    } else {
        fprintf(stream, "NULL%d\n", lvl);
    }

    return error;
}

MVPError mvptree_print(FILE *stream, MVPTree *tree) {
    if (stream == NULL || tree == NULL) {
        return MVP_ARGERR;
    }

    MVPError err = _mvptree_print(stream, tree, tree->node, 0);
    if (err != MVP_SUCCESS) {
        fprintf(stream, "malformed tree: %s\n", mvp_errstr(err));
    }

    return err;
}
