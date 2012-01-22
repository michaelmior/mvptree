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
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "CImg.h"

extern "C" {
    #include "mvptree.h"
}

#define MVP_BRANCHFACTOR 2
#define MVP_PATHLENGTH   5
#define MVP_LEAFCAP     25
#define PI 3.141592

static unsigned long long nbcalcs = 0;

using namespace cimg_library;

CImg<float>* ph_dct_matrix(const int N){
    CImg<float> *ptr_matrix = new CImg<float>(N, N, 1, 1, 1/sqrt((float) N));
    const float c1 = sqrt(2.0/(float)N);
    for (int x=0; x<N; x++){
        for (int y=1; y<N; y++){
            *ptr_matrix->data(x,y) = c1*cos((PI/2/N)*y*(2*x+1));
        }
    }
    return ptr_matrix;
}

int ph_dct_imagehash(const char* file,uint64_t &hash){

    if (!file){
        return -1;
    }

    CImg<uint8_t> src;
    try {
        src.load(file);
    } catch (CImgIOException ex){
        return -1;
    }

    CImg<float> meanfilter(7, 7, 1, 1, 1);
    CImg<float> img;
    if (src.spectrum() == 3){
        img = src.RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else if (src.spectrum() == 4){
        int width = img.width();
        int height = img.height();
        int depth = img.depth();

        img = src.crop(0,0,0,0,width-1,height-1,depth-1,2).RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else {
        img = src.channel(0).get_convolve(meanfilter);
    }

    img.resize(32,32);
    CImg<float> *C  = ph_dct_matrix(32);
    CImg<float> Ctransp = C->get_transpose();

    CImg<float> dctImage = (*C)*img*Ctransp;

    CImg<float> subsec = dctImage.crop(1,1,8,8).unroll('x');;

    float median = subsec.median();
    uint64_t one = 0x0000000000000001;
    hash = 0x0000000000000000;
    for (int i=0;i< 64;i++){
        float current = subsec(i);
        if (current > median) {
            hash |= one;
        }
        one = one << 1;
    }

    delete C;
    return 0;
}

float hamming_distance(MVPDP *pointA, MVPDP *pointB){
    if (!pointA || !pointB || pointA->datalen != pointB->datalen) {
        return -1.0f;
    }

    uint64_t a = *((uint64_t*) pointA->data);
    uint64_t b = *((uint64_t*) pointB->data);

    uint64_t x = a^b;
    const uint64_t m1  = 0x5555555555555555ULL;
    const uint64_t m2  = 0x3333333333333333ULL;
    const uint64_t h01 = 0x0101010101010101ULL;
    const uint64_t m4  = 0x0f0f0f0f0f0f0f0fULL;
    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4)) & m4;

    float result = (float)((x*h01)>>56);
    nbcalcs++;
    return result;
}


char** ph_readfilenames(const char *dirname, int &count) {
    count = 0;
    struct dirent *dir_entry;

    if (!dirname) { return NULL; }
    DIR *dir = opendir(dirname);
    if (!dir) { return NULL; }

    /*count files */
    while ((dir_entry = readdir(dir)) != NULL){
	if (strcmp(dir_entry->d_name, ".") && strcmp(dir_entry->d_name, ".."))
	    count++;
    }

    /* alloc list of files */
    char **files = (char**) malloc(count*sizeof(*files));
    if (!files) { return NULL; }

    errno = 0;
    int index = 0;
    char path[1024];
    path[0] = '\0';
    rewinddir(dir);
    while ((dir_entry = readdir(dir)) != 0) {
        if (strcmp(dir_entry->d_name, ".") && strcmp(dir_entry->d_name, "..")) {
            strcat(path, dirname);
            strcat(path, "/");
            strcat(path, dir_entry->d_name);
            files[index++] = strdup(path);
        }
        path[0] = '\0';
    }

    if (errno) { return NULL; }

    closedir(dir);
    return files;
}


int main(int argc, char **argv){
    if (argc < 3) {
        printf("not enough input args\n");
        printf(" %s  command filename [directory] [radius]\n\n", argv[0]);
        printf("  command    - command - e.g. 'add', 'query' or 'print'\n");
        printf("  filename   - file from which to read the tree\n");
        printf("  directory  - directory (for add and query)\n");
        printf("  radius     - radius for query operation (default = 21.0)\n");
        return 0;
    }

    const char *command  = argv[1];
    const char *filename = argv[2];
    const char *dirname  = argv[3];
    const float radius = argv[4] ? atof(argv[4]) : 21.0f;

    const int knearest = 5;

    printf("command  - %s\n", command);
    printf("filename - %s\n", filename);
    printf("dir      - %s\n", dirname);
    printf("radius   - %f\n", radius);
    printf("knearest - %d\n", knearest);

    CmpFunc distance_func = hamming_distance;

    MVPError err;
    MVPTree *tree = mvptree_read(
        filename,
        distance_func,
        MVP_BRANCHFACTOR,
        MVP_PATHLENGTH,
        MVP_LEAFCAP,
        &err
    );
    assert(tree);

    if (!strncasecmp(command, "add", 3) || !strncasecmp(command, "query", 3)) {
        int nbfiles;
        char **files = ph_readfilenames(dirname, nbfiles);
        assert(files);

        fprintf(stdout, "\n %d files in %s\n\n", nbfiles, dirname);

        MVPDP **points = (MVPDP**) malloc(nbfiles*sizeof(MVPDP*));
        assert(points);

        int count = 0;
        uint64_t hashvalue;
        for (int i=0; i < nbfiles; i++){
            char *name = strrchr(files[i], '/') + 1;

            if (ph_dct_imagehash(files[i], hashvalue) < 0){
                printf("Unable to get hash value.\n");
                continue;
            }
            printf("(%d) %llx %s\n", i, (unsigned long long) hashvalue, files[i]);

            points[count] = dp_alloc(UINT64ARRAY);
            points[count]->id = strdup(name);
            points[count]->data = malloc(1*UINT64ARRAY);
            points[count]->datalen = 1;
            memcpy(points[count]->data, &hashvalue, UINT64ARRAY);
            count++;
        }

        printf("\n");

        if (!strncasecmp(command, "add", 3)) {
            printf("Add %d hashes to tree.\n", count);
            MVPError error = mvptree_add(tree, points, count);
            if (error != MVP_SUCCESS){
                printf("Unable to add hash values to tree - %s\n", mvp_errstr(error));
                goto cleanup;
            }

            printf("Save file.\n");
            error = mvptree_write(tree, filename, 00755);
            if (error != MVP_SUCCESS){
                printf("Unable to save file - %s\n", mvp_errstr(error));
                goto cleanup;
            }
        } else if (!strncasecmp(command, "query", 3)){
            unsigned int nbresults;
            for (int i=0; i<count; i++){
                printf("(%d) looking up %s ...\n", i, files[i]);
                nbcalcs = 0;
                MVPDP **results = mvptree_retrieve(tree, points[i], knearest,\
                                                                       radius, &nbresults, &err);
                printf("-----------%u results (%llu calcs)--------------\n", nbresults, nbcalcs);
                for (unsigned int j=0; j < nbresults; j++) {
                    printf("(%d) %s\n", j, results[j]->id);
                }
                printf("-----------------------------------------------\n");
                printf("Hit enter key.\n");
                getchar();
                free(results);
            }
        }

cleanup:
        for (int i=0; i<nbfiles; i++){
            free(files[i]);
        }
        free(files);
    } else if (!strncasecmp(command, "print", 3)) {
        printf("-----------------------print-------------------------\n");
        mvptree_print(stdout, tree);
        printf("-----------------------------------------------------\n\n");
    }

    mvptree_clear(tree, free);

    return 0;
}
