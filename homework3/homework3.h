//
// INF 503 - Homework 3
// Jadon Fowler 3/23/22
//

#ifndef INF503_MAIN_H
#define INF503_MAIN_H

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <fstream>
#include <ctime>

typedef char *GenomeLine;

class Bacillus {
public:
    char *data;
    int length;

    Bacillus(char *fileName);

    char *randomStart(int size);

    void randomize(int percent);

    ~Bacillus();
};


class FASTAreadset_LL {
public:
    GenomeLine line;
    FASTAreadset_LL *next;

    FASTAreadset_LL(GenomeLine line, FASTAreadset_LL *parent);

    void insert(GenomeLine line, int length);

    ~FASTAreadset_LL();
};

class FASTAreadset_HT {
private:
    unsigned int capacity;
    FASTAreadset_LL **chains;
    unsigned int collisions;
public:
    FASTAreadset_HT(unsigned int m);

    unsigned int radixHash(GenomeLine line, int length);

    bool search(GenomeLine line, int length);

    void insert(GenomeLine line, int length);

    ~FASTAreadset_HT();

    double readFile(Bacillus *bacillus, int fragmentLength);

    unsigned int getCollisions();
};

GenomeLine generateGenomeLine(int length);

#endif //INF503_MAIN_H
