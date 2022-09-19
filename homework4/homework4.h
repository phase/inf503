//
// Created by jado on 4/14/22.
//

#ifndef INF503_HOMEWORK4_H
#define INF503_HOMEWORK4_H


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
    unsigned int collisions;
public:
    FASTAreadset_LL **chains;
    unsigned int capacity;

    FASTAreadset_HT(unsigned int m);

    unsigned int radixHash(GenomeLine line, int length);

    bool search(GenomeLine line, int length);

    void insert(GenomeLine line, int length);

    ~FASTAreadset_HT();

    double readFile(Bacillus *bacillus, int fragmentLength);

    unsigned int getCollisions();
};

GenomeLine generateGenomeLine(int length);


class BLAST_DB {
public:
    FASTAreadset_HT *table;

    BLAST_DB(int capacity);

    ~BLAST_DB();

    double readFile(Bacillus *bacillus);

    int search(char *line, int fragmentLength);

    int nwScore(char *a, char *b, int lengthA, int lengthB);
};


#endif //INF503_HOMEWORK4_H
