//
// INF 503 - Homework 5 + 6
// Jadon Fowler
//

#ifndef INF503_HOMEWORK5_H
#define INF503_HOMEWORK5_H

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <fstream>
#include <ctime>

// ACTG
const int keys = 4;

int indexFrom(char c);

/**
 * Dataset of genomic lines such as SARS-COV2
 */
class Dataset {
public:
    char *data;
    int length;

    Dataset(char *fileName);

    char *randomStart(int size);

    void randomize(int percent);

    ~Dataset();
};

class trie_node {
public:
    char value;
    bool terminal;
    trie_node *children[keys];

    trie_node(char v);

    ~trie_node();

    void insert(char *line, int length, int offset);

    bool search(char *line, int length, int offset, int fuzz);

    int countNodes();

    int countTerminalNodes();
};

class prefix_trie {
public:
    trie_node *root;

    prefix_trie();

    ~prefix_trie();

    void insert(char *line, int length);

    bool search(char *line, int length);

    int countNodes();
};

class suffix_tree {
public:
    trie_node *root;

    suffix_tree();

    suffix_tree(Dataset *dataset, int fragmentLength);

    ~suffix_tree();

    void insert(char *line, int length);

    bool search(char *line, int length);

};

#endif //INF503_HOMEWORK5_H
