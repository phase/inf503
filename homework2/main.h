//
// Created by jadon on 3/6/22.
//

#ifndef INF503_MAIN_H
#define INF503_MAIN_H

// main linked list data structure for FASTA data
class FASTA_readset_LL {
private:
    char *data;
    FASTA_readset_LL *next;

public:
    // empty constructor
    FASTA_readset_LL();

    // read a file and create nodes for every line
    void readFile(char *file_name);

    // constructor that calls readFile
    FASTA_readset_LL(char *file_name);

    // destructor safely deletes data & next
    ~FASTA_readset_LL();

    // iterate over the list in search of the needle and return a (null) pointer of the node
    FASTA_readset_LL *find(char *needle, int length);
};


#endif //INF503_MAIN_H
