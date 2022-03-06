//
// Created by jadon on 3/6/22.
// [jaf582@wind ~ ]$ srun --mem=5000 ./a.out problem1A /common/contrib/classroom/inf503/hw_dataset.fa /common/contrib/classroom/inf503/test_genome.fasta
//

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <fstream>
#include <ctime>

using namespace std;

// limit for reading lines - good for testing
// setting this to -1 will remove the limit
constexpr int MAX_LINES = -1;

class FASTA_readset_LL {
private:
    char *data;
    FASTA_readset_LL *next;

public:
    FASTA_readset_LL() {
        this->data = nullptr;
        this->next = nullptr;
    }

    void readFile(char *file_name) {
        // create filestream to read the file
        ifstream input;
        // initialize the filestream by pointing it to the right file
        input.open(file_name);

        // temporary buffer for genomic data... 1000 characters will do
        char *line = new char[1000];
        FASTA_readset_LL *head = this;
        int linesProcessed = 0;
        while (!input.eof()) {
            input >> line; // read in the read line
            if (line[0] == '>') {
                continue;
            }
            size_t len = strlen(line);
            char *copy = new char[len];
            strncpy(copy, line, len);

            head->data = copy;
            head->next = new FASTA_readset_LL();
            head = head->next;

            linesProcessed++;
            if (linesProcessed % 10000000 == 0) {
                cout << "...loaded " << linesProcessed << " fragments so far..." << endl;
            }
            if (MAX_LINES != -1 && linesProcessed >= MAX_LINES) {
                break;
            }
        }
        delete[] line;
        input.close();
    }

    FASTA_readset_LL(char *file_name) {
        readFile(file_name);
    }

    ~FASTA_readset_LL() {
        if (data != nullptr) {
            delete[] data;
        }
        if (next != nullptr) {
            delete next;
        }
    }

    FASTA_readset_LL *find(char *needle, int length) {
        auto head = this;
        // traverse the list
        while (head != nullptr) {
            if (head->data != nullptr) {
                if (strncmp(head->data, needle, length) == 0) {
                    // we got a hit
                    return head;
                }
            }
            head = head->next;
        }
        // we couldn't find it
        return nullptr;
    }
};

int main(int argc, char **argv) {
    if (argc < 3) {  // unexpected program call
        cout << endl << endl << "==========================" << endl
             << "INF 503 Homework 2 by Jadon Fowler" << endl
             << "Proper usage is:" << endl
             << "- Homework 1A:" << endl
             << "   ./homework2 problem1A /common/contrib/classroom/inf503/hw_dataset.fa" << endl
             << "- Homework 1B:" << endl
             << "   ./homework2 problem1B /common/contrib/classroom/inf503/hw_dataset.fa /common/contrib/classroom/inf503/test_genome.fasta"
             << endl << "==========================" << endl << endl;
        exit(-1);
    }

    char *problem = argv[1];
    char *file_name = argv[2];
    cout << "Reading in FASTA data..." << endl;
    auto fasta = new FASTA_readset_LL(file_name);
    cout << "Done reading in FASTA data." << endl;

    if (strcmp("problem1A", problem) == 0) {
        auto needles = {
                "CTAGGTACATCCACACACAGCAGCGCATTATGTATTTATTGGATTTATTT",
                "GCGCGATCAGCTTCGCGCGCACCGCGAGCGCCGATTGCACGAAATGGCGC",
                "CGATGATCAGGGGCGTTGCGTAATAGAAACTGCGAAGCCGCTCTATCGCC",
                "CGTTGGGAGTGCTTGGTTTAGCGCAAATGAGTTTTCGAGGCTATCAAAAA",
                "ACTGTAGAAGAAAAAAGTGAGGCTGCTCTTTTACAAGAAAAAGTNNNNNN"
        };
        for (auto needle: needles) {
            if (needle == nullptr) break;
            auto result = fasta->find((char *) needle, 50);
            cout << needle;
            if (result != nullptr) {
                cout << " IS in the dataset.";
            } else {
                cout << " is NOT in the dataset.";
            }
            cout << endl;
        }
    } else if (strcmp("problem1B", problem) == 0) {
        if (argc > 3) {
            char *bacillus = argv[3];

            // create filestream to read the file
            ifstream input;
            // initialize the filestream by pointing it to the right file
            input.open(bacillus);

            int totalChars = 0;
            char *line = new char[1000];
            int linesProcessed = 0;
            while (!input.eof()) {
                input >> line;
                // ignore
                if (line[0] == '>') {
                    cout << "...skipping line: " << line << endl;
                    continue;
                }
                totalChars += strlen(line);

                linesProcessed++;
                if (linesProcessed % 10000000 == 0) {
                    cout << "...loaded " << linesProcessed << " fragments so far..." << endl;
                }
                if (MAX_LINES != -1 && linesProcessed >= MAX_LINES) {
                    break;
                }
            }

            input.close();
            input.open(bacillus);
            char *dataset = new char[totalChars];
            linesProcessed = 0;
            int currentIndex = 0;
            while (!input.eof()) {
                input >> line;
                // ignore
                if (line[0] == '>') continue;
                strcpy(dataset + currentIndex, line);
                currentIndex += strlen(line);

                linesProcessed++;
                // debug
                if (linesProcessed % 10000000 == 0) {
                    cout << "...loaded " << linesProcessed << " fragments so far..." << endl;
                }
                if (MAX_LINES != -1 && linesProcessed >= MAX_LINES) {
                    break;
                }
            }
            int fragmentLength = 50;
            int totalSubstrings = totalChars - fragmentLength;
            cout << totalSubstrings << " fragments can be made from the Bacillus dataset." << endl;
            int sharedFragments = 0;
            clock_t t = clock();
            for (int i = 0; i < totalSubstrings; i++) {
                // start at the ith character in the dataset
                char *currentFragment = dataset + i;
                auto result = fasta->find(currentFragment, fragmentLength);
                if (result != nullptr) {
                    sharedFragments++;
                }
                if (i % 50 == 0) {
                    auto end = clock() - t;
                    auto time_taken = ((double) end) / CLOCKS_PER_SEC;
                    cout << "...processed " << i << " fragments in " << (time_taken) <<
                         "s with " << sharedFragments << " shared fragments..." << endl;
                }
            }

            cout << "Shared fragments: " << sharedFragments << endl;
        }
    }
    cout << "Done!" << endl;
    delete fasta;

    return 0;
}
