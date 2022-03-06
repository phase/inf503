//
// Created by jadon on 2/13/22.
//

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <fstream>
#include <ctime>

using namespace std;

class FASTA_readset {
private:

    /**
     * Resizes the arrays to double the capacity they were.
     */
    void resize() {
        int oldCapacity = capacity;
        capacity += capacity / 2;
        // cout << "Resizing to " << capacity << endl;

        // allocate a new array that is double the size
        int *newCounts = new int[capacity];
        std::copy(counts, counts + oldCapacity, newCounts);

        // replace count with the new array
        delete counts;
        counts = newCounts;

        // do the same for the lines
        char **newLines = new char *[capacity];
        std::copy(lines, lines + oldCapacity, newLines);

        delete lines;
        lines = newLines;
    }

    /**
     * Generic function to swap two pointers
     */
    template<typename Type>
    void swap(Type *a, Type *b) {
        Type temp = *a;
        *a = *b;
        *b = temp;
    }

    /**
     * Normal quicksort that sorts the counts and lines
     */
    void sort(char **lines, int *counts, int size) {
        int pivot = 0;
        if (size <= 1) {
            return;
        }

        for (int i = 0; i < size; i++) {
            if (strcmp(lines[i], lines[size - 1]) < 0) {
                swap(lines + i, lines + pivot);
                swap(counts + i, counts + pivot);
                pivot++;
            }
        }

        swap(lines + pivot, lines + size - 1);
        swap(counts + pivot, counts + size - 1);

        sort(lines, counts, pivot++);
        sort(lines + pivot, counts + pivot, size - pivot);
    }

public:
    int size;
    int capacity;
    int *counts;
    char **lines;

    FASTA_readset() {
        size = 0;
        capacity = 1000;
        counts = new int[capacity];
        lines = new char *[capacity];
    }

    FASTA_readset(char *file_name, int dataset_num) {
        ifstream input;
        input.open(file_name);
        int totalLines = 0;

        // temporary buffer for header... 1000 characters will do
        char *header = new char[1000];
        // temporary buffer for genomic data... 1000 characters will do
        char *line = new char[1000];

        int linesProcessed = 0;
        while (!input.eof()) {
            input >> header; // read in the header line
            if (input.eof()) break;
            input >> line; // read in the read line
            linesProcessed += 2;

            // parse header using a state machine
// parse header using a state machine

            bool r = false;
            char c; // current char
            int currentDataset = 0;
            for (int h = 0; c = header[h], c != '\0'; h++) {
                if (h == 0) {
                    if (c != '>') cerr << "INVALID SYNTAX or BAD PARSER" << endl;
                } else if (h == 1) {
                    if (c != 'R') cerr << "INVALID SYNTAX or BAD PARSER" << endl;
                    r = true;
                } else if (c == '_') {
                    r = false;
                    // about to start or in the middle of parsing the fragment copies
                } else if (r) {
                    // we just saw an R, ignore for now
                } else if (c >= '0' && c <= '9') {
                    // this is a number, can be multiple digits
                    // this is how many copies of the fragment is in each dataset
                    char *numBuf = new char[5];
                    char w;
                    int i = 0;
                    for (int n = h; w = header[n], w >= '0' && w <= '9'; n++) {
                        numBuf[i++] = w;
                    }
                    h += i - 1;
                    int copies = atoi(numBuf);
                    delete[] numBuf;
                    if (copies > 0 && currentDataset == dataset_num) {
                        totalLines++;
                    }
                    currentDataset++;
                }
            }

            if (linesProcessed % 10000000 == 0) {
                cout << "...found " << linesProcessed << " fragments for dataset " << dataset_num << " so far (found "
                     << totalLines << ")..." << endl;
            }
        }
        input.close();
        delete[] header;
        delete[] line;

        size = 0;
        capacity = totalLines;
        counts = new int[capacity];
        lines = new char *[capacity];
    }

    ~FASTA_readset() {
        delete counts;
        for (int i = 0; i < size; i++) {
            delete lines[i];
        }
        delete lines;
    }

    /**
     * Adds a line to this set
     * @param count the amount of times this line is in the set
     * @param line the line of characters
     */
    void add(int count, char *line) {
        if (size + 1 >= capacity) {
            // if we're trying to add to full buffers, resize them
            resize();
        }
        counts[size] = count;
        lines[size] = line;
        size++;
    }

    int getUniqueLines() {
        return size;
    }

    int getTotalLines() {
        // add up all the copy counts
        int total = 0;
        for (int i = 0; i < size; i++) {
            total += counts[i];
        }
        return total;
    }

    int findLine(char *needle) {
        // go over every line
        for (int i = 0; i < size; i++) {
            char *hay = lines[i];
            if (strcmp(needle, hay) == 0) {
                return i;
            }
        }
        // we couldn't find the string
        return -1;
    }

    int binarySearch(char *needle) {
        int l = 0;
        int r = size - 1;
        while (l <= r) {
            // C standard says integer div truncates towards 0
            int m = (l + r) / 2;
            char *hay = lines[m];
            // compare the strings lexicographically
            int compare = strcmp(needle, hay);
            if (compare > 0) {
                // the needle is larger than the hay
                // move the lower bound up
                l = m + 1;
            } else if (compare < 0) {
                // the hay is larger than the needle
                // move the upper bound up
                r = m - 1;
            } else {
                // the strings are equal, we found the index
                return m;
            }
        }
        // we couldn't find the string
        return -1;
    }

    void sort() {
        sort(lines, counts, size);
    }

};

int main(int argc, char **argv) {
    if (argc != 3) {  // unexpected program call
        cout << endl << endl << "==========================" << endl;
        cout << "Error: 2 input parameters expected" << endl;
        cout << "Proper usage is:" << endl;
        cout << "./homework <problem-flag> <filepath>" << endl;
        cout << "Example:" << endl;
        cout << "./homework problem1A /scratch/vyf2/HW1/sample_hw_dataset.fa" << endl;
        cout << "==========================" << endl << endl;
        cout << "exiting..." << endl;
        exit(-1);
    } else {
        cout << "The number of arguments passed: " << argc << endl;
        cout << "The first argument is: " << argv[0] << endl;
        cout << "The second argument is: " << argv[1] << endl;
        cout << "The third argument is: " << argv[2] << endl;
    }

    char *problem = argv[1];
    char *file_name = argv[2];

    // allocate the datasets
    int datasetCount = 14;
    FASTA_readset *datasets[datasetCount];
    for (int i = 0; i < datasetCount; i++) {
        datasets[i] = new FASTA_readset(file_name, i);
    }

    // create filestream to read the file
    ifstream input;
    // initialize the filestream by pointing it to the right file
    input.open(file_name);

    // temporary buffer for header... 1000 characters will do
    char *header = new char[1000];
    // temporary buffer for genomic data... 1000 characters will do
    char *line = new char[1000];

    int linesProcessed = 0;
    while (!input.eof()) {
        input >> header; // read in the header line
        if (input.eof()) break;
        input >> line; // read in the read line
        linesProcessed += 2;

        // parse header using a state machine

        bool r = false;
        char c; // current char
        int *fragmentCopies = new int[datasetCount];
        for (int f = 0; f < datasetCount; f++) {
            fragmentCopies[f] = 0;
        }
        int currentDataset = 0;
        for (int h = 0; c = header[h], c != '\0'; h++) {
            if (h == 0) {
                if (c != '>') cerr << "INVALID SYNTAX or BAD PARSER" << endl;
            } else if (h == 1) {
                if (c != 'R') cerr << "INVALID SYNTAX or BAD PARSER" << endl;
                r = true;
            } else if (c == '_') {
                r = false;
                // about to start or in the middle of parsing the fragment copies
            } else if (r) {
                // we just saw an R, ignore for now
            } else if (c >= '0' && c <= '9') {
                // this is a number, can be multiple digits
                // this is how many copies of the fragment is in each dataset
                char *numBuf = new char[5];
                char w;
                int i = 0;
                for (int n = h; w = header[n], w >= '0' && w <= '9'; n++) {
                    numBuf[i++] = w;
                }
                h += i - 1;
                int copies = atoi(numBuf);
                fragmentCopies[currentDataset] = copies;
                currentDataset++;
            }
        }

        for (int f = 0; f < datasetCount; f++) {
            FASTA_readset *dataset = datasets[f];
            int copies = fragmentCopies[f];
            if (copies > 0) {
                int lineLength;
                for (lineLength = 0; line[lineLength] != '\0'; lineLength++);
                char *lineCopy = new char[lineLength];
                for (lineLength = 0; line[lineLength] != '\0'; lineLength++) {
                    lineCopy[lineLength] = line[lineLength];
                }
                dataset->add(copies, lineCopy);
            }
            /*cout << "Dataset " << f << " has " << dataset->getUniqueLines() << " unique lines and "
                 << dataset->getTotalLines() << " total lines.\n";*/
        }
        delete[] fragmentCopies;

        if (linesProcessed % 10000000 == 0) {
            cout << "...loaded " << linesProcessed << " fragments so far..." << endl;
        }
    }
    cout << "Done loading fragments into dataset objects." << endl;

    /*
    // testing code
    for (int f = 0; f < datasetCount; f++) {
        FASTA_readset *dataset = datasets[f];
        cout << "Dataset " << f << " has " << dataset->getUniqueLines() << " unique lines and "
             << dataset->getTotalLines() << " total lines.\n" << endl;
        dataset->sort();
        for (int i = 0; i < dataset->size; i++) {
            cout << dataset->lines[i] << endl;
        }
    }

    for (int j = 0; j < datasetCount - 1; j++) {
        FASTA_readset *first = datasets[j];
        FASTA_readset *second = datasets[j + 1];
        for (int i = 0; i < first->size; i++) {
            char *needle = first->lines[i];
            int result = second->binarySearch(needle);
            if (result != -1) {
                cout << "Found " << needle << endl;
                cout << "  index in dataset " << j << ": " << i << endl;
                cout << "  index in dataset " << (j+1) << ": " << result << endl;
            } else {
                //cout << "Couldn't find " << needle << endl;
            }
        }
    }
    */
    cout << "\n" << endl;
    if (strcmp("problem1A", problem) == 0) {
        cout << "Problem 1A:" << endl;
        cout << " How many unique sequence fragments are in each of the 14 datasets?" << endl;
        cout << " How many total sequence fragments are in each dataset (i.e. when you consider copy numbers)?\n"
             << endl;
        for (int f = 0; f < datasetCount; f++) {
            FASTA_readset *dataset = datasets[f];
            cout << "Dataset " << (f + 1) << " has " << dataset->getUniqueLines() << " unique fragments and "
                 << dataset->getTotalLines() << " total fragments." << endl;
        }
    } else if (strcmp("problem1B", problem) == 0) {
        cout << "Problem 1B:" << endl;
        cout << " How many sequence fragments in dataset 1 are also in dataset 2?\n" << endl;

        clock_t t;
        t = clock();

        FASTA_readset *first = datasets[0];
        FASTA_readset *second = datasets[1];
        int sharedFragments = 0;
        double lastPercent = 0;
        for (int i = 0; i < first->size; i++) {
            char *needle = first->lines[i];
            int result = second->findLine(needle);
            if (result != -1) {
                sharedFragments += second->counts[result];
            }

            // debug
            double percentDone = floor(((double) i / (double) first->size) * 100.0);
            if (percentDone > lastPercent) {
                clock_t soFar = clock() - t;
                double time_taken = ((double) soFar) / CLOCKS_PER_SEC;
                cout << "..." << percentDone << "% (" << i << "/" << first->size << ") done in " << time_taken << "s..."
                     << endl;
                lastPercent = floor(percentDone) + 5.0;
            }
        }

        t = clock() - t;
        double time_taken = ((double) t) / CLOCKS_PER_SEC;
        cout << "There are " << sharedFragments << " fragments from dataset 1 in dataset 2." << endl;
        cout << "Searching took " << time_taken << "s." << endl;
    } else if (strcmp("problem1C", problem) == 0) {
        cout << "Problem 1C:" << endl;
        cout << " How many sequence fragments in dataset 1 are also in dataset 2?\n" << endl;

        clock_t t;
        t = clock();

        FASTA_readset *first = datasets[0];
        FASTA_readset *second = datasets[1];

        cout << "...sorting the first dataset..." << endl;
        first->sort();
        cout << "...sorting the second dataset..." << endl;
        second->sort();
        int sharedFragments = 0;
        double lastPercent = 0;
        for (int i = 0; i < first->size; i++) {
            char *needle = first->lines[i];
            int result = second->binarySearch(needle);
            if (result != -1) {
                sharedFragments += second->counts[result];
            }

            // debug
            double percentDone = floor(((double) i / (double) first->size) * 100.0);
            if (percentDone > lastPercent) {
                clock_t soFar = clock() - t;
                double time_taken = ((double) soFar) / CLOCKS_PER_SEC;
                cout << "..." << percentDone << "% (" << i << "/" << first->size << ") done in " << time_taken
                     << "s with " << sharedFragments << " shared fragments..."
                     << endl;
                lastPercent = floor(percentDone) + 5.0;
            }
        }

        t = clock() - t;
        double time_taken = ((double) t) / CLOCKS_PER_SEC;
        cout << "There are " << sharedFragments << " fragments from dataset 1 in dataset 2." << endl;
        cout << "Searching took " << time_taken << "s." << endl;
    }

    return 0;
}
