//
// Created by jado on 4/14/22.
//

#include "homework4.h"

using namespace std;

// limit for reading lines - good for testing
// setting this to -1 will remove the limit
constexpr int MAX_LINES = -1;

Bacillus::Bacillus(char *fileName) {
    clock_t totalTimer = clock();
    // create filestream to read the file
    ifstream input;
    // initialize the filestream by pointing it to the right file
    input.open(fileName);

    // read every line to determine how large the genome is
    length = 0;
    char *line = new char[1000];
    int linesProcessed = 0;
    while (!input.eof()) {
        input >> line;
        // ignore
        if (line[0] == '>') {
            cout << "...skipping line: " << line << endl;
            continue;
        }
        length += strlen(line);

        linesProcessed++;
        if (linesProcessed % 10000000 == 0) {
            cout << "...loaded " << linesProcessed << " fragments so far..." << endl;
        }
        if (MAX_LINES != -1 && linesProcessed >= MAX_LINES) {
            break;
        }
    }

    input.close();

    // copy all lines into the dataset array
    input.open(fileName);
    data = new char[length];
    linesProcessed = 0;
    int currentIndex = 0;
    while (!input.eof()) {
        input >> line;
        // ignore
        if (line[0] == '>') continue;
        strcpy(data + currentIndex, line);
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
    delete[] line;
    input.close();

    // report time
    auto end = clock() - totalTimer;
    auto timeTaken = ((double) end) / CLOCKS_PER_SEC;
    cout << "Loaded Bacillus in " << timeTaken << " seconds." << endl;
}

char *Bacillus::randomStart(int size) {
    return data + (rand() % (length - size));
}

void Bacillus::randomize(int percent) {
    for (int i = 0; i < length; i++) {
        if ((rand() % 100) < percent) {
            switch (rand() % 4) {
                case 0:
                    data[i] = 'A';
                    break;
                case 1:
                    data[i] = 'C';
                    break;
                case 2:
                    data[i] = 'T';
                    break;
                case 3:
                default:
                    data[i] = 'G';
                    break;
            }
        }
    }
}

Bacillus::~Bacillus() {
    delete data;
}

FASTAreadset_LL::FASTAreadset_LL(GenomeLine line, FASTAreadset_LL *parent) {
    this->line = line;
    this->next = nullptr;
    if (parent != nullptr) {
        parent->next = this;
    }
}

void FASTAreadset_LL::insert(GenomeLine line, int length) {
    auto current = this;
    while (current != nullptr) {
        if (strncmp(line, current->line, length) == 0) {
            // we found it, no need to add it again
            return;
        }
        if (current->next == nullptr) {
            // we found the last one, exit early
            break;
        }
        current = current->next;
    }
    // current now holds the last item in this LL
    // create a new node and set its parent to current
    new FASTAreadset_LL(line, current);
}

FASTAreadset_LL::~FASTAreadset_LL() {
    if (line != nullptr) {
        delete line;
    }
    if (next != nullptr) {
        delete next;
    }
}

FASTAreadset_HT::FASTAreadset_HT(unsigned int m) {
    this->capacity = m;
    this->chains = new FASTAreadset_LL *[m]{};
}

unsigned int FASTAreadset_HT::radixHash(GenomeLine line, int length) {
    unsigned int total = 1;
    int base = 291;
    unsigned int power = 1;
    for (int i = 0; i < length; i++) {
        // add the ascii value of every char together
        unsigned int charValue = line[i];
        total += charValue * power;
        power *= base;
    }
    return total % capacity;
}

bool FASTAreadset_HT::search(GenomeLine line, int length) {
    // get the proper chain for this line
    auto h = radixHash(line, length);
    auto chain = this->chains[h];
    if (chain != nullptr) {
        // check every node in the LL to see if the line is in the chain
        auto current = chain;
        while (current != nullptr) {
            if (strncmp(line, current->line, length) == 0) {
                // we found it!
                return true;
            }
            current = current->next;
        }
    }
    return false;
}

void FASTAreadset_HT::insert(GenomeLine line, int length) {
    auto h = radixHash(line, length);
    auto chain = this->chains[h];
    if (chain == nullptr) {
        // we don't yet have a chain for this hash
        // create one with the line to insert
        this->chains[h] = new FASTAreadset_LL(line, nullptr);
    } else {
        // collision! insert to the end of the list
        chain->insert(line, length);
        collisions++;
    }
}

FASTAreadset_HT::~FASTAreadset_HT() {
    delete[] chains;
}

double FASTAreadset_HT::readFile(Bacillus *bacillus, int fragmentLength) {
    clock_t totalTimer = clock();

    int totalSubstrings = bacillus->length - fragmentLength;

    cout << totalSubstrings << " fragments can be made from the Bacillus dataset." << endl;
    clock_t t = clock();
    for (int i = 0; i < totalSubstrings; i++) {
        // start at the ith character in the dataset
        char *currentFragment = bacillus->data + i;
        // insert a copy of this n-mer
        char *copy = new char[fragmentLength];
        strncpy(copy, currentFragment, fragmentLength);
        insert(copy, fragmentLength);
        if (i % 1000000 == 0) {
            auto end = clock() - t;
            auto timeTaken = ((double) end) / CLOCKS_PER_SEC;
            cout << "...processed " << i << " fragments in " << (timeTaken) << "s with " << collisions
                 << " collisions" << endl;
        }
    }

    auto end = clock() - totalTimer;
    auto timeTaken = ((double) end) / CLOCKS_PER_SEC;
    return timeTaken;
}

unsigned int FASTAreadset_HT::getCollisions() {
    return collisions;
}

GenomeLine generateGenomeLine(int length) {
    auto line = new char[length];
    for (int i = 0; i < length; i++) {
        switch (rand() % 4) {
            case 0:
                line[i] = 'A';
                break;
            case 1:
                line[i] = 'C';
                break;
            case 2:
                line[i] = 'T';
                break;
            case 3:
            default:
                line[i] = 'G';
                break;
        }
    }
    return line;
}

BLAST_DB::BLAST_DB(int capacity) {
    this->table = new FASTAreadset_HT(capacity);
}

BLAST_DB::~BLAST_DB() {
    delete table;
}

double BLAST_DB::readFile(Bacillus *bacillus) {
    return table->readFile(bacillus, 11);
}

int BLAST_DB::search(char *line, int fragmentLength) {
    int best = -fragmentLength;
    for (int i = 0; i < table->capacity; i++) {
        auto chain = table->chains[i];
        while (chain != nullptr) {
            auto s = nwScore(line, chain->line, fragmentLength, 11);
            best = max(s, best);
            chain = chain->next;
        }
    }
    return best;
}

int BLAST_DB::nwScore(char *a, char *b, int lengthA, int lengthB) {
    // constants
    int gapPenalty = -1;
    int mismatchPenalty = -1;
    int matchScore = 2;

    // init F matrix
    int **F = new int *[lengthA + 1];
    for (int i = 0; i <= lengthA; i++) F[i] = new int[lengthB + 1];

    for (int i = 0; i <= lengthA; i++) {
        F[i][0] = gapPenalty * i;
    }
    for (int j = 0; j <= lengthB; j++) {
        F[0][j] = gapPenalty * j;
    }
    for (int i = 1; i <= lengthA; i++) {
        for (int j = 1; j <= lengthB; j++) {
            char ac = a[i];
            char bc = b[j];

            // Match ← F(i−1, j−1) + S(Ai, Bj)
            int match = (ac == bc ? matchScore : mismatchPenalty) + F[i - 1][j - 1];
            // Delete ← F(i−1, j) + d
            int del = F[i - 1][j] + gapPenalty;
            // Insert ← F(i, j−1) + d
            int ins = F[i][j - 1] + gapPenalty;
            F[i][j] = max(match, max(del, ins));
        }
    }

    int lenMin = min(lengthA, lengthB);
    int s = F[lenMin][lenMin];

    if (false) { // debug
        for (int i = 0; i <= lengthA; i++) {
            for (int j = 0; j <= lengthB; j++) {
                int c = F[i][j];
                cout << (c >= 0 && c <= 9 ? " " : "") << c << " ";
            }
            cout << endl;
        }
        cout << endl;
        cout << endl;
    }

    // cleanup
    for (int i = 0; i < lengthA; i++) delete[] F[i];
    delete[] F;

    return s;
}

int main(int argc, char **argv) {
    if (argc < 4) {  // unexpected program call
        cout << endl << endl
             << "==========================" << endl
             << "INF 503 Homework 4 by Jadon Fowler" << endl
             << "Proper usage is:" << endl
             << "   ./homework4 problem1<A|B|C> <randomSampleSize> <input_file>" << endl
             << "Example:" << endl
             << "   ./homework4 problem1A 10000 /common/contrib/classroom/inf503/test_genome.fasta" << endl
             << "==========================" << endl << endl;
        exit(-1);
    }

    char *problem = argv[1];
    char *randomSizeString = argv[2];
    char *file_name = argv[3];
    cout << "Reading in data..." << endl;
    int randomFragmentAmount = atoi(randomSizeString);
    auto blast = new BLAST_DB(10000000);
    auto bacillus = new Bacillus(file_name);
    int fragmentLength = 50;
    double timeTaken = blast->readFile(bacillus);
    cout << "Done reading data in " << timeTaken << "s" << endl;

    if (strcmp("problem1A", problem) == 0) {
        clock_t timer = clock();
        int foundLines = 0;
        bacillus->randomize(5);
        cout << "Using " << randomFragmentAmount << " fragments from the Bacillus dataset as keys..." << endl;
        for (int i = 0; i < randomFragmentAmount; i++) {
            // find it in our table
            auto fragment = bacillus->randomStart(fragmentLength);
            int score = blast->search(fragment, fragmentLength);
            if (score >= 19/*magic?*/) {
                foundLines++;
                cout << foundLines << "/" << (i + 1) << endl;
            }
        }
        auto end = clock() - timer;
        auto firstSearchTime = ((double) end) / CLOCKS_PER_SEC;
        cout << "Found " << foundLines << " Bacillus lines in the dataset in " << firstSearchTime << "s."
             << endl;

        // search for 1 million random lines and search for them in the table
        timer = clock();
        foundLines = 0;
        cout << "Applying 5% error to Bacillus..." << endl;
        bacillus->randomize(5);

        cout << "Using " << randomFragmentAmount << " fragments from the Bacillus dataset as keys..." << endl;
        for (int i = 0; i < randomFragmentAmount; i++) {
            // find it in our table
            auto fragment = bacillus->randomStart(fragmentLength);
            int score = blast->search(fragment, fragmentLength);
            if (score >= 19/*magic?*/) {
                foundLines++;
                cout << foundLines << "/" << (i + 1) << endl;
            }
        }
        end = clock() - timer;
        auto secondSearchTime = ((double) end) / CLOCKS_PER_SEC;
        cout << "Found " << foundLines << " Bacillus lines in the dataset in " << secondSearchTime << "s."
             << endl;
    } else {
        cout << "Unknown problem: " << problem << endl;
    }
    cout << "Done!" << endl;
    delete blast;

    return 0;
}
