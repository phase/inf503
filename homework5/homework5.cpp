//
// INF 503 - Homework 5 + 6
// Jadon Fowler
//

#include "homework5.h"

using namespace std;

// how many characters can be different for a search to still return true
const int FUZZY_MAX = 1;

// limit for reading lines - good for testing
// setting this to -1 will remove the limit
constexpr int MAX_LINES = -1;

int indexFrom(char c) {
    switch (c) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'T':
            return 2;
        case 'G':
            return 3;
    }
    cout << "tried to convert " << c << " to index???? err" << endl;
    return -1;
}


Dataset::Dataset(char *fileName) {
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
    cout << "Loaded dataset " << fileName << " in " << timeTaken << " seconds." << endl;
}

char *Dataset::randomStart(int size) {
    int offset = (rand() % (length - size));
    return data + offset;
}

void Dataset::randomize(int percent) {
    for (int i = 0; i < length; i++) {
        int p = rand() % 100;
        if (p < percent) {
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
                case 3: // fallthrough
                default:
                    data[i] = 'G';
                    break;
            }
        }
    }
}

Dataset::~Dataset() {
    delete data;
}

trie_node::trie_node(char v) {
    this->value = v;
    this->terminal = false;
    for (int i = 0; i < keys; i++) {
        children[i] = nullptr;
    }
}

trie_node::~trie_node() {
    for (int i = 0; i < keys; i++) {
        auto ptr = children[i];
        if (ptr != nullptr) {
            delete ptr;
        }
    }
}

void trie_node::insert(char *line, int length, int offset) {
    auto ch = line[offset];
    if (ch == '\0') {
        // we've reached the end of the line
        return;
    }
    auto childIndex = indexFrom(ch);
    auto child = children[childIndex];
    if (child == nullptr) {
        // allocate a new node for this path
        child = new trie_node(ch);
        children[childIndex] = child;
    }
    if (offset + 1 == length) {
        // last one
        child->terminal = true;
    } else {
        // continue down the string
        child->insert(line, length, offset + 1);
    }
}

bool trie_node::search(char *line, int length, int offset, int fuzz) {
    if (offset >= length) {
        return false;
    }
    auto ch = line[offset];
    if (ch == '\0') {
        return false;
    }
    auto childIndex = indexFrom(ch);
    auto child = children[childIndex];
    if (child == nullptr) {
        // do fuzzy
        if (fuzz < FUZZY_MAX) {
            // check every other path
            for (int i = 0; i < keys; i++) {
                if (i == childIndex) {
                    // don't recheck the proper path
                    continue;
                }
                auto possiblePath = children[i];
                if (possiblePath != nullptr) {
                    // search other paths, only return early if we found it
                    if (possiblePath->search(line, length, offset + 1, fuzz + 1)) {
                        return true;
                    }
                }
            }
        }
        return false;
    } else {
        if (offset + 1 == length) {
            // last one
            return child->terminal;
        }
        return child->search(line, length, offset + 1, fuzz + 1);
    }
}

int trie_node::countNodes() {
    int total = 1;
    for (int i = 0; i < keys; i++) {
        auto ptr = children[i];
        if (ptr != nullptr) {
            total += ptr->countNodes();
        }
    }
    return total;
}

int trie_node::countTerminalNodes() {
    int total = this->terminal ? 1 : 0;
    for (int i = 0; i < keys; i++) {
        auto ptr = children[i];
        if (ptr != nullptr) {
            total += ptr->countTerminalNodes();
        }
    }
    return total;
}

prefix_trie::prefix_trie() {
    root = new trie_node('\0');
}

prefix_trie::~prefix_trie() {
    delete root;
}

void prefix_trie::insert(char *line, int length) {
    root->insert(line, length, 0);
}

bool prefix_trie::search(char *line, int length) {
    return root->search(line, length, 0, 0);
}

int prefix_trie::countNodes() {
    return root->countNodes();
}

suffix_tree::suffix_tree() {
    root = new trie_node('\0');
}

suffix_tree::suffix_tree(Dataset *dataset, int fragmentLength) {
    root = new trie_node('\0');
    int totalSubstrings = dataset->length - fragmentLength;
    for (int i = 0; i < totalSubstrings; i++) {
        // start at the ith character in the dataset
        char *currentFragment = dataset->data + i;
        // insert it in our table
        insert(currentFragment, fragmentLength);
    }
}

suffix_tree::~suffix_tree() {
    delete root;
}

void suffix_tree::insert(char *line, int length) {
    for (int i = 0; i < length; i++) {
        root->insert(line + i, length - i, 0);
    }
}

bool suffix_tree::search(char *line, int length) {
    return root->search(line, length, 0, 10 /* disable fuzzy matching */);
}


char *generateGenomeLine(int length) {
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

int main(int argc, char **argv) {
    srand(time(NULL)); // seed
    if (argc < 4) {  // unexpected program call
        cout << endl << endl
             << "==========================" << endl
             << "INF 503 Homework 5 (& 6) by Jadon Fowler" << endl
             << "Proper usage is:" << endl
             << "   ./homework5 <5|6><A|B> <capacity> <input_file>" << endl
             << "Example:" << endl
             << "   ./homework5 5A 10000 /common/contrib/classroom/inf503/test_genome.fasta" << endl
             << "==========================" << endl << endl;
        exit(-1);
    }

    char *problem = argv[1];
    char *capacity_string = argv[2];
    char *file_name = argv[3];
    cout << "Reading in SARS-COV2 data..." << endl;
    int capacity = atoi(capacity_string);
    auto dataset = new Dataset(file_name);

    // B is the same as A with randomization
    auto randomizedDataset = new Dataset(file_name);
    if (strcmp("5B", problem) == 0) {
        cout << "Randomizing dataset..." << endl;
        randomizedDataset->randomize(5);
    }

    int fragmentLength = 36;

    if (strcmp("5A", problem) == 0 || strcmp("5B", problem) == 0) {
        cout << "Building trie from SARS-COV2..." << endl;
        auto trie = new prefix_trie();
        for (int i = 0; i < capacity; i++) {
            auto fragment = randomizedDataset->randomStart(fragmentLength);
            trie->insert(fragment, fragmentLength);
        }
        auto nodeCount = trie->countNodes();
        cout << "The trie contains " << nodeCount << " nodes." << endl;
        cout << "The trie contains " << trie->root->countTerminalNodes() << " terminal nodes." << endl;

        int totalSubstrings = dataset->length - fragmentLength;

        clock_t timer = clock();
        int foundLines = 0;
        cout << "Using " << totalSubstrings << " fragments from the SARS-COV2 dataset as keys..." << endl;
        for (int i = 0; i < totalSubstrings; i++) {
            // start at the ith character in the dataset
            char *currentFragment = dataset->data + i;
            // find it in our table
            if (trie->search(currentFragment, fragmentLength)) {
                foundLines++;
            }
        }
        auto end = clock() - timer;
        auto firstSearchTime = ((double) end) / CLOCKS_PER_SEC;
        cout << "Found " << foundLines << " SARS-COV2 lines in the dataset in " << firstSearchTime << "s." << endl;
    } else if (strcmp("6A", problem) == 0) {
        cout << "Building suffix tree from SARS-COV2..." << endl;
        int totalSubstrings = dataset->length - fragmentLength;
        cout << "Inserting " << totalSubstrings << " fragments from the SARS-COV2 dataset into the tree..." << endl;

        clock_t timer = clock();
        auto tree = new suffix_tree(dataset, fragmentLength);
        auto end = clock() - timer;
        auto firstSearchTime = ((double) end) / CLOCKS_PER_SEC;
        cout << "Loaded SARS-COV2 data in " << firstSearchTime << "s." << endl;
        cout << "The suffix tree contains " << tree->root->countNodes() << " nodes." << endl;
        cout << "The suffix tree contains " << tree->root->countTerminalNodes() << " terminal nodes." << endl;

        timer = clock();
        int foundLines = 0;
        for (int i = 0; i < capacity; i++) {
            auto fragment = randomizedDataset->randomStart(fragmentLength);
            if (tree->search(fragment, fragmentLength)) {
                foundLines++;
            }
        }
        end = clock() - timer;
        auto secondSearchTime = ((double) end) / CLOCKS_PER_SEC;
        cout << "Found " << foundLines << " SARS-COV2 lines in the dataset in " << secondSearchTime << "s." << endl;
    } else {
        cout << "Unknown problem: " << problem << endl;
    }
    cout << "Done!" << endl;

    return 0;
}
