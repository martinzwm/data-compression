#include "bits.h"
#include "treenode.h"
#include "huffman.h"
#include "map.h"
#include "vector.h"
#include "priorityqueue.h"
#include "strlib.h"
#include "testing/SimpleTest.h"
using namespace std;


void deallocateTree(EncodingTreeNode* tree)
{
    // base case
    if (tree == nullptr) {
        return;
    }

    // recursive case
    deallocateTree(tree->zero);
    deallocateTree(tree->one);
    delete tree;
}

/**
 * Given a Queue<Bit> containing the compressed message bits and the encoding tree
 * used to encode those bits, decode the bits back to the original message text.
 *
 * You can assume that tree is a well-formed non-empty encoding tree and
 * bits queue contains a valid sequence of encoded bits.
 *
 */
bool isLeaf(EncodingTreeNode*& node);

string decodeText(EncodingTreeNode* tree, Queue<Bit>& messageBits)
{
    string message = "";
    EncodingTreeNode* curr = tree;

    while (!messageBits.isEmpty()) {
        Bit op = messageBits.dequeue();
        // step left or right depending on the bit
        if (op == 0) {
            curr = curr->zero;
        } else {
            curr = curr->one;
        }

        // once we reached to a character in the coding tree
        if (isLeaf(curr)) {
            message += curr->ch; // update message
            curr = tree; // reset curr to the root of the encoding tree
        }
    }
    return message;
}

/**
 * Reconstruct encoding tree from flattened form Queue<Bit> and Queue<char>.
 *
 * You can assume that the input Queues are well-formed and represent
 * a valid encoding tree.
 */
// Note that this function is deconstructive on treeBits and treeLeaves, that is okay
// because in reality we only need to read treeBits and treeLeaves once.
EncodingTreeNode* unflattenTree(Queue<Bit>& treeBits, Queue<char>& treeLeaves)
{
    EncodingTreeNode* root = new EncodingTreeNode(' ');

    // Sanity check (we are done)
    // Note that for legal inputs, treeBits should not be exhausted before treeLeaves
    // so we only need to check for treeBits.
    if (treeBits.isEmpty()) {
        return root;
    }

    Bit bit = treeBits.dequeue();
    // base case
    if (bit == 0) {
        root->ch = treeLeaves.dequeue();
        return root;
    }

    // recursive cases
    root->zero = unflattenTree(treeBits, treeLeaves);
    root->one = unflattenTree(treeBits, treeLeaves);
    return root;
}

/**
 * Decompress the given EncodedData and return the original text.
 *
 * You can assume the input data is well-formed and was created by a correct
 * implementation of compress.
 *
 * Your implementation may change the data parameter however you like. There
 * are no requirements about what it should look like after this function
 * returns.
 */
string decompress(EncodedData& data)
{
    EncodingTreeNode* coding_tree = unflattenTree(data.treeBits, data.treeLeaves);
    string message = decodeText(coding_tree, data.messageBits);

    deallocateTree(coding_tree);
    return message;
}


/**
 * Constructs an optimal Huffman coding tree for the given text, using
 * the algorithm described in lecture.
 *
 * Reports an error if the input text does not contain at least
 * two distinct characters.
 *
 * When assembling larger trees out of smaller ones, make sure to set the first
 * tree dequeued from the queue to be the zero subtree of the new tree and the
 * second tree as the one subtree.
 */

PriorityQueue<EncodingTreeNode*> getFrequency(string text);
EncodingTreeNode* convert_to_tree(PriorityQueue<EncodingTreeNode*>& freq);

EncodingTreeNode* buildHuffmanTree(string text)
{
    int N = text.length();
    if (N < 2) {
        error("The text to compress should contain at least 2 words.");
    }

    PriorityQueue<EncodingTreeNode*> freq = getFrequency(text);
    return convert_to_tree(freq);
}

// Helper function to buildHuffmanTree()
EncodingTreeNode* convert_to_tree(PriorityQueue<EncodingTreeNode*>& freq) {
    while (freq.size() > 1) {
        int left_num = freq.peekPriority();
        EncodingTreeNode* left = freq.dequeue();
        int right_num = freq.peekPriority();
        EncodingTreeNode* right = freq.dequeue();

        EncodingTreeNode* combined = new EncodingTreeNode(left, right);
        freq.enqueue(combined, left_num+right_num);
    }
    return freq.dequeue();
}

// Helper function to buildHuffmanTree()
PriorityQueue<EncodingTreeNode*> getFrequency(string text) {
    // Create frequency data using an hashmap O(n)
    unordered_map<char, int> freq;
    for (char character : text) {
        freq[character] += 1;
    }

    // Organize frequency data based on their frequency O(nlogn)
    PriorityQueue<EncodingTreeNode*> pq;
    for (auto pair : freq) {
        EncodingTreeNode* dummy = new EncodingTreeNode(pair.first);
        pq.add(dummy, pair.second);
    }
    return pq;
}

/**
 * Given a string and an encoding tree, encode the text using the tree
 * and return a Queue<Bit> of the encoded bit sequence.
 *
 * You can assume tree is a valid non-empty encoding tree and contains an
 * encoding for every character in the text.
 */
void createDict(EncodingTreeNode*& node, Vector<Bit> curr, unordered_map<char, Vector<Bit>>& encode_dict);

Queue<Bit> encodeText(EncodingTreeNode* tree, string text)
{
    unordered_map<char, Vector<Bit>> encode_dict;
    Vector<Bit> curr;
    createDict(tree, curr, encode_dict);

    Queue<Bit> code;
    for (char ch : text) {
        if (encode_dict.find(ch)==encode_dict.end()) {
            error("Character to decipher doesn't exist in the dictionary.");
        }
        for (Bit bit : encode_dict[ch]) {
            code.enqueue(bit);
        }
    }

    return code;
}

// Helper function to encodeText.
// Populate encode_dict based on the encoding tree.
void createDict(EncodingTreeNode*& node, Vector<Bit> curr, unordered_map<char, Vector<Bit>>& encode_dict) {
    if (isLeaf(node)) {
        encode_dict.insert({node->ch, curr});
        return;
    }

    curr.add(Bit(0));
    createDict(node->zero, curr, encode_dict);
    curr.remove(curr.size()-1); // undo

    curr.add(Bit(1));
    createDict(node->one, curr, encode_dict);
    curr.remove(curr.size()-1); // undo
}

/**
 * Flatten the given tree into a Queue<Bit> and Queue<char> in the manner
 * specified in the assignment writeup.
 *
 * You can assume the input Queues are empty on entry to this function.
 *
 * You can assume tree is a valid well-formed encoding tree.
 */
void flattenTree(EncodingTreeNode* tree, Queue<Bit>& treeBits, Queue<char>& treeLeaves)
{
    // base case (leaf)
    if (isLeaf(tree)) {
        treeLeaves.enqueue(tree->ch);
        treeBits.enqueue(Bit(0));
        return;
    }

    // recursive cases (non-leaf)
    treeBits.enqueue(Bit(1));
    flattenTree(tree->zero, treeBits, treeLeaves);
    flattenTree(tree->one, treeBits, treeLeaves);
}

/**
 * Compress the input text using Huffman coding, producing as output
 * an EncodedData containing the encoded message and encoding tree used.
 *
 * Reports an error if the message text does not contain at least
 * two distinct characters.
 */
EncodedData compress(string messageText)
{
    EncodedData file;
    EncodingTreeNode* coding_tree = buildHuffmanTree(messageText);
    file.messageBits = encodeText(coding_tree, messageText);

    Queue<Bit>  treeBits;
    Queue<char> treeLeaves;
    flattenTree(coding_tree, treeBits, treeLeaves);
    file.treeBits = treeBits;
    file.treeLeaves = treeLeaves;

    deallocateTree(coding_tree);
    return file;
}

/* * * * * * Test Cases Below This Point * * * * * */

// TODO: Write your own student tests



/* * * * * Provided Tests Below This Point * * * * */

bool isLeaf(EncodingTreeNode*& node) {
    if (node == nullptr) {
        error("should never reach here.");
    }
    return node->zero == nullptr && node->one == nullptr;
}

bool areEqual(EncodingTreeNode* a, EncodingTreeNode* b) {
    // base case
    if (a == nullptr && b == nullptr) {
        return true;
    } else if (a == nullptr || b == nullptr) {
        return false;
    }
    if (isLeaf(a) && isLeaf(b)) {
        if (a->ch != b->ch) {
            return false;
        }
    }

    // recursive case
    return areEqual(a->zero, b->zero) && areEqual(a->one, b->one);
}

EncodingTreeNode* createExampleTree()
{
    /* Example encoding tree used in multiple test cases:
     *                *
     *              /   \
     *             T     *
     *                  / \
     *                 *   E
     *                / \
     *               R   S
     */
    EncodingTreeNode* T = new EncodingTreeNode('T');
    EncodingTreeNode* R = new EncodingTreeNode('R');
    EncodingTreeNode* S = new EncodingTreeNode('S');
    EncodingTreeNode* E = new EncodingTreeNode('E');

    EncodingTreeNode* right_left = new EncodingTreeNode(R, S);
    EncodingTreeNode* right = new EncodingTreeNode(right_left, E);
    EncodingTreeNode* root = new EncodingTreeNode(T, right);

    return root;
}

STUDENT_TEST("Test functions in warmup") {
    EncodingTreeNode* tree1 = createExampleTree();
    EncodingTreeNode* tree2 = createExampleTree();
    EXPECT(areEqual(tree1, tree2));

    EXPECT(areEqual(tree1->one, tree2->one));

    EXPECT(!areEqual(tree1->zero, tree2->one));

    // After the below two lines, there should be no memory leak.
    deallocateTree(tree1);
    deallocateTree(tree2);
}

PROVIDED_TEST("decodeText, small example encoding tree") {
    EncodingTreeNode* tree = createExampleTree(); // see diagram above

    Queue<Bit> messageBits = { 1, 1 }; // E
    EXPECT_EQUAL(decodeText(tree, messageBits), "E");

    messageBits = { 1, 0, 1, 1, 1, 0 }; // SET
    EXPECT_EQUAL(decodeText(tree, messageBits), "SET");

    messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1}; // STREETS
    EXPECT_EQUAL(decodeText(tree, messageBits), "STREETS");

    deallocateTree(tree);
}

PROVIDED_TEST("unflattenTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  treeBits   = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'T', 'R', 'S', 'E' };
    EncodingTreeNode* tree = unflattenTree(treeBits, treeLeaves);

    EXPECT(areEqual(tree, reference));

    deallocateTree(tree);
    deallocateTree(reference);
}

PROVIDED_TEST("decompress, small example input") {
    EncodedData data = {
        { 1, 0, 1, 1, 0, 0, 0 }, // treeBits
        { 'T', 'R', 'S', 'E' },  // treeLeaves
        { 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1 } // messageBits
    };

    EXPECT_EQUAL(decompress(data), "TRESS");
}

PROVIDED_TEST("buildHuffmanTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    EncodingTreeNode* tree = buildHuffmanTree("STREETTEST");
    EXPECT(areEqual(tree, reference));

    deallocateTree(reference);
    deallocateTree(tree);
}

PROVIDED_TEST("encodeText, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above

    Queue<Bit> messageBits = { 1, 1 }; // E
    EXPECT_EQUAL(encodeText(reference, "E"), messageBits);

    messageBits = { 1, 0, 1, 1, 1, 0 }; // SET
    EXPECT_EQUAL(encodeText(reference, "SET"), messageBits);

    messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1 }; // STREETS
    EXPECT_EQUAL(encodeText(reference, "STREETS"), messageBits);

    deallocateTree(reference);
}

PROVIDED_TEST("flattenTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  expectedBits   = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> expectedLeaves = { 'T', 'R', 'S', 'E' };

    Queue<Bit>  treeBits;
    Queue<char> treeLeaves;
    flattenTree(reference, treeBits, treeLeaves);

    EXPECT_EQUAL(treeBits,   expectedBits);
    EXPECT_EQUAL(treeLeaves, expectedLeaves);

    deallocateTree(reference);
}

PROVIDED_TEST("compress, small example input") {
    EncodedData data = compress("STREETTEST");
    Queue<Bit>  treeBits    = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> treeChars   = { 'T', 'R', 'S', 'E' };
    Queue<Bit>  messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0 };

    EXPECT_EQUAL(data.treeBits, treeBits);
    EXPECT_EQUAL(data.treeLeaves, treeChars);
    EXPECT_EQUAL(data.messageBits, messageBits);
}

PROVIDED_TEST("Test end-to-end compress -> decompress") {
    Vector<string> inputs = {
        "HAPPY HIP HOP",
        "The job requires extra pluck and zeal from every young wage earner.",
        ":-) :-D XD <(^_^)>",
    };

    for (string input: inputs) {
        EncodedData data = compress(input);
        string output = decompress(data);

        EXPECT_EQUAL(output.size(), input.size());

        /* Don't clobber the output with a huge string if there's a mismatch. */
        EXPECT(input == output);
    }
}
