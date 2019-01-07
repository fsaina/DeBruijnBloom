#ifndef DEBRUIJNBLOOM_DeBruijn_H
#define DEBRUIJNBLOOM_DeBruijn_H

#include <unordered_map>
#include <string>
#include <vector>

using namespace std;

class DeBruijn {

private:
    struct Edge;

    struct Node {
        vector<Edge *> ingoingEdges;
        vector<unique_ptr<Edge>> outgoingEdges;
        string value;
    };

    struct Edge {
        Node *next;
        int weight;
    };

    unordered_map<string, unique_ptr<Node>> kmerNodeMap;

public:
    DeBruijn(const vector<string> &kmers, int k);
    Node *addNode(const string &kmer);
    Edge *getEdge(Node *n1, Node *n2);
    Node *getNode(const string &kmer);
    void connect(Node *n1, Node *n2);
    void print();
};

#endif //DEBRUIJNBLOOM_DeBruijn_H
