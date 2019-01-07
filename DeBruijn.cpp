#include "DeBruijn.h"
#include <iostream>

using namespace std;

DeBruijn::DeBruijn(const vector<string> &kmers, int k) {
    // construct a graph from kmers
    for (const auto & kmer : kmers) {
        auto prefix = kmer.substr(0, k - 1); // take first k-1
        auto prefixNode = getNode(prefix);

        if (prefixNode == nullptr) {
            prefixNode = addNode(prefix);
        }

        auto suffix = kmer.substr(1, k - 1); // take last k-1
        auto suffixNode = getNode(suffix);

        if (suffixNode == nullptr) {
            suffixNode = addNode(suffix);
        }

        connect(prefixNode, suffixNode);
    }
}

DeBruijn::Node* DeBruijn::addNode(const string &kmer) {
    unique_ptr<Node> node = make_unique<Node>();
    node->value = kmer;
    Node * node_raw = node.get();
    kmerNodeMap[kmer] = move(node);
    return node_raw;
}

DeBruijn::Node* DeBruijn::getNode(const string &kmer) {
    auto it = kmerNodeMap.find(kmer);

    if (it != kmerNodeMap.end()) {
        return it->second.get();
    } else {
        return nullptr;
    }
}

DeBruijn::Edge* DeBruijn::getEdge(DeBruijn::Node *n1, DeBruijn::Node *n2) {
    for (auto & edge : n1->outgoingEdges) {
        if (edge->next == n2) {
            return edge.get();
        }
    }
    return nullptr;
}

void DeBruijn::connect(DeBruijn::Node *n1, DeBruijn::Node *n2) {
    auto edge = getEdge(n1, n2);

    if (edge == nullptr) {
        auto edge = make_unique<Edge>();
        edge->next = n2;
        edge->weight = 1;

        n2->ingoingEdges.push_back(edge.get());
        n1->outgoingEdges.push_back(move(edge));
    } else {
        (edge->weight)++;
    }
}

void DeBruijn::printGraph() {
    for (auto & kmerNodePair : kmerNodeMap) {
        for (auto & out_edge : kmerNodePair.second->outgoingEdges) {
            cout
            << kmerNodePair.second->value
            << " -> "
            << out_edge->next->value
            << " "
            << out_edge->weight
            << endl;
        }
    }
}
