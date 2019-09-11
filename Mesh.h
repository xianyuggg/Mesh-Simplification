//
// Created by Xianyu on 2019-06-24.
//

#ifndef MESHSIMPLIFICATION_REFINED_MESH_H
#define MESHSIMPLIFICATION_REFINED_MESH_H


#include "Const.h"
#include "FaceVertex.h"
#include "PairHeap.h"

class Mesh {
public:
    int fcnt = 0, vcnt = 0;
    std::vector<Vertex> vertexes;
    std::vector<std::set<int> > edges;
    Heap heap;

    void calculateQ();
    // P stands for [a, b, c, d]T
    // Kp P dot PT ~ 4*1 dot 1*4: 4*4
    void addKp(const Face &face, double *kp);


    // Using 3 vertexes' coordinates to calculate P;
    void getP(const Face &face, double *p);


    void selectPairs();

    void addPair(int i, int j);
    void calculateV(double* q, double* nd, double* ad, double* bd);
    double calculateCost(double* dim, double* cost);

    void aggregation();


    Mesh() {
        vertexes.clear();
        edges.clear();
    };
    ~Mesh() = default;


    bool loadObj(std::string path);

    void writeObj(std::string path);


};

#endif //MESHSIMPLIFICATION_REFINED_MESH_H
