//
// Created by Xianyu on 2019-06-12.
//

#include "Mesh.h"
  

void Mesh::aggregation() {

    Pair pair = heap.top();
    heap.pop();
    int v0 = pair.v0, v1 = pair.v1;

    /* Print warning info */
    if (vertexes[v1].faces.size() == 0 || v0 == v1 || vertexes[v0].id == -1 || vertexes[v1].id == -1) {
        std::cout << "size0, V0 == V1, id-1 Warning!" << std::endl;
        return;
    }

    /*Remove v1 */
    vertexes[v1].id = -1;

    /* Find faces containing v0 v1 from v1
     * And remove these faces
     * */
    std::vector<Face> findings = vertexes[v1].findFaces(v0);
    for (auto &face:findings) {
        int v = face.find1Diff(v0, v1);
        vertexes[v].removeFaceDirectly(face);
        vertexes[v0].removeFaceDirectly(face);
        vertexes[v1].removeFaceDirectly(face);
        --fcnt;
    }

    /* Connect the connected points of v1 to v0 */
    for (Face face:vertexes[v1].faces) {
        /* The faces containing v0 and v1 have been moved
         * i.e. Update faces with v1;
         * */

        /* 找到一个三角形面片的另外两个点 */
        std::pair<int, int> vpair = face.find2Diff(v1);
        vertexes[vpair.first].removeFaceDirectly(face);
        vertexes[vpair.second].removeFaceDirectly(face);
        /* 将v1修改为v0 */
        for(int i = 0; i < 3; ++i)
            if(face.ids[i] == v1) {
                face.ids[i] = v0;
            }

        vertexes[v0].addFace(face);
        vertexes[vpair.first].addFace(face);
        vertexes[vpair.second].addFace(face);

    }
    /*delete old face in v1*/
    vertexes[v1].faces.clear();
    /*Replace with new node*/
    vertexes[v0].replace(pair.dim, vertexes[v1].q);
    /*id may change, so we don't want to delete q0*/
    vcnt--;
    /*Remove all the pairs containing q1, and q0, and reheap */


    /* Remove all the pairs containing v1 in the queue, and connect them to v0 */
    std:: set<int> pairs0 = heap.pairs[v0];
    std:: set<int> pairs1 = heap.pairs[v1];
    for(auto v: pairs1) {
        heap.remove(Pair(v, v1));
        if(!pairs0.count(v) && v != v0)
            addPair(v, v0);
    }
    /* Update the value containing v0 in the queue. */
    for(auto v: pairs0) if(v != v1)
            addPair(v, v0);

    return;

}
bool Mesh::loadObj(std::string path) {

    FILE *filein = fopen(path.c_str(), "r");
    if(filein == NULL)
    {
        std::cout << "Can't open the input file." << std::endl;
        return false;
    }
    std::cout << "Reading obj file " + path + "..." << std::endl;
    char info;
    char cmtLine[1000];

    while(fscanf(filein, "%c", &info) != EOF)
    {
        if(info == 'v')
        {
            Vertex vertex(vcnt++);
            fscanf(filein, "%lf%lf%lf", &vertex.dim[0], &vertex.dim[1], &vertex.dim[2]);
            vertexes.push_back(vertex);
        }
        else if(info == 'f')
        {
            int index[3];
            ++fcnt;
            fscanf(filein, "%d %d %d", &index[0], &index[1], &index[2]);
            index[0]--,index[1]--,index[2]--;
            Face face(index);
            for (auto i:index)
                vertexes[i].addFace(face);
        }
        else if(info == '#')
            fgets(cmtLine, 1000, filein);
    }
    fclose(filein);
    std::cout << "Vertex count : " << vcnt << std::endl;
    std::cout << "Face count : " << fcnt << std::endl;
    this->heap.pairs.resize(vcnt);
    return true;

}


void Mesh::writeObj(std::string path) {


    std::cout << "Writing obj file " + path + "..." << std::endl;
    std::ofstream output(path);
    int index = 1;
    double error = 0;

    /*Assign new id for each vertex*/
    for (int i = 0; i < vertexes.size(); ++i)
        if (vertexes[i].id != -1) {
            vertexes[i].id = index++;
            error += calculateCost(vertexes[i].dim, vertexes[i].q);
            output << vertexes[i] << std::endl;
        }

    std::set<Face> faceset;
    for (auto &vertex:vertexes) {
        for (auto face:vertex.faces) {
            bool insert = true;
            for (int i = 0; i < 3; ++i) {
                face.ids[i] = vertexes[face.ids[i]].id;
                if (face.ids[i] == -1) {
                    insert = false;
                    break;
                }
            }
            if (insert)
                faceset.insert(face);
        }
    }
    for (auto &face:faceset) {
        output << face << std::endl;
    }
    std::cout << "Avg Error: " << error / (index-1) << std::endl;
    std::cout << "Vertex count: " << vcnt << std::endl;
    std::cout << "Face count: " << faceset.size() << std::endl;

}





void Mesh::calculateQ() {
    std::cout << "Calculating Q for each vertex..." << std::endl;
    for (auto &v: vertexes) {
        memset(v.q, 0, sizeof(double) * 16);
        /* Q equals to SUM of Kp(P*Pt) , P = [a,b,c,d]T */
        for (auto &face: v.faces)
            addKp(face, v.q);
    }
    std::cout << "Calculating Q complete!" << std::endl;

}

void Mesh::addKp(const Face &face, double *kp) {

    double p[4];
    getP(face, p);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            kp[i * 4 + j] += p[i] * p[j];

}

void Mesh::getP(const Face &face, double *p) {
    /*
     * A = (y3 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1);
     * B = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
     * C = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
     */
    /* Ax + Bx + Cx + D = 0 */
    /* A,B,C,D should be normalized */
    double *p1 = vertexes[face.ids[0]].dim; // 0 ~ x, 1 ~ y, 2 ~ z;
    double *p2 = vertexes[face.ids[1]].dim;
    double *p3 = vertexes[face.ids[2]].dim;

    p[0] = (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p3[1] - p1[1]) * (p2[2] - p1[2]);
    p[1] = (p2[2] - p1[2]) * (p3[0] - p1[0]) - (p3[2] - p1[2]) * (p2[0] - p1[0]);
    p[2] = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
    double div = sqrt(SQR(p[0]) + SQR(p[1]) + SQR(p[2]));
    p[0] /= div, p[1] /= div, p[2] /= div;
    p[3] = -p[0] * p1[0] - p[1] * p1[1] - p[2] * p1[2];

}

void Mesh::addPair(int i, int j) {
    assert(i != j);
    Pair pair(i, j);
    double *a = vertexes[i].dim;
    double *b = vertexes[j].dim;
    double *aq = vertexes[i].q;
    double *bq = vertexes[j].q;

    double q[16];
    for (int i = 0; i < 16; ++i)
        q[i] = aq[i] + bq[i];

    calculateV(q, pair.dim, a, b);
    pair.cost = calculateCost(pair.dim, q);
    heap.insert(pair);
}


/*
 * new q, new dim, a dim, b dim
 */
void Mesh::calculateV(double *q, double *nd, double *ad, double *bd) {
    double det =
            q[0] * q[5] * q[10] - q[0] * q[6] * q[6] - q[1] * q[1] * q[10] + q[1] * q[6] * q[2] + q[2] * q[1] * q[6] -
            q[2] * q[5] * q[2];
    /*
     * If the matrix is not invertible
     * we use the mid-point
     */
    if (ABS(det) <= EPSILON) {
        nd[0] = (ad[0] + bd[0]) / 2., nd[1] = (ad[1] + bd[1]) / 2., nd[2] = (ad[2] + bd[2]) / 2.;
    } else {
        double x = q[3] * q[5] * q[10] - q[3] * q[6] * q[6] - q[7] * q[1] * q[10] + q[7] * q[6] * q[2] +
                   q[11] * q[1] * q[6] - q[11] * q[5] * q[2];
        double y = q[0] * q[7] * q[10] - q[0] * q[11] * q[6] - q[1] * q[3] * q[10] + q[1] * q[11] * q[2] +
                   q[2] * q[3] * q[6] - q[2] * q[7] * q[2];
        double z = q[0] * q[5] * q[11] - q[0] * q[6] * q[7] - q[1] * q[1] * q[11] + q[1] * q[6] * q[3] +
                   q[2] * q[1] * q[7] - q[2] * q[5] * q[3];
        nd[0] = -x / det, nd[1] = -y / det, nd[2] = -z / det;
    }
    return;
}

double Mesh::calculateCost(double *dim, double *q) {
    double tmp[4];
    memset(tmp, 0, sizeof(double) * 4);
    /* Expand the original dim */
    /* [Vx,Vy,Vz,1] */
    double vector[4] = {dim[0], dim[1], dim[2], 1};

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            tmp[i] += vector[j] * q[j * 4 + i];
    double cost = 0;
    for (int i = 0; i < 4; ++i)
        cost += (tmp[i] * vector[i]);
    return cost;

}

void Mesh::selectPairs() {
    std::cout << "Selecting Pairs ~ t = " << THRESHOLD << "..." << std::endl;
    edges.resize(vertexes.size());

    for (int i = 0; i < vertexes.size(); ++i)
        for (auto face: vertexes[i].faces)
            for (int j = 0; j < 3; ++j)
                if (face.ids[j] != i)
                    edges[i].insert(face.ids[j]);

    /*
     * 排序用于剪枝
     */
    std:: vector<int> index;
    index.resize(vertexes.size());
    for(int i = 0; i < index.size(); ++ i) index[i] = i;
    auto cmp = [this] (int i, int j) -> bool { return vertexes[i].dim[0] < vertexes[j].dim[0]; };
    std:: sort(index.begin(), index.end(), cmp);

    /* t equals to THRESHOLD */
    /* Adding all edges to pair */
    for(int i = 0; i < vertexes.size(); ++ i)
        for(int j: edges[i]) if(i < j)
                addPair(i, j);
    /*adding all pairs smaller than THRESHOLD or Edge(a,b) exists*/
    for (int i = 0; i < vertexes.size(); ++i) {
        int a = index[i], b;
        for (int j = i + 1; j < vertexes.size(); ++j) {
            b = index[j];
            if (vertexes[a].distance(vertexes[b]) < THRESHOLD)
                addPair(a, b);
            if(vertexes[b].dim[0] - vertexes[a].dim[0] > THRESHOLD)
                break;
        }

    }
    std::cout << "Pairs size: " << heap.queue.size() << std::endl;

}
