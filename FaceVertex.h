//
// Created by Xianyu on 2019-06-24.
//

#ifndef MESHSIMPLIFICATION_REFINED_FACEVERTEX_H
#define MESHSIMPLIFICATION_REFINED_FACEVERTEX_H
#include "Const.h"

class Face {
public:
    int ids[3];

    Face(int *v) {
        ids[0] = v[0], ids[1] = v[1], ids[2] = v[2];
        rerank();
    }

    Face(int v0, int v1, int v2) {
        ids[0] = v0, ids[1] = v1, ids[2] = v2;
        rerank();

    }
    void rerank(){
        // 使得保持顺时针的关系（保持法向量）
        // 使得最小值为ids0
        // 保证face唯一
        int min = std:: min(std:: min(ids[0], ids[1]), ids[2]);
        if(min == ids[0])
            return;
        else if(min == ids[1]) /* 012 -> 120 */
            std:: swap(ids[0], ids[1]), std:: swap(ids[1], ids[2]);
        else if(min == ids[2]) /* 012 -> 201 */
            std:: swap(ids[0], ids[1]), std:: swap(ids[0], ids[2]);
        return;
    }
        // 找到与v0，v1不同的值
    inline int find1Diff(int v0, int v1) {
        if (ids[0] != v0 && ids[0] != v1) return ids[0];
        if (ids[1] != v0 && ids[1] != v1) return ids[1];
        if (ids[2] != v0 && ids[2] != v1) return ids[2];
        assert(0);
    }
        //找到两个不同的值，并返回一个pair
    inline std::pair<int, int> find2Diff(int v0) {
        if (ids[0] == v0)
            return std::make_pair(ids[1], ids[2]);
        if (ids[1] == v0)
            return std::make_pair(ids[0], ids[2]);
        if (ids[2] == v0)
            return std::make_pair(ids[0], ids[1]);

        assert(0);
    }
        //检查是否有重复的边
    inline bool checkIfDuplicate(){
        if(ids[0] == ids[1] || ids[1] == ids[2] || ids[0] == ids[2])
            return true;
        else
            return false;
    }
        //判断是否包含2个
    inline bool ifcontain2(int v0, int v1) {

        assert(v0!=v1);
        if ((ids[0] == v0 || ids[1] == v0 || ids[2] == v0) && (ids[0] == v1 || ids[1] == v1 || ids[2] == v1))
            return true;
        return false;
    }
        //判断是否包含1个
    inline bool ifcontain1(int v) {
        if (ids[0] == v || ids[1] == v || ids[2] == v)
            return true;
        return false;
    }

    /*Find a way to avoid duplicate in set*/
    /*Reload the function of operator < */
    friend bool operator<(const Face &a, const Face &b) {
        if( a.ids[0] < b.ids[0] || (a.ids[0] == b.ids[0] && a.ids[1] < b.ids[1]) || (a.ids[0] == b.ids[0] && a.ids[1] == b.ids[1] && a.ids[2] < b.ids[2]))
            return true;
        return false;
    }

    friend std::ostream &operator<<(std::ostream &os, const Face &face) {
        os << "f " << face.ids[0] << " " << face.ids[1] << " " << face.ids[2];
        return os;
    }
    bool operator==(const Face& right) const{
        if(this->ids[0] == right.ids[0] || this->ids[0] == right.ids[1] || this->ids[0] == right.ids[2])
            if(this->ids[1] == right.ids[0] || this->ids[1] == right.ids[1] || this->ids[1] == right.ids[2])
                if(this->ids[2] == right.ids[0] || this->ids[2] == right.ids[1] || this->ids[2] == right.ids[2])
                    return true;
        return false;
    }

};


class Vertex {
public:
    explicit Vertex(int _id) {
        id = _id;
        memset(dim, 0, sizeof(double)* 3);
        memset(q, 0, sizeof(double)*16);
    }

    int id;
    double dim[3], q[16];
    std::set<Face> faces;

    inline friend std::ostream &operator<<(std::ostream &os, const Vertex &v) {
        os << "v " << v.dim[0] << " " << v.dim[1] << " " << v.dim[2];
        return os;
    }
    inline void addFace(const Face face) {
        faces.insert(face);
    }
    inline Face removeFaceWithReturn(const Face &face) {

        /*Linear complexity*/
        /*TODO: Try to use the default set remove and find api*/
        for(auto i = faces.begin(); i != faces.end(); i++)
            if(i->ids[0] == face.ids[0] || i->ids[0] == face.ids[1] || i->ids[0] == face.ids[2])
                if(i->ids[1] == face.ids[0] || i->ids[1] == face.ids[1] || i->ids[1] == face.ids[2])
                    if(i->ids[2] == face.ids[0] || i->ids[2] == face.ids[1] || i->ids[2] == face.ids[2]){
                        Face returnface(i->ids[0], i->ids[1], i->ids[2]);

                        faces.erase(i);
                        return returnface;
                    }
        return Face(0,0,0);

    }

    inline void removeFaceDirectly(const Face& face){
        faces.erase(face);

    }

    inline double distance(const Vertex &v) {
        return sqrt(SQR(dim[0] - v.dim[0]) + SQR(dim[1] - v.dim[1]) + SQR(dim[2] - v.dim[2]));
    }

    /*find face containing v0 and v1*/
    inline std::vector<Face> findFaces(int v0) {
        std::vector<Face> findings;
        for (auto face:faces) {
            if (face.ifcontain1(v0))
                findings.push_back(face);
        }
        return findings;
    }

    inline void replace(double *nd, double *nq) {
        for (int i = 0; i < 3; ++i)
            dim[i] = nd[i];

        for (int i = 0; i < 16; ++i)
            q[i] += nq[i];
        return;
    }

};


#endif //MESHSIMPLIFICATION_REFINED_FACEVERTEX_H
