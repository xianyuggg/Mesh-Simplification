//
// Created by Xianyu on 2019-06-24.
//

#ifndef MESHSIMPLIFICATION_REFINED_PAIRHEAP_H
#define MESHSIMPLIFICATION_REFINED_PAIRHEAP_H

#include "Const.h"

class Pair {
public:
    int time;
    int v0, v1;
    double cost;
    double dim[3]; // v_overline


    Pair(int _v0, int _v1) : v0(_v0), v1(_v1) {
        if (v0 > v1)
            std::swap(v0, v1);
    }

    inline bool operator<(const Pair &b) const {
        return cost > b.cost;

    }
    inline std::pair<int,int> getValue() const {
        return std::make_pair(v0,v1);
    }

    friend bool operator==(const Pair &left, const Pair &right) {
        return (left.v0 == right.v0 && left.v1 == right.v1);
    }

    friend std::ostream &operator<<(std::ostream &os, const Pair &pair) {
        os << "Pair: " << pair.v0 << " " << pair.v1;
        return os;
    }
};


class Heap {
public:
    std::map<std::pair<int,int >, int> mapper;
    std::priority_queue<Pair> queue;
    std::vector<std::set<int>> pairs;
    int stamp = 0;

    Heap() = default;

    inline bool count(std::pair<int, int> inquiry) {
        return mapper.count(inquiry) ? (mapper[inquiry] > 0) : false;
    }
    inline void refresh(){
        while(queue.size() > 0){
            const Pair& pair = queue.top();
            int time = mapper[pair.getValue()];
            /*检验top的pair是否有效*/
            if(time == pair.time)
                break;
            queue.pop();
        }
    }
    inline Pair top(){
        refresh();
        return queue.top();
    }
    inline void pop() {
        refresh();
        queue.pop();
        return;
    }
    inline void insert(Pair& pair) {
        pair.time = ++ stamp;
        mapper[pair.getValue()] = stamp;
        queue.push(pair);
        pairs[pair.v0].insert(pair.v1);
        pairs[pair.v1].insert(pair.v0);
        return;
    }
    inline void remove(const Pair &pair) {
        mapper[pair.getValue()] = -1;
        pairs[pair.v0].erase(pair.v1);
        pairs[pair.v1].erase(pair.v0);
        return;
    }

};


#endif //MESHSIMPLIFICATION_REFINED_PAIRHEAP_H
