#include <iostream>
#include "Mesh.h"

int main(int argc, char** argv) {


    std::string output;
    std::string input;
    double ratio = 0.1;

    if(argc < 4) {
        std:: cout << "Usage: simplify <input> <output> <ratio>" << std:: endl;
        std:: cout << "Using default! " << std::endl;
        input = "Meshes/buddha.obj";
        output = "Output/dragon_0.1.obj";
    }
    else
        input = (argv[1]), output = argv[2], ratio = atof(argv[3]);


    time_t start = clock();
    time_t end = clock();

    Mesh mesh;
    if(!mesh.loadObj(input))
        return 0;

    mesh.calculateQ();
    mesh.selectPairs();

    const int TARGET = mesh.fcnt * ratio;
    for (int i = 0; mesh.fcnt > TARGET ; ++i) {
        mesh.aggregation();
    }

    mesh.writeObj(output);


    end = clock();
    std::cout << "calculating time cost: " << double(end - start)/CLOCKS_PER_SEC << " s"<< std::endl;


    return 0;

}