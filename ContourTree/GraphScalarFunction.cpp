#include "GraphScalarFunction.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <sstream>

namespace contourtree {

GraphScalarFunction::GraphScalarFunction() {
}

int GraphScalarFunction::getMaxDegree() {
    return maxStar;
}

int GraphScalarFunction::getVertexCount() {
    return nv;
}

int GraphScalarFunction::getStar(int64_t v, std::vector<int64_t> &star) {
    int ct = 0;
    for(uint32_t vv: vertices[v].adj) {
        star[ct] = vv;
        ct ++;
    }
    return ct;
}

bool GraphScalarFunction::lessThan(int64_t v1, int64_t v2) {
    if(fnVals[v1] < fnVals[v2]) {
        return true;
    } else if(fnVals[v1] == fnVals[v2]) {
        return (v1 < v2);
    }
    return false;
}

scalar_t GraphScalarFunction::getFunctionValue(int64_t v) {
    return this->fnVals[v];
}

void GraphScalarFunction::initialize(uint32_t noNodes) {
    this->nv = noNodes;

    this->vertices.clear();
    this->vertices.resize(nv);

    this->fnVals.clear();
    this->fnVals.resize(nv);

}

void GraphScalarFunction::loadGraphFromAdjList(const std::vector<std::vector<int64_t>> adjList) {
    
    initialize(adjList.size());

    for (size_t v1 = 0; v1 < adjList.size(); ++v1) {
        for (size_t j = 0; j < adjList[v1].size(); ++j) {
            int64_t v2 = adjList[v1][j];
            if (v1 != v2) {
                vertices[v1].adj.insert(static_cast<int64_t>(v2));
                vertices[v2].adj.insert(static_cast<int64_t>(v1));
            }
        }
    }

    maxStar = 0;
    for (size_t v = 0; v < vertices.size(); ++v) {
        int starSize = vertices[v].adj.size();
        if (starSize > maxStar) {
            maxStar = starSize;
        }
    }
}


inline std::vector<std::string> splitString(std::string s, char delim) {
    std::vector<std::string> ret;
    std::string t;
    std::stringstream ss(s);
    while (std::getline(ss, t, delim)) {
        ret.push_back(t);
    }
    return ret;
}

void GraphScalarFunction::loadGraph(std::string edgeFile, char sep) {
    std::ifstream ip(edgeFile);
    std::string s;

    std::vector<std::vector<int64_t>> adjList;

    while(std::getline(ip,s)) {
        
        if (s[0] == '#') continue; // skip comment lines
        if (s.empty()) continue;   // skip empty lines

        adjList.push_back(std::vector<int64_t>());
        auto &adj = adjList.back();

        std::vector<std::string> parts = splitString(s, sep);
        for(std::string &e: parts) {
            int o = std::stoi(e);

            // ignore self-connected edges
            if (o != (adjList.size() - 1)) {
                adj.push_back(o);
            }
        }
    }

    loadGraphFromAdjList(adjList);
}

void GraphScalarFunction::updateFnValues(const std::vector<scalar_t> &fn) {
    assert(fn.size() == nv);

    for(int i = 0;i < nv;i ++) {
        this->fnVals[i] = fn[i];
    }
}

}
