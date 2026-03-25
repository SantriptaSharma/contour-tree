#include "TopologicalFeatures.hpp"

#include <fstream>
#include <iostream>
#include <cassert>

#include "constants.h"
#include <iomanip>
#include <filesystem>

namespace contourtree {

TopologicalFeatures::TopologicalFeatures() {}

void TopologicalFeatures::loadData(std::string dataLocation) {
    ctdata = ContourTreeData();
    ctdata.loadBinFile(dataLocation);

    // read order file
    // read meta data
    std::ifstream ip(dataLocation + ".order.dat");
    int64_t orderSize;
    ip >> orderSize;
    ip.close();

    order.resize(orderSize);
    wts.resize(orderSize);

    std::string binFile = dataLocation + ".order.bin";
    std::ifstream bin(binFile, std::ios::binary);
    bin.read((char*)order.data(), order.size() * sizeof(uint32_t));
    bin.read((char*)wts.data(), wts.size() * sizeof(float));
    bin.close();

    gsim.setInput(&ctdata);
    gsim.simplify(order, 1, 0, wts);
}

void TopologicalFeatures::addFeature(SimplifyCT& sim, uint32_t bno, std::vector<Feature>& features,
                                     std::set<size_t>& featureSet) {
    Branch b1 = sim.branches.at(bno);
    Feature f;
    f.from = ctdata.nodeVerts[b1.from];
    f.to = ctdata.nodeVerts[b1.from];

    std::deque<size_t> queue;
    queue.push_back(bno);
    while (queue.size() > 0) {
        size_t b = queue.front();
        queue.pop_front();
        if (b != bno && featureSet.find(b) != featureSet.end()) {
            // this cannot happen
            assert(false);
        }
        featureSet.insert(b);
        Branch br = sim.branches.at(b);
        f.arcs.insert(f.arcs.end(), br.arcs.data(), br.arcs.data() + br.arcs.size());
        for (int i = 0; i < br.children.size(); i++) {
            uint32_t bc = br.children.at(i);
            queue.push_back(bc);
        }
    }
    features.push_back(f);
}

inline bool isPathPresent(ContourTreeData &ctdata, uint32_t from, uint32_t to) {
    std::deque<uint32_t> queue;
    queue.push_back(from);
    while(queue.size() > 0) {
        uint32_t v = queue.front();
        queue.pop_front();
        if(v == to) {
            return true;
        }
        if(ctdata.fnVals[v] <= ctdata.fnVals[to]) {
            for(auto ano: ctdata.nodes[v].next) {
                assert(ctdata.arcs[ano].from == v);
                queue.push_back(ctdata.arcs[ano].to);
            }
        }
    }
    return false;
}

std::vector<Feature> TopologicalFeatures::getPartitionedExtremaFeatures(int &topk, float th) {
    std::vector<Feature> features;

    std::set<size_t> featureSet;
    if (topk == -1) {
        topk = 0;
        for (int i = order.size() - 1; i >= 0; i--) {
            if (wts[i] > th) {
                topk++;
                featureSet.insert(order[i]);
            } else {
                break;
            }
        }
    }
    if (topk == 0) topk = 1;

    for (int _i = 0; _i < topk; _i++) {
        size_t i = order.size() - _i - 1;
        Branch b1 = gsim.branches.at(order[i]);
        Feature f;
        f.from = ctdata.nodeVerts[b1.from];
        f.to = ctdata.nodeVerts[b1.to];

        size_t bno = order[i];
        std::deque<size_t> queue;
        queue.push_back(bno);
        while (queue.size() > 0) {
            size_t b = queue.front();
            queue.pop_front();
            if (b != bno && featureSet.find(b) != featureSet.end()) {
                continue;
            }
            Branch br = gsim.branches.at(b);
            f.arcs.insert(f.arcs.end(), br.arcs.data(), br.arcs.data() + br.arcs.size());
            for (int i = 0; i < br.children.size(); i++) {
                int bc = br.children.at(i);
                queue.push_back(bc);
            }
        }
        features.push_back(f);
    }
    return features;
}

inline void addArcs(size_t bno, Feature &f, const SimplifyCT &sim) {
    std::deque<size_t> queue;
    queue.push_back(bno);
    while (queue.size() > 0) {
        size_t b = queue.front();
        queue.pop_front();
        Branch br = sim.branches.at(b);
        f.arcs.insert(f.arcs.end(), br.arcs.data(), br.arcs.data() + br.arcs.size());
        for (int i = 0; i < br.children.size(); i++) {
            int bc = br.children.at(i);
            queue.push_back(bc);
        }
    }
}

std::vector<Feature> TopologicalFeatures::getArcFeatures(int &topk, float th) {
    SimplifyCT sim;
    sim.setInput(&ctdata);

    topk = sim.simplify(order, topk, th, wts);

    std::vector<Feature> features;
    std::set<size_t> featureSet;
    for (size_t _i = 0; _i < sim.branches.size(); _i++) {
        if (sim.removed[_i]) {
            continue;
        }
        featureSet.insert(_i);
    }
    for (size_t _i = 0; _i < sim.branches.size(); _i++) {
        size_t i = _i;
        if (sim.removed[i]) {
            continue;
        }
        Branch b1 = sim.branches.at(i);
        Feature f;
        f.from = ctdata.nodeVerts[b1.from];
        f.to = ctdata.nodeVerts[b1.to];

        addArcs(i,f,sim);
        // add arcs from multi saddles
        if(sim.arcArrayUpper[b1.from].size() > 0) {
            // find if b1.from to uv shares the path from b1.from to b1.to
            // for this it is sufficient to check whether there exists a non-decreasing path from b1.to to uv
            uint32_t uv = gsim.vArrayNext[b1.from];
            if(isPathPresent(this->ctdata,b1.to,uv)) {
                for(auto a : sim.arcArrayUpper[b1.from]) {
                    addArcs(a,f,sim);
                }
            }
        }
        if(sim.arcArrayLower[b1.to].size() > 0) {
            // find if b1.to to lv shares the path from b1.to to b1.from
            // for this it is sufficient to check whether there exists a non-increasing path from b1.from to uv
            uint32_t lv = gsim.vArrayPrev[b1.to];
            if(isPathPresent(this->ctdata,lv,b1.from)) {
                for(auto a : sim.arcArrayLower[b1.to]) {
                    addArcs(a,f,sim);
                }
            }
        }
        features.push_back(f);
    }
    return features;
}

void generatePersistenceDiagram(std::string filePrefix) {
    TopologicalFeatures topoFeatures;
    topoFeatures.loadData(filePrefix);
    int topk = -1;
    std::string opfile = filePrefix+"-persistence.csv";
    std::ofstream op(opfile);
    auto features = topoFeatures.getPartitionedExtremaFeatures(topk, 0);
    for(auto f: features) {
        int32_t from = topoFeatures.ctdata.nodeMap[f.from];
        int32_t to = topoFeatures.ctdata.nodeMap[f.to];
        double birth = topoFeatures.ctdata.fnVals[from];
        double death = topoFeatures.ctdata.fnVals[to];
        op << std::setprecision (15) << birth << "," << death << "\n";
    }
    op.close();
}

void generatePersistenceDiagrams(std::string folder) {
    namespace fs = std::filesystem;
    for (const auto& entry : fs::recursive_directory_iterator(folder)) {
        if (entry.is_regular_file() && entry.path().extension() == ".bin") {
            std::string filename = entry.path().stem().string();
            if (filename.length() >= 3 && filename.substr(filename.length() - 3) == ".rg") {
                std::string prefix = entry.path().parent_path().string() + "/" + filename.substr(0, filename.length() - 3);
                std::cout << prefix << "\n";
                generatePersistenceDiagram(prefix);
            }
        }
    }
}

}  // namespace contourtree
