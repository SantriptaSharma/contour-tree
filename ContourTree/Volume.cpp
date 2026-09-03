#include "Volume.hpp"

#include <fstream>
#include <iostream>

namespace contourtree {

Volume::Volume(const ContourTreeData& ctData, std::string partFile) {
    std::ifstream bin(partFile, std::ios::binary | std::ios::ate);
    uint32_t size = bin.tellg();
    std::cout << "part size: " << size << std::endl;
    bin.close();

    std::vector<uint32_t> cols(size / sizeof(uint32_t));
    bin.open(partFile, std::ios::binary);
    bin.read((char*)cols.data(), size);
    bin.close();

    vol.resize(ctData.noArcs, 0);
    brVol.resize(ctData.noArcs, 0);
    initVolumes(cols);
}

void Volume::initVolumes(const std::vector<uint32_t>& cols) {
    for (size_t i = 0; i < cols.size(); i++) {
        if (cols[i] == -1) {
            continue;
        }
        vol[cols[i]]++;
    }
}

void Volume::init(std::vector<float>& fn, std::vector<Branch>& br) {
    this->fn = fn.data();
    for (int i = 0; i < fn.size(); i++) {
        this->update(br, i);
    }
}

void Volume::update(const std::vector<Branch>& br, uint32_t brNo) {
    brVol[brNo] = 0;
    for (int i = 0; i < br[brNo].arcs.size(); i++) {
        brVol[brNo] += vol[br[brNo].arcs.at(i)];
    }
    for (int i = 0; i < br[brNo].children.size(); i++) {
        int child = br[brNo].children.at(i);
        brVol[brNo] += volume(br, child);
    }
    fn[brNo] = brVol[brNo];
}

float Volume::volume(const std::vector<Branch>& br, int brNo) {
    float val = 0;
    for (int i = 0; i < br[brNo].arcs.size(); i++) {
        val += vol[br[brNo].arcs.at(i)];
    }
    for (int i = 0; i < br[brNo].children.size(); i++) {
        int child = br[brNo].children.at(i);
        val += volume(br, child);
    }
    return val;
}

void Volume::branchRemoved(std::vector<Branch>&, uint32_t, std::vector<bool>&) {}

float Volume::getBranchWeight(uint32_t brNo) { return fn[brNo]; }

}  // namespace contourtree
