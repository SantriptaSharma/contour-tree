#ifndef SIMPLIFYCT_HPP
#define SIMPLIFYCT_HPP

#include "ContourTreeData.hpp"
#include "SimFunction.hpp"
#include <queue>
#include <vector>
#include <set>

#if defined(WIN32)
#include <functional>
#endif

namespace contourtree {

class SimplifyCT;

struct BranchCompare {
    BranchCompare() {}
    BranchCompare(const SimplifyCT* simct) : sim(simct) {}
    bool operator()(uint32_t v1, uint32_t v2);

    const SimplifyCT* sim;
};

class SimplifyCT {
public:
    SimplifyCT();

    void setInput(ContourTreeData* data);
    void simplify(SimFunction* simFn);
    int simplify(const std::vector<uint32_t>& order, int topk = -1, float th = 0,
                  const std::vector<float>& wts = std::vector<float>());
    void outputOrder(std::string fileName, bool normalize);
    void getSimplificationPlot(const std::vector<uint32_t> &order, const std::vector<float> &wts, std::vector<float> &fns,
                                std::vector<int32_t> &minct, std::vector<int32_t> &maxct);
    void getFilteredSimplificationPlot(const std::vector<uint32_t> &order, const std::vector<float> &wts, std::vector<float> &fns, 
                                        float minfnstart, float maxfnstart, float minfnend, float maxfnend, const std::set<char> &types,
                                        std::vector<int32_t> &filteredct);
    void getHomoValleyPlot(const std::vector<uint32_t>& order, const std::vector<float>& wts, std::vector<float>& fns, 
                                 std::vector<int32_t>& remainingct, const std::vector<uint32_t> &labels, 
                                 float homogeneity_threshold, const std::vector<uint32_t> &partition);
    void getHomoValleyPlotPlusCoverages(const std::vector<uint32_t>& order, const std::vector<float>& wts, std::vector<float>& fns, 
                                 std::vector<int32_t>& remainingct, std::vector<int32_t>& remaininghomoct, const std::vector<uint32_t> &labels, float homogeneity_threshold, const std::vector<uint32_t> &partition,
                                 std::vector<std::vector<double>>& maj_class_homo_coverages, std::vector<std::vector<int>>& maj_class_homo_counts,
                                 std::vector<std::vector<double>>& class_homo_coverages, std::vector<std::vector<double>>& class_coverages);
    // void getFilteredSimplificationPlotHomogeneity(const std::vector<uint32_t> &order, const std::vector<float> &wts, std::vector<float> &fns, 
    //                                     float minfnstart, float maxfnstart, float minfnend, float maxfnend, const std::set<char> &types,
    //                                     std::vector<int32_t> &filteredct, const std::vector<uint32_t> &labels, float homogeneity_threshold,
    //                                     const std::vector<uint32_t> &partition, std::vector<uint32_t> &branch_total_sizes,
    //                                     std::vector<uint32_t> &branch_majority_labels, std::vector<uint32_t> &branch_majority_sizes,
    //                                     std::vector<bool> &branch_was_homogeneous, std::vector<bool> &caused_homogeneous_destruction);

protected:
    void initSimplification(SimFunction* f);
    void addToQueue(uint32_t ano);
    bool isCandidate(const Branch& br);
    void removeArc(uint32_t ano);
    void mergeVertex(uint32_t v);

public:
    bool compare(uint32_t b1, uint32_t b2) const;

public:
    const ContourTreeData* data;
    std::vector<Branch> branches;
    std::vector<Node> nodes;

    std::vector<float> fn;
    std::vector<float> fnv;
    std::vector<bool> invalid;
    std::vector<bool> removed;
    std::vector<bool> inq;
    SimFunction* simFn;

    std::priority_queue<uint32_t, std::vector<uint32_t>, BranchCompare> queue;
    std::vector<uint32_t> order;
    std::vector<std::vector<uint32_t>> arcArrayUpper, arcArrayLower;
    std::vector<uint32_t> vArrayPrev, vArrayNext;
};

}  // namespace contourtree

#endif  // SIMPLIFYCT_HPP
