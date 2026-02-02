#include "SimplifyCT.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>

namespace contourtree {

bool BranchCompare::operator()(uint32_t v1, uint32_t v2) { return sim->compare(v1, v2); }

SimplifyCT::SimplifyCT() {
    queue = std::priority_queue<uint32_t, std::vector<uint32_t>, BranchCompare>(BranchCompare(this));
    order.clear();
}

void SimplifyCT::setInput(ContourTreeData* data) { this->data = data; }

void SimplifyCT::addToQueue(uint32_t ano) {
    if (isCandidate(branches[ano])) {
        queue.push(ano);
        inq[ano] = true;
    }
}

bool SimplifyCT::isCandidate(const Branch& br) {
    uint32_t from = br.from;
    uint32_t to = br.to;
    if (nodes[from].prev.size() == 0) {
        // minimum
        if (nodes[to].prev.size() > 1) {
            return true;
        } else {
            return false;
        }
    }
    if (nodes[to].next.size() == 0) {
        // maximum
        if (nodes[from].next.size() > 1) {
            return true;
        } else {
            return false;
        }
    }
    return false;
}

void SimplifyCT::initSimplification(SimFunction* f) {
    branches.resize(data->noArcs);
    nodes.resize(data->noNodes);
    for (uint32_t i = 0; i < branches.size(); i++) {
        branches[i].from = data->arcs[i].from;
        branches[i].to = data->arcs[i].to;
        branches[i].parent = -1;
        branches[i].arcs.push_back(i);

        nodes[branches[i].from].next.push_back(i);
        nodes[branches[i].to].prev.push_back(i);
    }

    fn.resize(branches.size());
    removed.resize(branches.size(), false);
    invalid.resize(branches.size(), false);
    inq.resize(branches.size(), false);

    arcArrayUpper.resize(nodes.size());
    arcArrayLower.resize(nodes.size());

    vArrayPrev.resize(nodes.size(),-1);
    vArrayNext.resize(nodes.size(),-1);

    simFn = f;
    if (f != NULL) {
        simFn->init(fn, branches);

        for (uint32_t i = 0; i < branches.size(); i++) {
            addToQueue(i);
        }
    }
}

bool SimplifyCT::compare(uint32_t b1, uint32_t b2) const {
    // If I want smallest weight on top, I need to return true if b1 > b2 (sort in descending order)
    if (fn[b1] > fn[b2]) {
        return true;
    }
    if (fn[b1] < fn[b2]) {
        return false;
    }
    float p1 = data->fnVals[branches[b1].to] - data->fnVals[branches[b1].from];
    float p2 = data->fnVals[branches[b2].to] - data->fnVals[branches[b2].from];
    if (p1 > p2) {
        return true;
    }
    if (p1 < p2) {
        return false;
    }
    int diff1 = branches[b1].to - branches[b1].from;
    int diff2 = branches[b2].to - branches[b2].from;
    if (diff1 > diff2) {
        return true;
    }
    if (diff1 < diff2) {
        return false;
    }
    return (branches[b2].from > branches[b1].from);
}

void SimplifyCT::removeArc(uint32_t ano) {
    Branch br = branches[ano];
    uint32_t from = br.from;
    uint32_t to = br.to;
    uint32_t mergedVertex = -1;
    if (nodes[from].prev.size() == 0) {
        // minimum
        mergedVertex = to;
        arcArrayLower[mergedVertex].push_back(ano);
    }
    if (nodes[to].next.size() == 0) {
        // maximum
        mergedVertex = from;
        arcArrayUpper[mergedVertex].push_back(ano);
    }
    nodes[from].next.erase(std::remove(nodes[from].next.begin(), nodes[from].next.end(), ano),
                           nodes[from].next.end());
    nodes[to].prev.erase(std::remove(nodes[to].prev.begin(), nodes[to].prev.end(), ano),
                         nodes[to].prev.end());
    removed[ano] = true;

    if (nodes[mergedVertex].prev.size() == 1 && nodes[mergedVertex].next.size() == 1) {
        mergeVertex(mergedVertex);
    }
    if (simFn != NULL) simFn->branchRemoved(branches, ano, invalid);
}

void SimplifyCT::mergeVertex(uint32_t v) {
    uint32_t prev = nodes[v].prev.at(0);
    uint32_t next = nodes[v].next.at(0);

    vArrayPrev[v] = branches[prev].from;
    vArrayNext[v] = branches[next].to;

    int a = -1;
    int rem = -1;
    if (inq[prev]) {
        invalid[prev] = true;
        removed[next] = true;
        branches[prev].to = branches[next].to;
        a = prev;
        rem = next;

        for (int i = 0; i < nodes[branches[prev].to].prev.size(); i++) {
            if (nodes[branches[prev].to].prev[i] == next) {
                nodes[branches[prev].to].prev[i] = prev;
            }
        }
    } else {
        invalid[next] = true;
        removed[prev] = true;
        branches[next].from = branches[prev].from;
        a = next;
        rem = prev;

        for (int i = 0; i < nodes[branches[next].from].next.size(); i++) {
            if (nodes[branches[next].from].next[i] == prev) {
                nodes[branches[next].from].next[i] = next;
            }
        }
        if (simFn != NULL && !inq[next]) {
            addToQueue(next);
        }
    }
    for (int i = 0; i < branches[rem].children.size(); i++) {
        int ch = branches[rem].children.at(i);
        branches[a].children.push_back(ch);
        assert(branches[ch].parent == rem);
        branches[ch].parent = a;
    }
    //    branches[a].arcs << branches[rem].arcs;
    branches[a].arcs.insert(branches[a].arcs.end(), branches[rem].arcs.begin(),
                            branches[rem].arcs.end());
    for (int i = 0; i < arcArrayLower[v].size(); i++) {
        uint32_t aa = arcArrayLower[v].at(i);
        branches[a].children.push_back(aa);
        branches[aa].parent = a;
    }
    for (int i = 0; i < arcArrayUpper[v].size(); i++) {
        uint32_t aa = arcArrayUpper[v].at(i);
        branches[a].children.push_back(aa);
        branches[aa].parent = a;
    }
    branches[rem].parent = -2;
}

void SimplifyCT::simplify(SimFunction* simFn) {
    std::cout << "init" << std::endl;
    initSimplification(simFn);

    std::cout << "going over priority queue" << std::endl;
    while (queue.size() > 0) {
        uint32_t ano = queue.top();
        queue.pop();
        inq[ano] = false;
        if (!removed[ano]) {
            if (invalid[ano]) {
                simFn->update(branches, ano);
                invalid[ano] = false;
                addToQueue(ano);
            } else {
                if (isCandidate(branches[ano])) {
                    removeArc(ano);
                    order.push_back(ano);
                }
            }
        }
    }
    std::cout << "pass over removed" << std::endl;
    int root = 0;
    for (int i = 0; i < removed.size(); i++) {
        if (!removed[i]) {
            //            assert(root == 0);
            order.push_back(i);
            root++;
        }
    }
}

int SimplifyCT::simplify(const std::vector<uint32_t>& order, int topk, float th,
                          const std::vector<float>& wts) {
    std::cout << "init" << std::endl;
    initSimplification(NULL);

    std::cout << "going over order queue" << std::endl;
    for (int i = 0; i < order.size(); i++) {
        inq[order.at(i)] = true;
    }
    int removed = 0;
    if (topk > 0) {
        size_t ct = order.size() - topk;
        for (int i = 0; i < ct; i++) {
            uint32_t ano = order.at(i);
            if (!isCandidate(branches[ano])) {
                std::cout << "failing candidate test" << std::endl;
                assert(false);
            }
            inq[ano] = false;
            removeArc(ano);
        }
    } else if (th != 0) {
        for (int i = 0; i < order.size() - 1; i++) {
            uint32_t ano = order.at(i);
            if (!isCandidate(branches[ano])) {
                std::cout << "failing candidate test" << std::endl;
                assert(false);
            }
            float fn = wts.at(i);
            if (fn > th) {
                break;
            }
            inq[ano] = false;
            removeArc(ano);
            removed ++;
        }
        topk = order.size() - removed;
    }
    return topk;
}

void SimplifyCT::outputOrder(std::string fileName, bool normalize) {
    std::cout << "Writing meta data" << std::endl;
    {
        std::ofstream pr(fileName + ".order.dat");
        pr << order.size() << "\n";
        pr.close();
    }
    std::vector<float> wts;
    for (size_t i = 0; i < order.size(); i++) {
        uint32_t ano = order.at(i);
        float val = this->simFn->getBranchWeight(ano);
        wts.push_back(val);
    }

    // normalize weights
    if(normalize) {
        float maxWt = wts.at(wts.size() - 1);
        if (maxWt == 0) maxWt = 1;
        for (int i = 0; i < wts.size(); i++) {
            wts[i] /= maxWt;
        }
    }

    std::cout << "writing tree output" << std::endl;
    std::string binFile = fileName + ".order.bin";
    std::ofstream of(binFile, std::ios::binary);
    of.write((char*)order.data(), order.size() * sizeof(uint32_t));
    of.write((char*)wts.data(), wts.size() * sizeof(float));
    of.close();
}

void SimplifyCT::getSimplificationPlot(const std::vector<uint32_t> &order, const std::vector<float> &wts, std::vector<float> &fns, 
                                        std::vector<int32_t> &minct, std::vector<int32_t> &maxct) {
    initSimplification(NULL);

    std::cout << "going over order queue" << std::endl;
    for (int i = 0; i < order.size(); i++) {
        inq[order.at(i)] = true;
    }

    fns.clear();
    minct.clear();
    maxct.clear();
    fns.push_back(0);
    minct.push_back(0);
    maxct.push_back(0);

    int minrem = 0;
    int maxrem = 0;
    for (int i = 0; i < order.size() - 1; i++) {
        uint32_t ano = order.at(i);
        if (!isCandidate(branches[ano])) {
            std::cout << "failing candidate test" << std::endl;
            assert(false);
        }
        float fn = wts.at(i);
        uint32_t from = branches[ano].from;
        if(data->type[from] == MINIMUM) {
            minrem ++;
        }
        uint32_t to = branches[ano].to;
        if(data->type[to] == MAXIMUM) {
            maxrem ++;
        }
        inq[ano] = false;
        removeArc(ano);

        fns.push_back(fn);
        minct.push_back(minrem);
        maxct.push_back(maxrem);
    }

    for(int i = 0;i < fns.size();i ++) {
        minct[i] = minrem - minct[i] + 1;
        maxct[i] = maxrem - maxct[i] + 1;
        if(i > 0) {
            fns[i] = std::max(fns[i],fns[i-1]);
        }
    }
}

// TODO: fix the saddle behaviour on this, seems incorrect 
void SimplifyCT::getFilteredSimplificationPlot(const std::vector<uint32_t>& order, const std::vector<float>& wts, std::vector<float>& fns, 
                                               float minfnstart, float maxfnstart, float minfnend, float maxfnend, const std::set<char> &types,
                                               std::vector<int32_t>& filteredct) {
        initSimplification(NULL);

        for (int i = 0; i < order.size(); i++) {
            inq[order.at(i)] = true;
        }

        fns.clear();
        filteredct.clear();
        fns.push_back(0);
        filteredct.push_back(0);

        int removed = 0;
        for (int i = 0; i < order.size() - 1; i++) {
            uint32_t ano = order.at(i);
            if (!isCandidate(branches[ano])) {
                std::cout << "failing candidate test" << std::endl;
                assert(false);
            }
            float fn = wts.at(i);
            uint32_t from = branches[ano].from;
            uint32_t to = branches[ano].to;

            char type_from = data->type[from];
            char type_to = data->type[to];

            bool is_valid_saddle = (type_from == SADDLE && type_to == SADDLE) && types.find(SADDLE) != types.end();
            bool is_valid_minimum = (type_from == MINIMUM) && types.find(MINIMUM) != types.end();
            bool is_valid_maximum = (type_to == MAXIMUM) && types.find(MAXIMUM) != types.end();
            bool is_within_bounds = (data->fnVals[from] >= minfnstart && data->fnVals[from] <= maxfnstart) && (data->fnVals[to] >= minfnend && data->fnVals[to] <= maxfnend);

            if((is_valid_minimum || is_valid_saddle || is_valid_maximum) && is_within_bounds) removed++;    

            inq[ano] = false;
            removeArc(ano);

            fns.push_back(fn);
            filteredct.push_back(removed);
        }

        for(int i = 0;i < fns.size();i ++) {
            filteredct[i] = removed - filteredct[i] + 1;
            if(i > 0) {
                fns[i] = std::max(fns[i],fns[i-1]);
            }
        }
    }

    static uint32_t getNumClasses(std::vector<uint32_t> const &labels) {
        uint32_t num_classes = 0;
        for (uint32_t label : labels) {
            if (label + 1 > num_classes) {
                num_classes = label + 1;
            }
        }
        return num_classes;
    }
    
    static void computeArcClassCounts(const std::vector<uint32_t> &partition, const std::vector<uint32_t> &labels, 
                                      uint32_t num_classes, uint32_t num_arcs, std::vector<std::vector<uint32_t>> &arcClassCounts) {
        arcClassCounts.resize(num_arcs, std::vector<uint32_t>(num_classes, 0));
        
        for (uint32_t i = 0; i < partition.size(); i++) {
            uint32_t arc_no = partition[i];
            uint32_t label = labels[i];
            if (arc_no < arcClassCounts.size() && label < num_classes) {
                arcClassCounts[arc_no][label]++;
            }
        }
    }

    static std::pair<uint32_t, float> getMajorityClass(uint32_t ano, std::vector<Branch> const & branches, 
                                                       const std::vector<std::vector<uint32_t>> &arcClassCounts) {
        if (ano >= arcClassCounts.size()) return {0, 0.0f};
        
        std::vector<uint32_t> branchClassCounts(arcClassCounts[0].size(), 0);
        
        Branch b = branches[ano];
        for (uint32_t arc_id : b.arcs) {
            for (uint32_t c = 0; c < arcClassCounts[arc_id].size(); c++) {
                branchClassCounts[c] += arcClassCounts[arc_id][c];
            }
        }
        
        uint32_t total_points = 0;
        uint32_t max_class_count = 0;
        uint32_t majority_class = 0;
        for (uint32_t c = 0; c < branchClassCounts.size(); c++) {
            uint32_t count = branchClassCounts[c];
            total_points += count;
            if (count > max_class_count) {
                max_class_count = count;
                majority_class = c;
            }
        }
        if (total_points == 0) return {0, 0.0f};
        float proportion = static_cast<float>(max_class_count) / static_cast<float>(total_points);
        return {majority_class, proportion};
    }

    void SimplifyCT::getHomoValleyPlot(const std::vector<uint32_t>& order, const std::vector<float>& wts, std::vector<float>& fns, 
                                             std::vector<int32_t>& remainingct, const std::vector<uint32_t> &labels, 
                                             float homogeneity_threshold, const std::vector<uint32_t> &partition) {
        initSimplification(NULL);

        for (int i = 0; i < order.size(); i++) {
            inq[order.at(i)] = true;
        }

        fns.clear();
        remainingct.clear();

        std::vector<std::vector<uint32_t>> arcClassCounts;
        computeArcClassCounts(partition, labels, getNumClasses(labels), branches.size(), arcClassCounts);

        auto countHomoValleys = [&]() {
            int count = 0;
            for (uint32_t b = 0; b < branches.size(); b++) {
                auto from = branches[b].from;
                if (inq[b] && data->type[from] == MINIMUM) {
                    auto [majority_class, proportion] = getMajorityClass(b, branches, arcClassCounts);
                    if (proportion >= homogeneity_threshold) {
                        count++;
                    }
                }
            }
            return count;
        };

        fns.push_back(0);
        remainingct.push_back(countHomoValleys());

        for (int i = 0; i < order.size() - 1; i++) {
            uint32_t ano = order.at(i);
            if (!isCandidate(branches[ano])) {
                std::cout << "failing candidate test" << std::endl;
                assert(false);
            }

            float fn = wts.at(i);
            
            inq[ano] = false;
            removeArc(ano);
            
            fns.push_back(fn);
            remainingct.push_back(countHomoValleys());
        }
    }

    // IGNORE
    // void SimplifyCT::getFilteredSimplificationPlotHomogeneity(const std::vector<uint32_t>& order, const std::vector<float>& wts, std::vector<float>& fns, 
    //                                            float minfnstart, float maxfnstart, float minfnend, float maxfnend, const std::set<char> &types,
    //                                            std::vector<int32_t>& filteredct, const std::vector<uint32_t> &labels, float homogeneity_threshold,
    //                                            const std::vector<uint32_t> &partition, std::vector<uint32_t> &branch_total_sizes,
    //                                            std::vector<uint32_t> &branch_majority_labels, std::vector<uint32_t> &branch_majority_sizes,
    //                                            std::vector<bool> &branch_was_homogeneous, std::vector<bool> &caused_homogeneous_destruction) {
    //     initSimplification(NULL);

    //     for (int i = 0; i < order.size(); i++) {
    //         inq[order.at(i)] = true;
    //     }

    //     uint32_t num_classes = getNumClasses(labels);
    //     std::vector<std::vector<uint32_t>> branchClassCounts = computeBranchClassCounts(partition, labels, num_classes, branches.size());

    //     fns.clear();
    //     filteredct.clear();
    //     branch_total_sizes.clear();
    //     branch_majority_labels.clear();
    //     branch_majority_sizes.clear();
    //     branch_was_homogeneous.clear();
    //     caused_homogeneous_destruction.clear();
    //     fns.push_back(0);
    //     filteredct.push_back(0);
    //     branch_total_sizes.push_back(0);
    //     branch_majority_labels.push_back(0);
    //     branch_majority_sizes.push_back(0);
    //     branch_was_homogeneous.push_back(false);
    //     caused_homogeneous_destruction.push_back(false);

    //     int removed = 0;
    //     for (int i = 0; i < order.size() - 1; i++) {
    //         uint32_t ano = order.at(i);
    //         if (!isCandidate(branches[ano])) {
    //             std::cout << "failing candidate test" << std::endl;
    //             assert(false);
    //         }
    //         float fn = wts.at(i);
    //         uint32_t from = branches[ano].from;
    //         uint32_t to = branches[ano].to;

    //         char type_from = data->type[from];
    //         char type_to = data->type[to];

    //         bool is_valid_saddle = (type_from == SADDLE && type_to == SADDLE) && types.find(SADDLE) != types.end();
    //         bool is_valid_minimum = (type_from == MINIMUM) && types.find(MINIMUM) != types.end();
    //         bool is_valid_maximum = (type_to == MAXIMUM) && types.find(MAXIMUM) != types.end();
    //         bool is_within_bounds = (data->fnVals[from] >= minfnstart && data->fnVals[from] <= maxfnstart) && (data->fnVals[to] >= minfnend && data->fnVals[to] <= maxfnend);

    //         // Record metadata for the branch being removed
    //         uint32_t total_size = 0;
    //         uint32_t majority_size = 0;
    //         uint32_t majority_label = 0;
    //         for (uint32_t c = 0; c < num_classes; c++) {
    //             uint32_t count = branchClassCounts[ano][c];
    //             total_size += count;
    //             if (count > majority_size) {
    //                 majority_size = count;
    //                 majority_label = c;
    //             }
    //         }
    //         float proportion = (total_size > 0) ? static_cast<float>(majority_size) / static_cast<float>(total_size) : 0.0f;
    //         bool was_homogeneous = (proportion > homogeneity_threshold);

    //         // Before removing the arc, track which branches need to receive this branch's class counts
    //         // We need to examine what removeArc and mergeVertex do
    //         uint32_t mergedVertex = static_cast<uint32_t>(-1);
    //         if (nodes[from].prev.size() == 0) {
    //             mergedVertex = to;
    //         }
    //         if (nodes[to].next.size() == 0) {
    //             mergedVertex = from;
    //         }

    //         // After removeArc is called, if mergeVertex happens, we need to transfer class counts
    //         // We'll track the potential merge before calling removeArc
    //         std::vector<uint32_t> transfer_from;
    //         std::vector<uint32_t> transfer_to;
            
    //         if (mergedVertex != static_cast<uint32_t>(-1) && 
    //             nodes[mergedVertex].prev.size() == 1 && nodes[mergedVertex].next.size() == 1) {
    //             uint32_t prev = nodes[mergedVertex].prev.at(0);
    //             uint32_t next = nodes[mergedVertex].next.at(0);
                
    //             // One of prev/next will remain active, the other will be removed and merged into it
    //             if (inq[prev]) {
    //                 // prev stays active (a), next is removed (rem)
    //                 // Children of next will get prev as parent
    //                 for (uint32_t ch : branches[next].children) {
    //                     transfer_from.push_back(ch);
    //                     transfer_to.push_back(prev);
    //                 }
    //                 // next itself merges into prev
    //                 transfer_from.push_back(next);
    //                 transfer_to.push_back(prev);
    //             } else {
    //                 // next stays active (a), prev is removed (rem)
    //                 // Children of prev will get next as parent
    //                 for (uint32_t ch : branches[prev].children) {
    //                     transfer_from.push_back(ch);
    //                     transfer_to.push_back(next);
    //                 }
    //                 // prev itself merges into next
    //                 transfer_from.push_back(prev);
    //                 transfer_to.push_back(next);
    //             }
    //         }

    //         inq[ano] = false;
    //         removeArc(ano);

    //         // Track whether this removal caused any homogeneous destruction
    //         bool this_removal_destroyed_homogeneity = false;
            
    //         // Now perform the class count transfers and check if homogeneity is destroyed
    //         for (size_t t = 0; t < transfer_from.size(); t++) {
    //             uint32_t from_branch = transfer_from[t];
    //             uint32_t to_branch = transfer_to[t];
    //             if (from_branch < arcClassCounts.size() && to_branch < arcClassCounts.size()) {
    //                 // Get the majority class and proportion before merging
    //                 auto [old_majority_class, old_proportion] = getMajorityClass(to_branch, branches, branchClassCounts);
    //                 bool was_homogeneous = (old_proportion > homogeneity_threshold);
                    
    //                 // Transfer class counts
    //                 for (uint32_t c = 0; c < num_classes; c++) {
    //                     branchClassCounts[to_branch][c] += branchClassCounts[from_branch][c];
    //                 }
                    
    //                 // Get the majority class and proportion after merging
    //                 auto [new_majority_class, new_proportion] = getMajorityClass(to_branch);
    //                 bool is_still_homogeneous_for_same_class = (new_proportion > homogeneity_threshold) && (new_majority_class == old_majority_class);
                    
    //                 // If it was homogeneous for a class but no longer homogeneous for that same class, count it
    //                 if (was_homogeneous && !is_still_homogeneous_for_same_class && 
    //                     (is_valid_minimum || is_valid_saddle || is_valid_maximum) && is_within_bounds) {
    //                     removed++;
    //                     this_removal_destroyed_homogeneity = true;
    //                 }
    //             }
    //         }

    //         fns.push_back(fn);
    //         filteredct.push_back(removed);
    //         branch_total_sizes.push_back(total_size);
    //         branch_majority_labels.push_back(majority_label);
    //         branch_majority_sizes.push_back(majority_size);
    //         branch_was_homogeneous.push_back(was_homogeneous);
    //         caused_homogeneous_destruction.push_back(this_removal_destroyed_homogeneity);
    //     }

    //     for(int i = 0;i < fns.size();i ++) {
    //         filteredct[i] = removed - filteredct[i] + 1;
    //         if(i > 0) {
    //             fns[i] = std::max(fns[i],fns[i-1]);
    //         }
    //     }
    // }    

}  // namespace contourtree
