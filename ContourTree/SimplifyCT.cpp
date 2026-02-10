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
                if (!removed[b] && data->type[from] == MINIMUM) {
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

    void SimplifyCT::getHomoValleyPlotPlusCoverages(const std::vector<uint32_t>& order, const std::vector<float>& wts, std::vector<float>& fns,
                                                    std::vector<int32_t>& remainingct, std::vector<int32_t>& remaininghomoct, const std::vector<uint32_t>& labels,
                                                    float homogeneity_threshold, const std::vector<uint32_t>& partition,
                                                    std::vector<std::vector<double>>& maj_class_homo_coverages, std::vector<std::vector<int>>& maj_class_homo_counts,
                                                    std::vector<std::vector<double>>& class_homo_coverages, std::vector<std::vector<double>>& class_coverages) {
        initSimplification(NULL);

        for (int i = 0; i < order.size(); i++) {
            inq[order.at(i)] = true;
        }

        fns.clear();
        remainingct.clear();
        maj_class_homo_coverages.clear();
        maj_class_homo_counts.clear();
        class_homo_coverages.clear();
        class_coverages.clear();

        SimplifyCT gsim;
        gsim.setInput(const_cast<ContourTreeData*>(data));
        gsim.simplify(order, 1, 0, wts);

        // Compute arc class counts
        std::vector<std::vector<uint32_t>> arcClassCounts;
        uint32_t num_classes = getNumClasses(labels);
        computeArcClassCounts(partition, labels, num_classes, branches.size(), arcClassCounts);

        // Count total points per class
        std::vector<uint32_t> class_totals(num_classes, 0);
        for (uint32_t label : labels) {
            if (label < num_classes) {
                class_totals[label]++;
            }
        }

        // Helper lambda to check if path exists in contour tree
        auto isPathPresent = [&](uint32_t from, uint32_t to) -> bool {
            std::deque<uint32_t> queue;
            queue.push_back(from);
            while(queue.size() > 0) {
                uint32_t v = queue.front();
                queue.pop_front();
                if(v == to) {
                    return true;
                }
                if(data->fnVals[v] <= data->fnVals[to]) {
                    for(auto ano: data->nodes[v].next) {
                        assert(data->arcs[ano].from == v);
                        queue.push_back(data->arcs[ano].to);
                    }
                }
            }
            return false;
        };

        // Helper lambda to recursively add all arcs from a branch and its children
        auto addArcsToSet = [&](size_t bno, std::set<uint32_t>& arc_set) {
            std::deque<size_t> queue;
            queue.push_back(bno);
            while (queue.size() > 0) {
                size_t b = queue.front();
                queue.pop_front();
                Branch br = branches[b];
                arc_set.insert(br.arcs.begin(), br.arcs.end());
                for (int i = 0; i < br.children.size(); i++) {
                    int bc = br.children[i];
                    queue.push_back(bc);
                }
            }
        };

        // Lambda to compute all statistics
        auto computeStats = [&](std::vector<uint32_t>& maj_class_points, std::vector<uint32_t>& valleys_per_class,
                               std::vector<uint32_t>& class_homo_points, std::vector<uint32_t>& class_all_points) {
            int count = 0;
            int homo_count = 0;

            maj_class_points.assign(num_classes, 0);
            valleys_per_class.assign(num_classes, 0);
            class_homo_points.assign(num_classes, 0);
            class_all_points.assign(num_classes, 0);
            
            for (uint32_t b = 0; b < branches.size(); b++) {
                if (removed[b]) continue;
                
                auto from = branches[b].from;
                bool is_minimum = (data->type[from] == MINIMUM);
                
                if (!is_minimum) continue;
                count++;
                
                // Collect all arcs for this branch
                std::set<uint32_t> all_arcs;
                addArcsToSet(b, all_arcs);
                
                // Add arcs from multi-saddles using arcArrayUpper/Lower
                Branch b1 = branches[b];
                if(arcArrayUpper[b1.from].size() > 0) {
                    uint32_t uv = gsim.vArrayNext[b1.from];
                    if(isPathPresent(b1.to, uv)) {
                        for(auto a : arcArrayUpper[b1.from]) {
                            addArcsToSet(a, all_arcs);
                        }
                    }
                }
                if(arcArrayLower[b1.to].size() > 0) {
                    uint32_t lv = gsim.vArrayPrev[b1.to];
                    if(isPathPresent(lv, b1.from)) {
                        for(auto a : arcArrayLower[b1.to]) {
                            addArcsToSet(a, all_arcs);
                        }
                    }
                }
                
                auto [majority_class, proportion] = getMajorityClass(b, branches, arcClassCounts);
                bool is_homogeneous = (proportion >= homogeneity_threshold);
                
                // Count points from all classes in all valleys
                for (uint32_t arc_id : all_arcs) {
                    if (arc_id < arcClassCounts.size()) {
                        for (uint32_t c = 0; c < num_classes; c++) {
                            class_all_points[c] += arcClassCounts[arc_id][c];
                        }
                    }

                    if (!is_homogeneous) continue;
                    
                    if (arc_id < arcClassCounts.size() && majority_class < arcClassCounts[arc_id].size()) {
                        maj_class_points[majority_class] += arcClassCounts[arc_id][majority_class];
                    }

                    if (arc_id < arcClassCounts.size()) {
                        for (uint32_t c = 0; c < num_classes; c++) {
                            class_homo_points[c] += arcClassCounts[arc_id][c];
                        }
                    }
                }
                
                if (is_homogeneous) {
                    homo_count++;
                    valleys_per_class[majority_class]++;
                }
            }
            return std::make_pair(count, homo_count);
        };

        std::vector<uint32_t> maj_class_points;
        std::vector<uint32_t> valleys_per_class;
        std::vector<uint32_t> class_homo_points;
        std::vector<uint32_t> class_all_points;
        
        // Initial state
        fns.push_back(0);
        auto [ct, homo_ct] = computeStats(maj_class_points, valleys_per_class, class_homo_points, class_all_points);
        remainingct.push_back(ct);
        remaininghomoct.push_back(homo_ct);
        
        // Compute coverages for initial state
        std::vector<double> current_maj_class_homo_cov(num_classes);
        std::vector<int> current_maj_class_homo_cnt(num_classes);
        std::vector<double> current_class_homo_cov(num_classes);
        std::vector<double> current_class_cov(num_classes);
        
        for (uint32_t c = 0; c < num_classes; c++) {
            current_maj_class_homo_cnt[c] = valleys_per_class[c];
            if (class_totals[c] > 0) {
                current_maj_class_homo_cov[c] = static_cast<double>(maj_class_points[c]) / static_cast<double>(class_totals[c]);
                current_class_homo_cov[c] = static_cast<double>(class_homo_points[c]) / static_cast<double>(class_totals[c]);
                current_class_cov[c] = static_cast<double>(class_all_points[c]) / static_cast<double>(class_totals[c]);
            } else {
                current_maj_class_homo_cov[c] = 0.0;
                current_class_homo_cov[c] = 0.0;
                current_class_cov[c] = 0.0;
            }
        }
        maj_class_homo_coverages.push_back(current_maj_class_homo_cov);
        maj_class_homo_counts.push_back(current_maj_class_homo_cnt);
        class_homo_coverages.push_back(current_class_homo_cov);
        class_coverages.push_back(current_class_cov);

        // Simplification loop
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
            auto [ct, homo_ct] = computeStats(maj_class_points, valleys_per_class, class_homo_points, class_all_points);
            remainingct.push_back(ct);
            remaininghomoct.push_back(homo_ct);
            
            // Compute coverages for this step
            for (uint32_t c = 0; c < num_classes; c++) {
                current_maj_class_homo_cnt[c] = valleys_per_class[c];
                if (class_totals[c] > 0) {
                    current_maj_class_homo_cov[c] = static_cast<double>(maj_class_points[c]) / static_cast<double>(class_totals[c]);
                    current_class_homo_cov[c] = static_cast<double>(class_homo_points[c]) / static_cast<double>(class_totals[c]);
                    current_class_cov[c] = static_cast<double>(class_all_points[c]) / static_cast<double>(class_totals[c]);
                } else {
                    current_maj_class_homo_cov[c] = 0.0;
                    current_class_homo_cov[c] = 0.0;
                    current_class_cov[c] = 0.0;
                }
            }
            maj_class_homo_coverages.push_back(current_maj_class_homo_cov);
            maj_class_homo_counts.push_back(current_maj_class_homo_cnt);
            class_homo_coverages.push_back(current_class_homo_cov);
            class_coverages.push_back(current_class_cov);
        }
    }

}  // namespace contourtree
