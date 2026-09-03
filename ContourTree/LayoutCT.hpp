#ifndef LAYOUTCT_HPP
#define LAYOUTCT_HPP

#include "TopologicalFeatures.hpp"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stdint.h>
namespace contourtree {

struct Point {
    float x,y,z;
};

struct BranchLocation {
    int id;
    float angle;
    int level;
    float x,z;

    BranchLocation() : level(-1) {
    }
};

class LayoutCT
{
public:
    LayoutCT(TopologicalFeatures* tf);

    void layoutTree(int simplifiedCount);
    std::unordered_map<uint32_t, Point> getNodeLocations();

protected:
    int countLeaves(int brNo, int level);
    void assignAngles(int brNo);
    void assignLocations();
    void assignBranchIds();

public:
    SimplifyCT *simct;
    std::vector<uint32_t> &order;

    std::vector<Point> nodeLocs;
    std::vector<BranchLocation> branches;

    int noNodes;
    int noBranches;
    std::vector<int> leafCount;
    std::vector<float> startAngle;
    int noLeaves;
    int maxLevel;
//    int lastIndex;
    std::vector<float> r;
    std::vector<int> brMap;
    std::unordered_set<int> brSet;
    std::unordered_map<int,int> invMap;
    std::unordered_map<uint32_t,uint32_t> branchIds;
    float layoutExtent;

    const float rootTwo;
};

void SaveLayoutToOFF(std::string dataName, int &topk, float thresh);

// Rich variant: also computes per-feature metadata (persistence, size, majority class,
// homogeneity, majority share) via RichFeature and appends it to each edge line in the
// OFF file, after the standard "2 v1 v2 val1 val2 type1 type2" columns, as:
// id persistence size majority_class major_class_size homogeneity majority_share from_node_id to_node_id majority_class_label
void SaveRichLayoutToOFF(std::string dataName, int &topk, float thresh,
                         const std::vector<uint32_t>& partition,
                         const std::vector<uint32_t>& labels,
                         const std::vector<uint32_t>& preds = std::vector<uint32_t>(),
                         const std::vector<uint32_t>& class_sizes = std::vector<uint32_t>(),
                         const std::vector<std::string>& class_labels = std::vector<std::string>());

} // namespace contourtree
#endif // LAYOUTCT_HPP