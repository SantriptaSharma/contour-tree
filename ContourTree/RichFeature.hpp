#ifndef RICHFEATURE_HPP
#define RICHFEATURE_HPP

#include "TopologicalFeatures.hpp"
#include <map>
#include <set>
#include <vector>
#include <string>
#include <limits>

namespace contourtree {

struct RichFeature {
    // Basic feature info
    uint32_t id;
    uint32_t from;
    uint32_t to;
    float fn_from;
    float fn_to;
    char type_from;
    char type_to;
    float persistence;
    std::vector<uint32_t> arcs;
    
    // Label-based metadata
    std::set<uint32_t> members;
    uint32_t size;
    std::map<uint32_t, uint32_t> class_counts;
    std::map<uint32_t, float> class_proportions;  // proportion of each class within this feature (class_count / feature_size)
    std::map<uint32_t, float> class_coverage;     // coverage of each class in dataset (class_count / total_class_size)
    uint32_t majority_class;
    uint32_t major_class_size;
    float homogeneity;
    
    // Prediction-based metadata
    std::map<uint32_t, uint32_t> pred_class_counts;
    std::map<std::pair<uint32_t, uint32_t>, uint32_t> confusion;
    uint32_t pred_correct;
    uint32_t pred_incorrect;
    float pred_accuracy;
    
    RichFeature() : id(0), from(0), to(0), fn_from(0), fn_to(0), 
                    type_from(0), type_to(0), persistence(0),
                    size(0), majority_class(0), major_class_size(0), homogeneity(0),
                    pred_correct(0), pred_incorrect(0), pred_accuracy(0) {}
};

// Compute rich features with all metadata
std::vector<RichFeature> computeRichFeatures(
    TopologicalFeatures& topo,
    int topk,
    float threshold,
    const std::vector<uint32_t>& partition,
    const std::vector<uint32_t>& labels,
    const std::vector<uint32_t>& preds,
    const std::vector<uint32_t>& class_sizes = std::vector<uint32_t>()
);

struct FilterCriteria {
    uint32_t min_size;
    uint32_t max_size;
    float min_fn_from;
    float max_fn_from;
    float min_homogeneity;
    float max_homogeneity;
    std::set<std::pair<char, char>> allowed_types;
    
    FilterCriteria() : min_size(0), max_size(UINT32_MAX),
                       min_fn_from(0.0f), max_fn_from(std::numeric_limits<float>::max()),
                       min_homogeneity(0.0f), max_homogeneity(1.0f) {}
};

// Filter features based on criteria, returns indices of matching features
std::vector<size_t> filterFeatures(
    const std::vector<RichFeature>& features,
    const FilterCriteria& criteria
);

}  // namespace contourtree

#endif  // RICHFEATURE_HPP
