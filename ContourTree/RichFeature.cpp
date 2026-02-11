#include "RichFeature.hpp"
#include <cassert>
#include <iostream>
#include <algorithm>

namespace contourtree {

std::vector<RichFeature> computeRichFeatures(
    TopologicalFeatures& topo,
    int topk,
    float threshold,
    const std::vector<uint32_t>& partition,
    const std::vector<uint32_t>& labels,
    const std::vector<uint32_t>& preds,
    const std::vector<uint32_t>& class_sizes
) {
    // Get basic features from TopologicalFeatures
    int topk_copy = topk;
    std::vector<Feature> basic_features = topo.getArcFeatures(topk_copy, threshold);
    
    std::vector<RichFeature> rich_features;
    rich_features.reserve(basic_features.size());
    
    // Build arc map: arc_id -> feature index
    std::map<uint32_t, size_t> arc_map;
    
    for (size_t i = 0; i < basic_features.size(); i++) {
        RichFeature rf;
        rf.id = i;
        rf.from = basic_features[i].from;
        rf.to = basic_features[i].to;
        rf.arcs = basic_features[i].arcs;
        
        // Get function values and types
        uint32_t from_corrected = topo.ctdata.nodeMap[rf.from];
        uint32_t to_corrected = topo.ctdata.nodeMap[rf.to];
        
        rf.fn_from = topo.ctdata.fnVals[from_corrected];
        rf.fn_to = topo.ctdata.fnVals[to_corrected];
        rf.type_from = topo.ctdata.type[from_corrected];
        rf.type_to = topo.ctdata.type[to_corrected];
        
        rf.persistence = rf.fn_to - rf.fn_from;
        
        // Build arc map
        for (uint32_t arc_id : rf.arcs) {
            arc_map[arc_id] = i;
        }
        
        rich_features.push_back(rf);
    }
    
    // Populate members and class counts from partition
    for (size_t i = 0; i < partition.size(); i++) {
        uint32_t arc_id = partition[i];
        
        auto it = arc_map.find(arc_id);
        if (it == arc_map.end()) {
            std::cerr << "Warning: Arc " << arc_id << " not found in feature map!" << std::endl;
            continue;
        }
        
        size_t feat_idx = it->second;
        RichFeature& feat = rich_features[feat_idx];
        
        feat.members.insert(i);
        feat.size++;
        
        uint32_t label = labels[i];
        feat.class_counts[label]++;
    }
    
    // Compute majority class and homogeneity for each feature
    for (RichFeature& feat : rich_features) {
        if (feat.class_counts.empty()) {
            // No points in this feature, use from label as default
            uint32_t from_corrected = topo.ctdata.nodeMap[feat.from];
            feat.majority_class = (from_corrected < labels.size()) ? labels[from_corrected] : 0;
            feat.major_class_size = 0;
            feat.homogeneity = 0.0;
        } else {
            // Find majority class
            uint32_t max_count = 0;
            for (const auto& pair : feat.class_counts) {
                if (pair.second > max_count) {
                    max_count = pair.second;
                    feat.majority_class = pair.first;
                }
            }
            feat.major_class_size = max_count;
            feat.homogeneity = feat.size > 0 ? static_cast<float>(max_count) / feat.size : 0.0f;
        }
        
        // Compute class proportions (within feature) and class coverage (within dataset)
        for (const auto& pair : feat.class_counts) {
            uint32_t class_idx = pair.first;
            uint32_t count = pair.second;
            
            feat.class_proportions[class_idx] = feat.size > 0 ? 
                static_cast<float>(count) / feat.size : 0.0f;
            
            uint32_t total_class_size = class_sizes[class_idx];
            feat.class_coverage[class_idx] = total_class_size > 0 ? 
                static_cast<float>(count) / total_class_size : 0.0f;
        }
    }
    
    // Compute prediction-based metadata
    bool has_preds = !preds.empty() && preds.size() == labels.size();
    
    for (RichFeature& feat : rich_features) {
        for (uint32_t dp : feat.members) {
            uint32_t true_label = labels[dp];
            uint32_t pred_label = has_preds ? preds[dp] : labels[dp];
            
            feat.pred_class_counts[pred_label]++;
            feat.confusion[{true_label, pred_label}]++;
            
            if (true_label == pred_label) {
                feat.pred_correct++;
            } else {
                feat.pred_incorrect++;
            }
        }
        
        feat.pred_accuracy = feat.size > 0 ? 
            static_cast<float>(feat.pred_correct) / feat.size : 0.0f;
    }
    
    return rich_features;
}

std::vector<size_t> filterFeatures(
    const std::vector<RichFeature>& features,
    const FilterCriteria& criteria
) {
    std::vector<size_t> filtered_indices;
    filtered_indices.reserve(features.size());
    
    for (size_t i = 0; i < features.size(); i++) {
        const RichFeature& feat = features[i];
        
        // Check size constraints
        if (feat.size < criteria.min_size || feat.size > criteria.max_size) {
            continue;
        }
        
        // Check fn_from constraints
        if (feat.fn_from < criteria.min_fn_from || feat.fn_from > criteria.max_fn_from) {
            continue;
        }
        
        // Check homogeneity constraints
        if (feat.homogeneity < criteria.min_homogeneity || feat.homogeneity > criteria.max_homogeneity) {
            continue;
        }
        
        // Check type constraints
        if (!criteria.allowed_types.empty()) {
            std::pair<char, char> feat_type = {feat.type_from, feat.type_to};
            if (criteria.allowed_types.find(feat_type) == criteria.allowed_types.end()) {
                continue;
            }
        }
        
        filtered_indices.push_back(i);
    }
    
    return filtered_indices;
}

}  // namespace contourtree
