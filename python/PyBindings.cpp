#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include <memory>
#include "GraphScalarFunction.hpp"
#include "constants.h"
#include "TriMesh.hpp"
#include "Grid3D.hpp"
#include "MergeTree.hpp"
#include "ContourTreeData.hpp"
#include "SimplifyCT.hpp"
#include "SimFunction.hpp"
#include "Persistence.hpp"
#include "HyperVolume.hpp"
#include "TopologicalFeatures.hpp"
#include "LayoutCT.hpp"
#include "RichFeature.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyct, m) {
    m.doc() = "Python bindings for Contour Tree library";
    
    // Expose constants
    m.attr("REGULAR") = contourtree::REGULAR;
    m.attr("MINIMUM") = contourtree::MINIMUM;
    m.attr("MAXIMUM") = contourtree::MAXIMUM;
    m.attr("SADDLE") = contourtree::SADDLE;
    
    // Expose TreeType enum
    py::enum_<contourtree::TreeType>(m, "TreeType")
        .value("TypeJoinTree", contourtree::TypeJoinTree)
        .value("TypeSplitTree", contourtree::TypeSplitTree)
        .value("TypeContourTree", contourtree::TypeContourTree);
    
    // Expose ScalarFunction base class
    py::class_<contourtree::ScalarFunction, std::shared_ptr<contourtree::ScalarFunction>>(m, "ScalarFunction")
        .def("getMaxDegree", &contourtree::ScalarFunction::getMaxDegree)
        .def("getVertexCount", &contourtree::ScalarFunction::getVertexCount)
        .def("getStar", [](contourtree::ScalarFunction& self, int64_t v) {
            std::vector<int64_t> star(self.getMaxDegree());
            int count = self.getStar(v, star);
            star.resize(count);
            return star;
        }, "Get star neighbors of vertex v")
        .def("lessThan", &contourtree::ScalarFunction::lessThan)
        .def("getFunctionValue", &contourtree::ScalarFunction::getFunctionValue);
    
    // Expose GraphScalarFunction class
    py::class_<contourtree::GraphScalarFunction, std::shared_ptr<contourtree::GraphScalarFunction>, contourtree::ScalarFunction>(m, "GraphScalarFunction")
        .def(py::init<>())
        .def("getVertexCount", &contourtree::GraphScalarFunction::getVertexCount)
        .def("getFunctionValue", &contourtree::GraphScalarFunction::getFunctionValue)
        .def("initialize", &contourtree::GraphScalarFunction::initialize,
             "Initialize graph with specified number of nodes")
        .def("loadGraph", &contourtree::GraphScalarFunction::loadGraph,
             "Load graph from edge file")
        .def("loadGraphFromAdjList", &contourtree::GraphScalarFunction::loadGraphFromAdjList, 
            "Load graph from an adjacency list")
        .def("updateFnValues", &contourtree::GraphScalarFunction::updateFnValues,
             "Update function values");

    // Expose TriMesh class
    py::class_<contourtree::TriMesh, std::shared_ptr<contourtree::TriMesh>, contourtree::ScalarFunction>(m, "TriMesh")
        .def(py::init<>())
        .def("getVertexCount", &contourtree::TriMesh::getVertexCount)
        .def("getFunctionValue", &contourtree::TriMesh::getFunctionValue)
        .def("loadData", &contourtree::TriMesh::loadData, "Load mesh data from file");

    // Grid3D (template) bindings: instantiate for float and double
    using Grid3DUint8 = contourtree::Grid3D<uint8_t>;
    py::class_<Grid3DUint8, std::shared_ptr<Grid3DUint8>, contourtree::ScalarFunction>(m, "Grid3DUint8")
        .def(py::init<int,int,int>(), py::arg("resx"), py::arg("resy"), py::arg("resz"))
        .def("getVertexCount", &Grid3DUint8::getVertexCount)
        .def("getFunctionValue", &Grid3DUint8::getFunctionValue)
        .def("loadGrid", &Grid3DUint8::loadGrid, py::arg("fileName"), "Load scalar field from binary file")
        .def_readonly("dimx", &Grid3DUint8::dimx)
        .def_readonly("dimy", &Grid3DUint8::dimy)
        .def_readonly("dimz", &Grid3DUint8::dimz);

    using Grid3DFloat = contourtree::Grid3D<float>;
    py::class_<Grid3DFloat, std::shared_ptr<Grid3DFloat>, contourtree::ScalarFunction>(m, "Grid3DFloat")
        .def(py::init<int,int,int>(), py::arg("resx"), py::arg("resy"), py::arg("resz"))
        .def("getVertexCount", &Grid3DFloat::getVertexCount)
        .def("getFunctionValue", &Grid3DFloat::getFunctionValue)
        .def("loadGrid", &Grid3DFloat::loadGrid, py::arg("fileName"), "Load scalar field from binary file")
        .def_readonly("dimx", &Grid3DFloat::dimx)
        .def_readonly("dimy", &Grid3DFloat::dimy)
        .def_readonly("dimz", &Grid3DFloat::dimz);

    using Grid3DDouble = contourtree::Grid3D<double>;
    py::class_<Grid3DDouble, std::shared_ptr<Grid3DDouble>, contourtree::ScalarFunction>(m, "Grid3DDouble")
        .def(py::init<int,int,int>(), py::arg("resx"), py::arg("resy"), py::arg("resz"))
        .def("getVertexCount", &Grid3DDouble::getVertexCount)
        .def("getFunctionValue", &Grid3DDouble::getFunctionValue)
        .def("loadGrid", &Grid3DDouble::loadGrid, py::arg("fileName"), "Load scalar field from binary file")
        .def_readonly("dimx", &Grid3DDouble::dimx)
        .def_readonly("dimy", &Grid3DDouble::dimy)
        .def_readonly("dimz", &Grid3DDouble::dimz);

    // Expose MergeTree class
    py::class_<contourtree::MergeTree, std::shared_ptr<contourtree::MergeTree>>(m, "MergeTree")
        .def(py::init<>())
        .def("computeTree", &contourtree::MergeTree::computeTree,
                py::arg("data"), py::arg("treeType"),
                "Compute the specified tree using the provided ScalarFunction")
        .def("output", &contourtree::MergeTree::output,
                py::arg("fileName"), py::arg("treeType"),
                "Write tree data to files with the given base name");

    py::class_<contourtree::Node, std::shared_ptr<contourtree::Node>>(m, "Node")
        .def(py::init<>())
        .def_readonly("next", &contourtree::Node::next)
        .def_readonly("prev", &contourtree::Node::prev);

    py::class_<contourtree::Arc, std::shared_ptr<contourtree::Arc>>(m, "Arc")
        .def(py::init<>())
        .def_readonly("frm", &contourtree::Arc::from)
        .def_readonly("to", &contourtree::Arc::to)
        .def_readonly("id", &contourtree::Arc::id);

    // Expose ContourTreeData class (minimal)
    py::class_<contourtree::ContourTreeData, std::shared_ptr<contourtree::ContourTreeData>>(m, "ContourTreeData")
        .def(py::init<>())
        .def("loadBinFile", &contourtree::ContourTreeData::loadBinFile, py::arg("filename"), "Load data from a binary file")
        .def_readonly("fnVals", &contourtree::ContourTreeData::fnVals)
        .def_readonly("type", &contourtree::ContourTreeData::type)
        .def_readonly("nodes", &contourtree::ContourTreeData::nodes)
        .def_readonly("nodeMap", &contourtree::ContourTreeData::nodeMap)
        .def_readonly("arcs", &contourtree::ContourTreeData::arcs);

    // Base class for similarity functions (no Python overrides needed right now)
    py::class_<contourtree::SimFunction, std::shared_ptr<contourtree::SimFunction>>(m, "SimFunction");

    // Persistence and HyperVolume constructors
    py::class_<contourtree::Persistence, std::shared_ptr<contourtree::Persistence>, contourtree::SimFunction>(m, "Persistence")
        .def(py::init<const contourtree::ContourTreeData&>(), py::arg("ctData"));

    py::class_<contourtree::HyperVolume, std::shared_ptr<contourtree::HyperVolume>, contourtree::SimFunction>(m, "HyperVolume")
        .def(py::init<const contourtree::ContourTreeData&, std::string>(),
             py::arg("ctData"), py::arg("partFile"));

    // Expose SimplifyCT class
    py::class_<contourtree::SimplifyCT, std::shared_ptr<contourtree::SimplifyCT>>(m, "SimplifyCT")
        .def(py::init<>())
        .def("setInput", &contourtree::SimplifyCT::setInput, py::arg("data"))
        .def(
            "simplify",
            static_cast<void (contourtree::SimplifyCT::*)(contourtree::SimFunction*)>(&contourtree::SimplifyCT::simplify),
            py::arg("simFn"),
            "Simplify using Persistence or HyperVolume")
        .def("outputOrder", &contourtree::SimplifyCT::outputOrder, py::arg("fileName"), py::arg("normalize"), "Write the branch removal order to disk")
        .def("getSimplificationPlot",  [](contourtree::SimplifyCT& self, const std::vector<uint32_t>& order, const std::vector<float>& wts) {
            std::vector<float> fns;
            std::vector<int32_t> minct, maxct;

            self.getSimplificationPlot(order, wts, fns, minct, maxct);
        
            return py::make_tuple(fns, minct, maxct);
        }, py::arg("order"), py::arg("wts"),
           "Get function values and counts of minima/maxima at each of them for simplification plot")
        .def("getFilteredSimplificationPlot",  [](contourtree::SimplifyCT& self, const std::vector<uint32_t>& order, const std::vector<float>& wts,
                                                  float minfnstart, float maxfnstart, float minfnend, float maxfnend, const std::vector<char> &types) {
            std::vector<float> fns;
            std::vector<int32_t> filteredct;

            std::set<char> type_set(types.begin(), types.end());

            self.getFilteredSimplificationPlot(order, wts, fns, minfnstart, maxfnstart, minfnend, maxfnend, type_set, filteredct);
        
            return py::make_tuple(fns, filteredct);
        }, py::arg("order"), py::arg("wts"), py::arg("minfnstart"), py::arg("maxfnstart"),
           py::arg("minfnend"), py::arg("maxfnend"), py::arg("types"),
           "Get function values and counts of filtered extrema at each of them for simplification plot")
        .def("getHomoValleyPlot",  [](contourtree::SimplifyCT& self, const std::vector<uint32_t>& order, const std::vector<float>& wts,
                                      const std::vector<uint32_t> &labels, float homogeneity_threshold, const std::vector<uint32_t> &partition) {
            std::vector<float> fns;
            std::vector<int32_t> remainingct;

            self.getHomoValleyPlot(order, wts, fns, remainingct, labels, homogeneity_threshold, partition);
        
            return py::make_tuple(fns, remainingct);
        }, py::arg("order"), py::arg("wts"), py::arg("labels"), py::arg("homogeneity_threshold"), py::arg("partition"),
           "Get function values, counts of remaining homogeneous valleys, class coverages, and class valley counts at each simplification step")
        .def("getHomoValleyPlotPlusCoverages",  [](contourtree::SimplifyCT& self, const std::vector<uint32_t>& order, const std::vector<float>& wts,
                                      const std::vector<uint32_t> &labels, float homogeneity_threshold, const std::vector<uint32_t> &partition) {
            std::vector<float> fns;
            std::vector<int32_t> remainingct;
            std::vector<int32_t> remaininghomoct;
            std::vector<std::vector<double>> maj_class_homo_coverages;
            std::vector<std::vector<int>> maj_class_homo_counts;
            std::vector<std::vector<double>> class_homo_coverages;
            std::vector<std::vector<double>> class_coverages;

            self.getHomoValleyPlotPlusCoverages(order, wts, fns, remainingct, remaininghomoct, labels, homogeneity_threshold, partition, 
                                                maj_class_homo_coverages, maj_class_homo_counts, class_homo_coverages, class_coverages);
        
            return py::make_tuple(fns, remainingct, remaininghomoct, maj_class_homo_coverages, maj_class_homo_counts, class_homo_coverages, class_coverages);
        }, py::arg("order"), py::arg("wts"), py::arg("labels"), py::arg("homogeneity_threshold"), py::arg("partition"),
           "Get function values, counts, majority class homogeneous coverages, majority class homogeneous counts, class homogeneous coverages, and all class coverages");
        // .def("getFilteredSimplificationPlotHomogeneity",  [](contourtree::SimplifyCT& self, const std::vector<uint32_t>& order, const std::vector<float>& wts,
        //                                           float minfnstart, float maxfnstart, float minfnend, float maxfnend, const std::vector<char> &types,
        //                                           const std::vector<uint32_t> &labels, float homogeneity_threshold, const std::vector<uint32_t> &partition) {
        //     std::vector<float> fns;
        //     std::vector<int32_t> filteredct;
        //     std::vector<uint32_t> branch_total_sizes;
        //     std::vector<uint32_t> branch_majority_labels;
        //     std::vector<uint32_t> branch_majority_sizes;
        //     std::vector<bool> branch_was_homogeneous;
        //     std::vector<bool> caused_homogeneous_destruction;

        //     std::set<char> type_set(types.begin(), types.end());

        //     self.getFilteredSimplificationPlotHomogeneity(order, wts, fns, minfnstart, maxfnstart, minfnend, maxfnend, type_set, filteredct, labels, homogeneity_threshold, partition,
        //                                                   branch_total_sizes, branch_majority_labels, branch_majority_sizes, branch_was_homogeneous, caused_homogeneous_destruction);
        
        //     return py::make_tuple(fns, filteredct, branch_total_sizes, branch_majority_labels, branch_majority_sizes, branch_was_homogeneous, caused_homogeneous_destruction);
        // }, py::arg("order"), py::arg("wts"), py::arg("minfnstart"), py::arg("maxfnstart"),
        //    py::arg("minfnend"), py::arg("maxfnend"), py::arg("types"), py::arg("labels"), py::arg("homogeneity_threshold"), py::arg("partition"),
        //    "Get function values, counts, and metadata (total_size, majority_label, majority_size, was_homogeneous, caused_homogeneous_destruction) for each removed branch");

    // Expose TopologicalFeatures::Feature class
    py::class_<contourtree::Feature, std::shared_ptr<contourtree::Feature>>(m, "Feature")
        .def(py::init<>())
        .def_readwrite("arcs", &contourtree::Feature::arcs)
        .def_readwrite("frm", &contourtree::Feature::from)
        .def_readwrite("to", &contourtree::Feature::to);

    // Expose RichFeature class
    py::class_<contourtree::RichFeature, std::shared_ptr<contourtree::RichFeature>>(m, "RichFeature")
        .def(py::init<>())
        .def_readwrite("id", &contourtree::RichFeature::id)
        .def_readwrite("frm", &contourtree::RichFeature::from)
        .def_readwrite("to", &contourtree::RichFeature::to)
        .def_readwrite("fn_frm", &contourtree::RichFeature::fn_from)
        .def_readwrite("fn_to", &contourtree::RichFeature::fn_to)
        .def_readwrite("type_frm", &contourtree::RichFeature::type_from)
        .def_readwrite("type_to", &contourtree::RichFeature::type_to)
        .def_readwrite("pers", &contourtree::RichFeature::persistence)
        .def_readwrite("arcs", &contourtree::RichFeature::arcs)
        .def_readwrite("members", &contourtree::RichFeature::members)
        .def_readwrite("size", &contourtree::RichFeature::size)
        .def_readwrite("class_counts", &contourtree::RichFeature::class_counts)
        .def_readwrite("class_proportions", &contourtree::RichFeature::class_proportions)
        .def_readwrite("class_coverage", &contourtree::RichFeature::class_coverage)
        .def_readwrite("majority_class", &contourtree::RichFeature::majority_class)
        .def_readwrite("major_class_size", &contourtree::RichFeature::major_class_size)
        .def_readwrite("homogeneity", &contourtree::RichFeature::homogeneity)
        .def_readwrite("pred_class_counts", &contourtree::RichFeature::pred_class_counts)
        .def_readwrite("confusion", &contourtree::RichFeature::confusion)
        .def_readwrite("pred_correct", &contourtree::RichFeature::pred_correct)
        .def_readwrite("pred_incorrect", &contourtree::RichFeature::pred_incorrect)
        .def_readwrite("pred_accuracy", &contourtree::RichFeature::pred_accuracy);

    // Expose TopologicalFeatures class
    py::class_<contourtree::TopologicalFeatures, std::shared_ptr<contourtree::TopologicalFeatures>>(m, "TopologicalFeatures")
        .def(py::init<>())
        .def("loadData", &contourtree::TopologicalFeatures::loadData, py::arg("filename"))
        .def("getArcFeatures", [](contourtree::TopologicalFeatures& self, int topk, float th) {
            int topk_copy = topk;
            auto features = self.getArcFeatures(topk_copy, th);
            return py::make_tuple(features, topk_copy);
        }, py::arg("topk"), py::arg("th") = 0.0f, "Returns tuple of (features, updated_topk)")
        .def("getPartitionedExtremaFeatures", [](contourtree::TopologicalFeatures& self, int topk, float th) {
            int topk_copy = topk;
            auto features = self.getPartitionedExtremaFeatures(topk_copy, th);
            return py::make_tuple(features, topk_copy);
        }, py::arg("topk"), py::arg("th") = 0.0f, "Returns tuple of (features, updated_topk)")
        .def_readonly("ctdata", &contourtree::TopologicalFeatures::ctdata);

    // Expose computeRichFeatures function
    m.def("computeRichFeatures", [](contourtree::TopologicalFeatures& topo, int topk, float threshold,
                                     const std::vector<uint32_t>& partition, const std::vector<uint32_t>& labels,
                                     const std::vector<uint32_t>& preds, const std::vector<uint32_t>& class_sizes) {
        return contourtree::computeRichFeatures(topo, topk, threshold, partition, labels, preds, class_sizes);
    }, py::arg("topo"), py::arg("topk"), py::arg("threshold"), py::arg("partition"), 
       py::arg("labels"), py::arg("preds") = std::vector<uint32_t>(),
       py::arg("class_sizes") = std::vector<uint32_t>(),
       "Compute rich features with label and prediction metadata");

    // Expose FilterCriteria struct
    py::class_<contourtree::FilterCriteria>(m, "FilterCriteria")
        .def(py::init<>())
        .def_readwrite("min_size", &contourtree::FilterCriteria::min_size)
        .def_readwrite("max_size", &contourtree::FilterCriteria::max_size)
        .def_readwrite("min_fn_from", &contourtree::FilterCriteria::min_fn_from)
        .def_readwrite("max_fn_from", &contourtree::FilterCriteria::max_fn_from)
        .def_readwrite("min_homogeneity", &contourtree::FilterCriteria::min_homogeneity)
        .def_readwrite("max_homogeneity", &contourtree::FilterCriteria::max_homogeneity)
        .def_readwrite("allowed_types", &contourtree::FilterCriteria::allowed_types);

    // Expose filterFeatures function
    m.def("filterFeatures", &contourtree::filterFeatures,
          py::arg("features"), py::arg("criteria"),
          "Filter features based on criteria, returns indices of matching features");

    // Expose Point struct
    py::class_<contourtree::Point, std::shared_ptr<contourtree::Point>>(m, "Point")
        .def(py::init<>())
        .def_readwrite("x", &contourtree::Point::x)
        .def_readwrite("y", &contourtree::Point::y)
        .def_readwrite("z", &contourtree::Point::z);

    // Expose LayoutCT class
    py::class_<contourtree::LayoutCT, std::shared_ptr<contourtree::LayoutCT>>(m, "LayoutCT")
        .def(py::init<contourtree::TopologicalFeatures*>(), py::arg("tf"))
        .def("layoutTree", &contourtree::LayoutCT::layoutTree, py::arg("simplifiedCount"))
        .def("getNodeLocations", &contourtree::LayoutCT::getNodeLocations);

    // Expose SaveLayoutToOFF function
    m.def("SaveLayoutToOFF", &contourtree::SaveLayoutToOFF, py::arg("dataName"), py::arg("topk"), py::arg("thresh"),
          "Save layout to OFF file");
}
