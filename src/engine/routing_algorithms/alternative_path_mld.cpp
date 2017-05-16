#include "engine/routing_algorithms/alternative_path.hpp"
#include "engine/routing_algorithms/routing_base_mld.hpp"

#include "util/integer_range.hpp"

#include <boost/assert.hpp>

#include <algorithm>
#include <iterator>
#include <memory>
#include <unordered_set>
#include <vector>

// TODO: debug, remove
#include <iostream>

namespace osrm
{
namespace engine
{
namespace routing_algorithms
{

// Unqualified calls below are from the mld namespace.
// This alternative implementation works only for mld.
using namespace mld;

// Alternative Routes for MLD.
//
// Start search from s and continue "for a while" when t was found. Save all vertices.
// Start search from t and continue "for a while" when s was found. Save all vertices.
// Intersect both vertex sets: these are the candidate vertices.
// For all candidate vertices c a (potentially arbitrarily bad) alternative route is (s, c, t).
// Apply heuristic to evaluate alternative route based on stretch, overlap, how reasonable it is.
//
// For MLD specifically we can pull off some tricks to make evaluating alternatives fast:
//   Only consider (s, c, t) with c border vertex: re-use MLD search steps.
//   Add meta data to border vertices: consider (s, c, t) only when c is e.g. on a highway.
//   Prune based on vertex cell id
//
// https://github.com/Project-OSRM/osrm-backend/issues/3905
InternalManyRoutesResult
alternativePathSearch(SearchEngineData<Algorithm> &search_engine_data,
                      const datafacade::ContiguousInternalMemoryDataFacade<Algorithm> &facade,
                      const PhantomNodes &phantom_node_pair)
{
    const auto &partition = facade.GetMultiLevelPartition();

    search_engine_data.InitializeOrClearFirstThreadLocalStorage(facade.GetNumberOfNodes());

    auto &forward_heap = *search_engine_data.forward_heap_1;
    auto &reverse_heap = *search_engine_data.reverse_heap_1;

    insertNodesInHeaps(forward_heap, reverse_heap, phantom_node_pair);

    //
    // Save nodes in the forward and backward search space overlap as candidates
    //

    NodeID middle = SPECIAL_NODEID;

    EdgeWeight forward_heap_min = forward_heap.MinKey();
    EdgeWeight reverse_heap_min = reverse_heap.MinKey();

    EdgeWeight path_weight = INVALID_EDGE_WEIGHT;
    EdgeWeight overlap_weight = INVALID_EDGE_WEIGHT;

    const auto overlap_factor = 1.66;

    const auto force_loop_forward = needsLoopForward(phantom_node_pair);
    const auto force_loop_backward = needsLoopBackwards(phantom_node_pair);

    std::vector<NodeID> candidates;

    // The structure below is a bit weird - here's why: we want to re-use the MLD routingStep for
    // stepping our search space from s and from t. We don't know how far to overlap until we have
    // the shortest path. Once we have the shortest path we can use its weight to terminate when
    // we're over factor * weight. We set the weight for the routingStep to INVALID_EDGE_WEIGHT
    // so that the routingStep will continue even after we reached the shortest path upper bound.

    while (forward_heap.Size() + reverse_heap.Size() > 0)
    {
        if (path_weight != INVALID_EDGE_WEIGHT)
            overlap_weight = path_weight * overlap_factor;

        // Termination criteria - when we have a shortest path this will guarantee for our overlap.
        const auto keep_going = forward_heap_min + reverse_heap_min < overlap_weight;

        if (!keep_going)
            break;

        // Force forward step to not break early when we reached the middle, continue for overlap.
        overlap_weight = INVALID_EDGE_WEIGHT;

        if (!forward_heap.Empty())
        {
            routingStep<FORWARD_DIRECTION>(facade,
                                           forward_heap,
                                           reverse_heap,
                                           middle,
                                           overlap_weight,
                                           force_loop_forward,
                                           force_loop_backward,
                                           phantom_node_pair);

            if (!forward_heap.Empty())
                forward_heap_min = forward_heap.MinKey();

            if (middle != SPECIAL_NODEID)
            {
                candidates.push_back(middle);
            }
        }

        // Adjusting upper bound for forward search
        path_weight = std::min(path_weight, overlap_weight);
        // Force reverse step to not break early when we reached the middle, continue for overlap.
        overlap_weight = INVALID_EDGE_WEIGHT;

        if (!reverse_heap.Empty())
        {
            routingStep<REVERSE_DIRECTION>(facade,
                                           reverse_heap,
                                           forward_heap,
                                           middle,
                                           overlap_weight,
                                           force_loop_forward,
                                           force_loop_backward,
                                           phantom_node_pair);

            if (!reverse_heap.Empty())
                reverse_heap_min = reverse_heap.MinKey();

            if (middle != SPECIAL_NODEID)
            {
                candidates.push_back(middle);
            }
        }

        // Adjusting upper bound for reverse search
        path_weight = std::min(path_weight, overlap_weight);
    }

    //
    // Todo: filter via candidate nodes with heuristics
    //

    std::cout << ">>> number of candidates: " << candidates.size() << std::endl;

    std::sort(begin(candidates), end(candidates));
    auto it = std::unique(begin(candidates), end(candidates));
    candidates.erase(it, end(candidates));

    std::cout << ">>> number of unique candidates: " << candidates.size() << std::endl;

    const auto route_found = path_weight != INVALID_EDGE_WEIGHT && middle != SPECIAL_NODEID;

    if (!route_found)
        return InternalManyRoutesResult{};

    //
    // Reconstruct paths
    //

    // Get packed path as edges {from node ID, to node ID, from_clique_arc}
    auto packed_path = retrievePackedPathFromHeap(forward_heap, reverse_heap, middle);

    // Beware the edge case when start, middle, end are all the same.
    // In this case we return a single node, no edges. We also don't unpack.
    const NodeID source_node = !packed_path.empty() ? std::get<0>(packed_path.front()) : middle;

    //
    // Todo: dup. code with mld::search except for level entry: we run a slight mld::search
    //       adaption here and then dispatch to mld::search for recursively descending down.
    //

    std::vector<NodeID> unpacked_nodes;
    std::vector<EdgeID> unpacked_edges;
    unpacked_nodes.reserve(packed_path.size());
    unpacked_edges.reserve(packed_path.size());

    unpacked_nodes.push_back(source_node);

    for (auto const &packed_edge : packed_path)
    {
        NodeID source, target;
        bool overlay_edge;
        std::tie(source, target, overlay_edge) = packed_edge;
        if (!overlay_edge)
        { // a base graph edge
            unpacked_nodes.push_back(target);
            unpacked_edges.push_back(facade.FindEdge(source, target));
        }
        else
        { // an overlay graph edge
            LevelID level = getNodeQueryLevel(partition, source, phantom_node_pair); // XXX
            CellID parent_cell_id = partition.GetCell(level, source);
            BOOST_ASSERT(parent_cell_id == partition.GetCell(level, target));

            LevelID sublevel = level - 1;

            // Here heaps can be reused, let's go deeper!
            forward_heap.Clear();
            reverse_heap.Clear();
            forward_heap.Insert(source, 0, {source});
            reverse_heap.Insert(target, 0, {target});

            // TODO: when structured bindings will be allowed change to
            // auto [subpath_weight, subpath_source, subpath_target, subpath] = ...
            EdgeWeight subpath_weight;
            std::vector<NodeID> subpath_nodes;
            std::vector<EdgeID> subpath_edges;
            std::tie(subpath_weight, subpath_nodes, subpath_edges) = search(search_engine_data,
                                                                            facade,
                                                                            forward_heap,
                                                                            reverse_heap,
                                                                            force_loop_forward,
                                                                            force_loop_backward,
                                                                            INVALID_EDGE_WEIGHT,
                                                                            sublevel,
                                                                            parent_cell_id);
            BOOST_ASSERT(!subpath_edges.empty());
            BOOST_ASSERT(subpath_nodes.size() > 1);
            BOOST_ASSERT(subpath_nodes.front() == source);
            BOOST_ASSERT(subpath_nodes.back() == target);
            unpacked_nodes.insert(
                unpacked_nodes.end(), std::next(subpath_nodes.begin()), subpath_nodes.end());
            unpacked_edges.insert(unpacked_edges.end(), subpath_edges.begin(), subpath_edges.end());
        }
    }

    auto primary_route =
        extractRoute(facade, path_weight, phantom_node_pair, unpacked_nodes, unpacked_edges);

    //
    // Todo: do the same for alternative paths
    //

    std::vector<InternalRouteResult> routes;
    routes.push_back(std::move(primary_route));

    return InternalManyRoutesResult{std::move(routes)};
}

} // namespace routing_algorithms
} // namespace engine
} // namespace osrm}
