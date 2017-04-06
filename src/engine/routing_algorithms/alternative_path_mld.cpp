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
    search_engine_data.InitializeOrClearFirstThreadLocalStorage(facade.GetNumberOfNodes());

    auto &forward_heap = *search_engine_data.forward_heap_1;
    auto &reverse_heap = *search_engine_data.reverse_heap_1;

    insertNodesInHeaps(forward_heap, reverse_heap, phantom_node_pair);

    // TODO: use for pruning
    // const auto &partition = facade.GetMultiLevelPartition();
    // const auto &cells = facade.GetCellStorage();

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

        // Force routeSteps to not break early when we reached the middle, continue for overlap.
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

        path_weight = overlap_weight;
    }

    //
    // Filter via candidate nodes with heuristics
    //

    std::cout << ">>> number of candidates: " << candidates.size() << std::endl;

    std::sort(begin(candidates), end(candidates));
    auto it = std::unique(begin(candidates), end(candidates));
    candidates.erase(it, end(candidates));

    std::cout << ">>> number of unique candidates: " << candidates.size() << std::endl;

    const auto route_found = path_weight != INVALID_EDGE_WEIGHT && middle != SPECIAL_NODEID;

    //
    // Reconstruct routes (if any)
    //

    InternalRouteResult primary_route;
    primary_route.segment_end_coordinates = {phantom_node_pair};

    if (!route_found)
        return InternalManyRoutesResult{std::move(primary_route)};

    std::vector<InternalRouteResult> routes;
    routes.push_back(std::move(primary_route));

    return InternalManyRoutesResult{std::move(routes)};
}

} // namespace routing_algorithms
} // namespace engine
} // namespace osrm}
