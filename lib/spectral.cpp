#ifndef MESH_SIMPLIFICATION_CPP
#define MESH_SIMPLIFICATION_CPP

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <cstdint> // for int64_t, uint32_t
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
// #include <glm/gtc/matrix_inverse.hpp> // Not strictly needed for this spectral version if not using QEM's matrix inversion

using namespace std;
using namespace glm;

const float COT_EPSILON = 1e-6f; // Epsilon for cotangent calculations and positive weight clamping
const int V_OPT_ITERATIONS = 5;  // Number of iterations for optimal vertex position
const float V_OPT_CONVERGENCE_THRESHOLD = 1e-5f;

// Edge structure definition
struct Edge
{
    int v1, v2;
    double cost;
    vec3 optimal_pos;
    bool operator>(const Edge &o) const { return cost > o.cost; }
};

// Global data
static vector<vec3> vertices;
static vector<ivec3> faces; // Current list of faces
// static vector<mat4> quadrics; // Removed for spectral method
static priority_queue<Edge, vector<Edge>, greater<Edge>> edge_queue;
static unordered_map<int64_t, Edge> edge_map;       // Using int64_t
static vector<unordered_set<int>> vertex_neighbors; // Adjacency list
static Edge best_edge;                              // For SimplifyStep

// Calculate unique edge key
inline int64_t edge_key(int a, int b)
{
    if (a > b)
        std::swap(a, b);
    return ((int64_t)a << 32) | (uint32_t)b;
}

// Robust cotangent of angle <apb (angle at p)
// v_pa = a - p, v_pb = b - p
float cotangent_robust(const vec3 &v_pa, const vec3 &v_pb)
{
    float dot_product = dot(v_pa, v_pb);
    float cross_product_mag = length(cross(v_pa, v_pb));

    if (cross_product_mag < COT_EPSILON)
    { // Collinear or nearly collinear
        // Check if angle is close to 0 or 180.
        // If dot_product is positive, angle is close to 0, cot is large positive.
        // If dot_product is negative, angle is close to 180, cot is large negative.
        // For weights, often clamped to positive.
        return (dot_product > 0) ? (1.0f / COT_EPSILON) : (-1.0f / COT_EPSILON); // Large magnitude
    }
    return dot_product / cross_product_mag;
}

// Find vertices that form triangles with edge (v_opt_pos, P_neighbor_pos) in the new configuration
// P_neighbor_idx is an index from the ring_indices
// v1_orig_idx, v2_orig_idx are the vertices of the edge being collapsed
// ring_indices contains all neighbors of v1_orig_idx or v2_orig_idx (excluding v1, v2 themselves)
std::vector<int> find_triangle_partners_for_new_edge(
    int P_neighbor_idx,
    int v1_orig_idx,
    int v2_orig_idx,
    const std::unordered_set<int> &ring_indices)
{
    std::vector<int> partners;
    std::unordered_set<int> distinct_partners; // To handle cases where a partner might be found via v1 and v2

    for (const auto &face : faces)
    {
        int face_verts[] = {face.x, face.y, face.z};
        bool has_P_neighbor = false;
        bool has_v1 = false;
        bool has_v2 = false;
        int other_v = -1;

        for (int v_idx : face_verts)
        {
            if (v_idx == P_neighbor_idx)
                has_P_neighbor = true;
            else if (v_idx == v1_orig_idx)
                has_v1 = true;
            else if (v_idx == v2_orig_idx)
                has_v2 = true;
            else
                other_v = v_idx;
        }

        if (has_P_neighbor)
        {
            if (has_v1 && !has_v2)
            { // Triangle (P_neighbor, v1_orig, other_v)
                // other_v must be in ring_indices to form new triangle (v_opt, P_neighbor, other_v)
                if (other_v != -1 && ring_indices.count(other_v))
                {
                    if (distinct_partners.find(other_v) == distinct_partners.end())
                    {
                        partners.push_back(other_v);
                        distinct_partners.insert(other_v);
                    }
                }
            }
            else if (has_v2 && !has_v1)
            { // Triangle (P_neighbor, v2_orig, other_v)
                if (other_v != -1 && ring_indices.count(other_v))
                {
                    if (distinct_partners.find(other_v) == distinct_partners.end())
                    {
                        partners.push_back(other_v);
                        distinct_partners.insert(other_v);
                    }
                }
            }
        }
        if (partners.size() >= 2)
            break; // Max 2 partners for a manifold edge
    }
    return partners;
}

// Calculate edge contraction optimal position and cost using a local Laplacian energy approach
pair<vec3, double> compute_edge_cost(int v1_idx, int v2_idx)
{
    vec3 p1 = vertices[v1_idx];
    vec3 p2 = vertices[v2_idx];

    // 1. Identify 1-ring neighbors of the edge (v1, v2)
    unordered_set<int> ring_indices;
    if (v1_idx < vertex_neighbors.size())
    { // Check bounds
        for (int neighbor : vertex_neighbors[v1_idx])
        {
            if (neighbor != v1_idx && neighbor != v2_idx)
            {
                ring_indices.insert(neighbor);
            }
        }
    }
    if (v2_idx < vertex_neighbors.size())
    { // Check bounds
        for (int neighbor : vertex_neighbors[v2_idx])
        {
            if (neighbor != v1_idx && neighbor != v2_idx)
            {
                ring_indices.insert(neighbor);
            }
        }
    }

    if (ring_indices.empty())
    {                                   // Edge is isolated or connects two components with no other connections
        return {(p1 + p2) * 0.5f, 0.0}; // Or a high cost if such collapses are undesirable
    }

    // 2. Iteratively compute optimal position v_opt
    vec3 v_opt_current = (p1 + p2) * 0.5f;
    vector<double> current_weights(vertices.size(), 0.0); // Store weights for final cost calculation

    for (int iter = 0; iter < V_OPT_ITERATIONS; ++iter)
    {
        vec3 sum_weighted_positions(0.0f);
        double sum_weights = 0.0;

        bool possible_to_calc_weights = true;

        for (int P_k_idx : ring_indices)
        {
            if (P_k_idx >= vertices.size())
                continue; // Should not happen if graph is consistent
            vec3 P_k_pos = vertices[P_k_idx];

            // Find triangle partners for the new edge (v_opt_current, P_k_pos)
            std::vector<int> partners_idx = find_triangle_partners_for_new_edge(P_k_idx, v1_idx, v2_idx, ring_indices);

            double cotan_sum = 0.0;
            if (partners_idx.size() >= 1)
            {
                if (partners_idx[0] >= vertices.size())
                    continue;
                vec3 P_partner1_pos = vertices[partners_idx[0]];
                // Angle at P_partner1 in triangle (v_opt_current, P_k_pos, P_partner1_pos)
                cotan_sum += cotangent_robust(v_opt_current - P_partner1_pos, P_k_pos - P_partner1_pos);
            }
            if (partners_idx.size() >= 2)
            {
                if (partners_idx[1] >= vertices.size())
                    continue;
                vec3 P_partner2_pos = vertices[partners_idx[1]];
                // Angle at P_partner2 in triangle (v_opt_current, P_k_pos, P_partner2_pos)
                cotan_sum += cotangent_robust(v_opt_current - P_partner2_pos, P_k_pos - P_partner2_pos);
            }

            double weight_k = 0.0;
            if (!partners_idx.empty())
            {                               // If boundary edge, might only have one partner
                weight_k = cotan_sum / 2.0; // Common factor
            }

            // Clamp weights to be positive (common practice)
            weight_k = std::max((double)COT_EPSILON, weight_k);
            current_weights[P_k_idx] = weight_k;

            sum_weighted_positions += (float)weight_k * P_k_pos;
            sum_weights += weight_k;
        }

        if (sum_weights < COT_EPSILON)
        { // Not enough constraint or problematic weights
            // Fallback: v_opt_current remains midpoint or previous value if iter > 0
            // No update possible, break or use midpoint as final
            if (iter == 0)
                v_opt_current = (p1 + p2) * 0.5f; // Ensure it's at least midpoint
            break;
        }

        vec3 v_opt_next = sum_weighted_positions / (float)sum_weights;

        if (distance(v_opt_current, v_opt_next) < V_OPT_CONVERGENCE_THRESHOLD && iter > 0)
        {
            v_opt_current = v_opt_next;
            break;
        }
        v_opt_current = v_opt_next;
    }

    // 3. Calculate final cost using the determined v_opt_current and final weights
    double cost = 0.0;
    for (int P_k_idx : ring_indices)
    {
        if (P_k_idx >= vertices.size())
            continue;
        vec3 P_k_pos = vertices[P_k_idx];
        double weight_k = current_weights[P_k_idx];
        if (weight_k == 0.0)
        {                                             // if weight was not set due to no partners etc., re-calculate or use default
                                                      // This indicates an issue or an edge case not fully handled in weight calculation for v_opt iter.
                                                      // For simplicity, let's assume weights were computed. If not, this might be a small default.
                                                      // A more robust way would be to ensure all `current_weights[P_k_idx]` are properly set during iterations.
                                                      // For now, let's use a minimal contribution if weight is zero.
            cost += distance(v_opt_current, P_k_pos); // Or a factor of it
        }
        else
        {
            cost += weight_k * distance(v_opt_current, P_k_pos);
        }
    }

    // Normalize cost? Or scale? For now, raw energy.
    // Cost can be very small if v_opt perfectly centers. Can be large if it causes high tension.
    // A small penalty for collapsing to avoid issues with perfectly flat regions if cost is 0
    cost = std::max(cost, (double)COT_EPSILON * distance(p1, p2));

    return {v_opt_current, cost};
}

// Build adjacency list
void build_adjacency()
{
    int nV = vertices.size();
    vertex_neighbors.assign(nV, {});
    for (const auto &f : faces)
    {
        if (f.x < nV && f.y < nV && f.z < nV && f.x >= 0 && f.y >= 0 && f.z >= 0)
        {
            vertex_neighbors[f.x].insert(f.y);
            vertex_neighbors[f.x].insert(f.z);
            vertex_neighbors[f.y].insert(f.x);
            vertex_neighbors[f.y].insert(f.z);
            vertex_neighbors[f.z].insert(f.x);
            vertex_neighbors[f.z].insert(f.y);
        }
        else
        {
            // cerr << "Warning: Face with out-of-bounds vertex index during adjacency build." << endl;
        }
    }
}

// Initialize edge priority queue
void init_edges()
{
    edge_map.clear();
    while (!edge_queue.empty())
        edge_queue.pop();

    // Temporary adjacency build if not up-to-date (e.g., first call or after major changes)
    // build_adjacency(); // Ensure vertex_neighbors is correct. Called by simplify/simplifyStep entry.

    for (int v = 0; v < vertices.size(); ++v)
    {
        if (v >= vertex_neighbors.size())
            continue; // Should not happen
        for (int nb : vertex_neighbors[v])
        {
            if (nb >= vertices.size())
                continue; // Should not happen
            if (v < nb)
            { // Process each edge once
                pair<vec3, double> edge_data = compute_edge_cost(v, nb);
                Edge e{v, nb, edge_data.second, edge_data.first};
                int64_t k = edge_key(v, nb);
                edge_map[k] = e;
                edge_queue.push(e);
            }
        }
    }
}

// Update edges related to a specific vertex v (after a contraction)
void update_edges(int v)
{
    if (v >= vertices.size() || v >= vertex_neighbors.size())
        return; // Safety check

    // First, remove old edges connected to v from the map (they will be re-added or are gone)
    // This step is tricky because the map stores edges, and we need to iterate neighbors
    // to find keys. Iterating `edge_map` and removing is safer if keys are complex.
    // However, the current Q uses edges from `vertex_neighbors`, so recomputing for those is standard.

    for (int nb : vertex_neighbors[v])
    {
        if (nb >= vertices.size())
            continue;

        // remove old edge (v, nb) if it was in map from previous state.
        // The priority queue might still have old versions; they'll be skipped.
        // The map update is key.
        int64_t old_k = edge_key(v, nb);
        // edge_map.erase(old_k); // Erase might be too broad if some edges remain valid conceptually
        // but cost changes. It's safer to overwrite or add.

        pair<vec3, double> edge_data = compute_edge_cost(v, nb);
        Edge e{v, nb, edge_data.second, edge_data.first};

        int64_t k = edge_key(v, nb); // Key might involve old nb index if nb itself was merged into v.
                                     // The edge is always between v and a current neighbor nb.
        edge_map[k] = e;             // Update or insert
        edge_queue.push(e);          // Add new one; old ones in PQ will be skipped due to cost mismatch or map check.
    }
}

// Contract edge e = (v1, v2) by moving v1 to e.optimal_pos and removing v2
bool contract_edge(const Edge &e)
{
    int v1 = e.v1; // Vertex to keep
    int v2 = e.v2; // Vertex to remove

    if (v1 >= vertices.size() || v2 >= vertices.size() ||
        v1 >= vertex_neighbors.size() || v2 >= vertex_neighbors.size())
    {
        return false; // Invalid indices
    }

    // Check if edge still exists (v2 is a neighbor of v1)
    if (!vertex_neighbors[v1].count(v2))
    {
        return false; // Edge was already removed by a previous contraction
    }

    vertices[v1] = e.optimal_pos;
    // Quadric accumulation removed: quadrics[v1] += quadrics[v2];

    // Update topology:
    // For each neighbor 'nb' of v2:
    //  - Remove v2 from nb's neighbors list.
    //  - Add v1 to nb's neighbors list (if nb != v1).
    //  - Add nb to v1's neighbors list (if nb != v1).
    // Clear v2's neighbors list.
    // Remove v2 from v1's neighbors list.

    unordered_set<int> v2_neighbors_copy = vertex_neighbors[v2]; // Iterate on a copy

    for (int nb_of_v2 : v2_neighbors_copy)
    {
        if (nb_of_v2 == v1)
            continue; // Skip v1 itself

        if (nb_of_v2 < vertex_neighbors.size())
        {
            vertex_neighbors[nb_of_v2].erase(v2);
            vertex_neighbors[nb_of_v2].insert(v1);
        }
        vertex_neighbors[v1].insert(nb_of_v2); // Add nb_of_v2 to v1's neighbors
    }

    vertex_neighbors[v1].erase(v2); // v1 is no longer a neighbor of itself implicitly via v2
    if (!vertex_neighbors[v2].empty())
    {                                 // Check if not already cleared
        vertex_neighbors[v2].clear(); // v2 is isolated
    }

    // Update faces: replace all occurrences of v2 with v1
    // Remove degenerate faces (e.g., where two vertices become the same)
    vector<ivec3> new_faces;
    new_faces.reserve(faces.size());
    for (const auto &f : faces)
    {
        ivec3 new_f = f;
        if (new_f.x == v2)
            new_f.x = v1;
        if (new_f.y == v2)
            new_f.y = v1;
        if (new_f.z == v2)
            new_f.z = v1;

        // Check for degenerate faces
        if (new_f.x != new_f.y && new_f.y != new_f.z && new_f.x != new_f.z)
        {
            new_faces.push_back(new_f);
        }
    }
    faces.swap(new_faces);

    // After modifying topology around v1, update edges incident to v1
    update_edges(v1);

    // Update edges for neighbors of v2 that are now connected to v1
    // These neighbors' edges costs might change due to v1's new position and connectivity
    for (int nb_of_v2 : v2_neighbors_copy)
    {
        if (nb_of_v2 != v1)
        {                           // if nb_of_v2 became a new neighbor of v1
            update_edges(nb_of_v2); // Recompute costs for edges incident to these neighbors
        }
    }

    return true;
}

// Cleanup unused vertices (those not part of any face)
void cleanup()
{
    if (vertices.empty())
        return;
    int nV_old = vertices.size();
    vector<bool> used(nV_old, false);
    int current_max_idx = -1;
    for (const auto &f : faces)
    {
        used[f.x] = true;
        used[f.y] = true;
        used[f.z] = true;
        current_max_idx = std::max({current_max_idx, f.x, f.y, f.z});
    }

    if (current_max_idx == -1 && !faces.empty())
    { // Should not happen if faces exist
        // cerr << "Warning: Max index is -1 but faces exist. Problem in face data." << endl;
        return;
    }
    if (current_max_idx == -1 && faces.empty() && !vertices.empty())
    { // No faces, clear all vertices
        vertices.clear();
        vertex_neighbors.clear(); // Also clear adjacency
        return;
    }

    vector<vec3> new_vertices_list;
    new_vertices_list.reserve(current_max_idx + 1);
    vector<int> remap(nV_old, -1);
    int new_idx_counter = 0;

    for (int i = 0; i <= current_max_idx; ++i)
    { // Iterate only up to max_idx seen in faces
        if (i < nV_old && used[i])
        { // Ensure 'i' is a valid old index
            remap[i] = new_idx_counter;
            new_vertices_list.push_back(vertices[i]);
            new_idx_counter++;
        }
    }

    for (auto &f : faces)
    {
        f.x = remap[f.x];
        f.y = remap[f.y];
        f.z = remap[f.z];
        if (f.x == -1 || f.y == -1 || f.z == -1)
        {
            // This would indicate an issue, face references a removed vertex not properly handled.
            // cerr << "Error during cleanup: Face references remapped -1 index." << endl;
        }
    }
    vertices.swap(new_vertices_list);

    // Rebuild adjacency for the cleaned mesh
    // build_adjacency(); // Important if simplify/simplifyStep doesn't call it first
}

// Main simplification function
void simplify(float targetRatio)
{
    // No init_quadrics() needed
    build_adjacency(); // Build initial adjacency list
    init_edges();      // Populate priority queue with initial edge costs

    int initial_faces = faces.size();
    if (initial_faces == 0)
        return;
    int target_num_faces = int(initial_faces * (1.0f - targetRatio));
    if (target_num_faces < 0)
        target_num_faces = 0;

    int contractions = 0;
    const int max_contractions = initial_faces; // Heuristic limit

    while (faces.size() > target_num_faces && !edge_queue.empty() && contractions < max_contractions)
    {
        Edge e = edge_queue.top();
        edge_queue.pop();

        int64_t k = edge_key(e.v1, e.v2);

        // Check if edge is still valid and its cost is up-to-date in the map
        auto it = edge_map.find(k);
        if (it == edge_map.end() || fabs(it->second.cost - e.cost) > 1e-9)
        {
            continue; // Stale edge from priority queue
        }
        // Additional check: ensure vertices still exist (not cleaned up unexpectedly)
        if (e.v1 >= vertices.size() || e.v2 >= vertices.size() ||
            e.v1 >= vertex_neighbors.size() || e.v2 >= vertex_neighbors.size() ||
            !vertex_neighbors[e.v1].count(e.v2))
        {                       // Ensure they are still neighbors
            edge_map.erase(it); // Remove from map as it's invalid
            continue;
        }

        if (contract_edge(e))
        {
            edge_map.erase(k); // Edge successfully contracted, remove from map
            contractions++;
        }
        else
        {
            // Contraction failed (e.g. edge became invalid between pop and contract)
            // Remove from map to prevent further attempts on this stale version
            edge_map.erase(it);
        }
    }
    cleanup();         // Remove unused vertices and remap indices
    build_adjacency(); // Rebuild adjacency for the final cleaned mesh
}

// Single step simplification
void simplifyStep()
{
    // No init_quadrics()
    build_adjacency(); // Ensure adjacency is up-to-date before computing costs
    init_edges();      // Recalculate all edge costs and repopulate queue

    if (!edge_queue.empty())
    {
        Edge e = edge_queue.top();
        // Popping here means if contract_edge fails, this edge is lost for this step.
        // This matches QEM's behavior. For interactive use, might peek instead.
        edge_queue.pop();

        int64_t k = edge_key(e.v1, e.v2);

        auto it = edge_map.find(k);
        if (it == edge_map.end() || fabs(it->second.cost - e.cost) > 1e-9)
        {
            // Stale edge, typically we'd loop in simplify() to find the next best.
            // For simplifyStep(), if the top is stale, we might not simplify this step.
            return;
        }
        if (e.v1 >= vertices.size() || e.v2 >= vertices.size() ||
            e.v1 >= vertex_neighbors.size() || e.v2 >= vertex_neighbors.size() ||
            !vertex_neighbors[e.v1].count(e.v2))
        {
            edge_map.erase(it);
            return;
        }

        if (contract_edge(e))
        {
            edge_map.erase(k);
            best_edge = e; // Store the contracted edge
        }
        else
        {
            edge_map.erase(it);
        }
    }
    cleanup();
    build_adjacency(); // Rebuild for next step
}

// --- C Interface ---
extern "C"
{
    int SetMesh(vec3 *inV, ivec3 *inF, int nV, int nF)
    {
        vertices.assign(inV, inV + nV);
        faces.assign(inF, inF + nF);
        // Clear any old state
        edge_map.clear();
        while (!edge_queue.empty())
            edge_queue.pop();
        vertex_neighbors.clear();
        best_edge = Edge{-1, -1, 0.0, vec3(0.0f)};
        return 0;
    }

    // Important: These Get* functions return pointers to static buffers.
    // This is not thread-safe and data is overwritten on subsequent calls.
    // Common for simple C interfaces but be aware.
    vec3 *GetMeshVertices(int *outV)
    {
        *outV = vertices.size();
        static vector<vec3> buf_v; // static buffer
        buf_v = vertices;
        return buf_v.data();
    }

    ivec3 *GetMeshFaces(int *outF)
    {
        *outF = faces.size();
        static vector<ivec3> buf_f; // static buffer
        buf_f = faces;
        return buf_f.data();
    }

    void Simplify(float targetRatio)
    {
        simplify(targetRatio);
    }

    void SimplifyStep()
    {
        simplifyStep();
    }

    // Returns pointer to a static buffer containing {v1, v2} of last best_edge.
    // v1 is the vertex that was kept and moved, v2 was removed.
    int *GetActiveVertices()
    {
        static vector<int> buf_active_v(2);
        buf_active_v[0] = best_edge.v1; // The vertex that 'absorbed' v2
        buf_active_v[1] = best_edge.v2; // The vertex that was removed
        return buf_active_v.data();
    }
}

#endif // MESH_SIMPLIFICATION_CPP