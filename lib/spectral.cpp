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

using namespace std;
using namespace glm;

const float COT_EPSILON = 1e-6f;
const int V_OPT_ITERATIONS = 5;
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
static vector<ivec3> faces;
static priority_queue<Edge, vector<Edge>, greater<Edge>> edge_queue;
static unordered_map<int64_t, Edge> edge_map;
static vector<unordered_set<int>> vertex_neighbors;    // vertex_idx -> set of neighbor vertex_indices
static vector<unordered_set<int>> vertex_face_indices; // vertex_idx -> set of incident face_indices (indices into `faces` vector)
static Edge best_edge;
static unordered_map<int, int> id_map; // <old_id, new_id>

// Global state flags
static bool s_mesh_structures_valid = false; // Are vertices, faces, vertex_neighbors, vertex_face_indices up-to-date?
static bool s_priority_queue_valid = false;  // Is edge_queue and edge_map reflecting current mesh state?

// Calculate unique edge key
inline int64_t edge_key(int a, int b)
{
    if (a > b)
        std::swap(a, b);
    return ((int64_t)a << 32) | (uint32_t)b;
}

// Robust cotangent of angle <apb (angle at p)
float cotangent_robust(const vec3 &v_pa, const vec3 &v_pb)
{
    float dot_product = dot(v_pa, v_pb);
    float cross_product_mag = length(cross(v_pa, v_pb));

    if (cross_product_mag < COT_EPSILON)
    {
        return (dot_product > 0) ? (1.0f / COT_EPSILON) : (-1.0f / COT_EPSILON);
    }
    return dot_product / cross_product_mag;
}

// Helper to get the other two vertices from a face given one vertex
// Returns {-1, -1} if v is not in the face.
std::pair<int, int> get_other_two_verts(const ivec3 &face, int v)
{
    if (face.x == v)
        return {face.y, face.z};
    if (face.y == v)
        return {face.x, face.z};
    if (face.z == v)
        return {face.x, face.y};
    return {-1, -1};
}

// Calculate edge contraction optimal position and cost
pair<vec3, double> compute_edge_cost(int v1_idx, int v2_idx)
{
    if (v1_idx >= vertices.size() || v2_idx >= vertices.size() ||
        v1_idx < 0 || v2_idx < 0)
    { // Basic sanity check
        return {vec3(0.0f), std::numeric_limits<double>::max()};
    }
    vec3 p1 = vertices[v1_idx];
    vec3 p2 = vertices[v2_idx];

    unordered_set<int> ring_indices;
    if (v1_idx < vertex_neighbors.size())
    {
        for (int neighbor : vertex_neighbors[v1_idx])
        {
            if (neighbor != v2_idx)
                ring_indices.insert(neighbor);
        }
    }
    if (v2_idx < vertex_neighbors.size())
    {
        for (int neighbor : vertex_neighbors[v2_idx])
        {
            if (neighbor != v1_idx)
                ring_indices.insert(neighbor);
        }
    }

    if (ring_indices.empty())
    {
        return {(p1 + p2) * 0.5f, 0.0}; // Or high cost
    }

    vec3 v_opt_current = (p1 + p2) * 0.5f;
    unordered_map<int, double> current_weights; // P_k_idx -> weight_k

    for (int iter = 0; iter < V_OPT_ITERATIONS; ++iter)
    {
        vec3 sum_weighted_positions(0.0f);
        double sum_weights = 0.0;
        current_weights.clear();

        for (int P_k_idx : ring_indices)
        {
            if (P_k_idx >= vertices.size() || P_k_idx < 0)
                continue;
            vec3 P_k_pos = vertices[P_k_idx];

            vector<int> partners_idx_list;
            unordered_set<int> distinct_partners_set;

            auto find_partners_for_pk = [&](int center_v_original_idx)
            {
                if (center_v_original_idx >= vertex_face_indices.size() || center_v_original_idx < 0)
                    return;
                for (int face_idx : vertex_face_indices[center_v_original_idx])
                {
                    if (face_idx >= faces.size() || face_idx < 0)
                        continue;
                    const ivec3 &f = faces[face_idx];
                    pair<int, int> others = get_other_two_verts(f, center_v_original_idx); // Other two verts in face f besides center_v

                    int potential_partner = -1;
                    if (others.first == P_k_idx && ring_indices.count(others.second))
                    {
                        potential_partner = others.second;
                    }
                    else if (others.second == P_k_idx && ring_indices.count(others.first))
                    {
                        potential_partner = others.first;
                    }

                    if (potential_partner != -1 && potential_partner != P_k_idx)
                    { // Partner must be in ring and not P_k itself
                        if (distinct_partners_set.find(potential_partner) == distinct_partners_set.end())
                        {
                            partners_idx_list.push_back(potential_partner);
                            distinct_partners_set.insert(potential_partner);
                        }
                    }
                }
            };

            find_partners_for_pk(v1_idx);
            find_partners_for_pk(v2_idx);

            double cotan_sum = 0.0;
            for (int partner_vertex_idx : partners_idx_list)
            {
                if (partner_vertex_idx >= vertices.size() || partner_vertex_idx < 0)
                    continue;
                vec3 P_partner_pos = vertices[partner_vertex_idx];
                // Angle at P_partner_pos in triangle (v_opt_current, P_k_pos, P_partner_pos)
                cotan_sum += cotangent_robust(v_opt_current - P_partner_pos, P_k_pos - P_partner_pos);
            }

            double weight_k = 0.0;
            // If P_k forms two triangles in the new 1-ring of v_opt, there are two cotan terms.
            // If P_k is on a boundary (of the 1-ring fan around v_opt), only one term.
            // The sum is used directly. Common factor of 0.5 often appears if formula is sum(0.5 * cot * ||...||^2)
            // Here we use w_k = sum_cotangents.
            if (!partners_idx_list.empty())
            {
                weight_k = cotan_sum; // Or cotan_sum / 2.0 depending on formulation. Let's try direct sum.
            }

            weight_k = std::max((double)COT_EPSILON, weight_k); // Clamp weights
            current_weights[P_k_idx] = weight_k;

            sum_weighted_positions += (float)weight_k * P_k_pos;
            sum_weights += weight_k;
        }

        if (sum_weights < COT_EPSILON)
        {
            if (iter == 0)
                v_opt_current = (p1 + p2) * 0.5f;
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

    double cost = 0.0;
    for (int P_k_idx : ring_indices)
    {
        if (P_k_idx >= vertices.size() || P_k_idx < 0)
            continue;
        vec3 P_k_pos = vertices[P_k_idx];
        auto it = current_weights.find(P_k_idx);
        double weight_k = (it != current_weights.end()) ? it->second : COT_EPSILON; // Use stored or default small
        cost += weight_k * distance(v_opt_current, P_k_pos);
    }
    cost = std::max(cost, (double)COT_EPSILON * distance(p1, p2));
    return {v_opt_current, cost};
}

// Build vertex-vertex and vertex-face adjacencies
void build_adjacency_and_face_indices()
{
    int nV = vertices.size();
    vertex_neighbors.assign(nV, {});
    vertex_face_indices.assign(nV, {});
    for (int i = 0; i < faces.size(); ++i)
    {
        const auto &f = faces[i];
        if (f.x < 0 || f.x >= nV || f.y < 0 || f.y >= nV || f.z < 0 || f.z >= nV)
            continue;

        vertex_neighbors[f.x].insert(f.y);
        vertex_neighbors[f.x].insert(f.z);
        vertex_neighbors[f.y].insert(f.x);
        vertex_neighbors[f.y].insert(f.z);
        vertex_neighbors[f.z].insert(f.x);
        vertex_neighbors[f.z].insert(f.y);

        vertex_face_indices[f.x].insert(i);
        vertex_face_indices[f.y].insert(i);
        vertex_face_indices[f.z].insert(i);
    }
    s_mesh_structures_valid = true;
}

// Initialize edge priority queue
void init_edges()
{
    if (!s_mesh_structures_valid)
    { // Should be called after build_adjacency_and_face_indices
        // cerr << "Warning: init_edges called with invalid mesh structures." << endl;
        build_adjacency_and_face_indices(); // Attempt to recover
    }

    edge_map.clear();
    while (!edge_queue.empty())
        edge_queue.pop();

    for (int v = 0; v < vertices.size(); ++v)
    {
        if (v >= vertex_neighbors.size())
            continue;
        for (int nb : vertex_neighbors[v])
        {
            if (nb >= vertices.size())
                continue;
            if (v < nb)
            {
                pair<vec3, double> edge_data = compute_edge_cost(v, nb);
                Edge e{v, nb, edge_data.second, edge_data.first};
                edge_map[edge_key(v, nb)] = e;
                edge_queue.push(e);
            }
        }
    }
    s_priority_queue_valid = true;
}

// Update edges related to vertex v and its neighbors
void update_edges_for_vertex_and_its_neighbors(int v)
{
    if (v >= vertices.size() || v < 0 || v >= vertex_neighbors.size())
        return;

    unordered_set<int> vertices_to_update;
    vertices_to_update.insert(v);
    for (int nb_of_v : vertex_neighbors[v])
    {
        vertices_to_update.insert(nb_of_v);
    }

    for (int current_v : vertices_to_update)
    {
        if (current_v >= vertices.size() || current_v < 0 || current_v >= vertex_neighbors.size())
            continue;
        for (int nb : vertex_neighbors[current_v])
        {
            if (nb >= vertices.size() || nb < 0)
                continue;
            // Only recompute for edge (current_v, nb) once, and current_v must be less than nb to avoid double counting
            // or re-evaluate all edges from current_v
            // if (current_v < nb) { // This might miss edges if nb was also affected, simpler to update all from current_v
            pair<vec3, double> edge_data = compute_edge_cost(current_v, nb);
            Edge e{current_v, nb, edge_data.second, edge_data.first};
            edge_map[edge_key(current_v, nb)] = e; // Update or insert
            edge_queue.push(e);
            // }
        }
    }
}

// Contract edge e = (v1, v2)
bool contract_edge(const Edge &e)
{
    int v_keep = e.v1;   // Vertex to keep
    int v_remove = e.v2; // Vertex to remove

    if (v_keep >= vertices.size() || v_remove >= vertices.size() || v_keep < 0 || v_remove < 0 ||
        v_keep >= vertex_neighbors.size() || v_remove >= vertex_neighbors.size() || // Adjacency might not be full size if sparse
        !vertex_neighbors[v_keep].count(v_remove))
    { // Check if edge still exists
        return false;
    }

    unordered_set<int> old_v_remove_neighbors = vertex_neighbors[v_remove]; // Copy before modifying

    vertices[v_keep] = e.optimal_pos;

    // Rebuild faces list:
    vector<ivec3> new_faces_list;
    new_faces_list.reserve(faces.size());
    for (const auto &f : faces)
    {
        ivec3 new_f = f;
        bool changed = false;
        if (new_f.x == v_remove)
        {
            new_f.x = v_keep;
            changed = true;
        }
        if (new_f.y == v_remove)
        {
            new_f.y = v_keep;
            changed = true;
        }
        if (new_f.z == v_remove)
        {
            new_f.z = v_keep;
            changed = true;
        }

        if (new_f.x != new_f.y && new_f.y != new_f.z && new_f.x != new_f.z)
        {
            new_faces_list.push_back(new_f);
        }
    }
    faces.swap(new_faces_list);

    // Adjacency structures are now stale. Rebuild them.
    build_adjacency_and_face_indices(); // This is O(N_faces_new)
    s_mesh_structures_valid = true;     // Rebuilt

    // Costs for edges around v_keep and former neighbors of v_remove are stale.
    // update_edges_for_vertex_and_its_neighbors will recompute and push to PQ.
    update_edges_for_vertex_and_its_neighbors(v_keep);
    for (int old_nb_of_v_remove : old_v_remove_neighbors)
    {
        if (old_nb_of_v_remove != v_keep && old_nb_of_v_remove < vertices.size())
        { // vertices.size check, as it might have been removed earlier by another op if cleanup not aggressive
            update_edges_for_vertex_and_its_neighbors(old_nb_of_v_remove);
        }
    }
    s_priority_queue_valid = true; // PQ has new items, stale ones will be ignored by map check

    return true;
}

// Cleanup unused vertices
void cleanup()
{
    id_map.clear();
    if (vertices.empty())
    {
        s_mesh_structures_valid = true; // Empty but valid
        s_priority_queue_valid = false; // No edges
        return;
    }

    int nV_old = vertices.size();
    vector<bool> used(nV_old, false);
    int max_idx_in_faces = 0;
    for (const auto &f : faces)
    {
        if (f.x < 0 || f.y < 0 || f.z < 0)
            continue; // Should not happen with valid faces
        max_idx_in_faces = std::max({max_idx_in_faces, f.x, f.y, f.z});
        if (f.x < nV_old)
            used[f.x] = true;
        if (f.y < nV_old)
            used[f.y] = true;
        if (f.z < nV_old)
            used[f.z] = true;
    }

    vector<vec3> new_vertices_list;
    vector<int> remap(max_idx_in_faces + 1, -1); // Ensure remap is large enough
                                                 // Or remap based on nV_old if safer
    if (nV_old > max_idx_in_faces + 1)
        remap.resize(nV_old, -1);

    int current_new_idx = 0;
    for (int i = 0; i < nV_old; ++i)
    { // Iterate old vertex indices
        if (used[i])
        {
            if (i < remap.size())
            {
                remap[i] = current_new_idx;
                id_map[i] = current_new_idx;
            }
            else
            { /* Error: remap array too small, vertex index out of expected bounds */
            }
            new_vertices_list.push_back(vertices[i]);
            current_new_idx++;
        }
    }

    for (auto &f : faces)
    {
        if (f.x < remap.size())
            f.x = remap[f.x];
        else
            f.x = -1; // Mark invalid if out of bound
        if (f.y < remap.size())
            f.y = remap[f.y];
        else
            f.y = -1;
        if (f.z < remap.size())
            f.z = remap[f.z];
        else
            f.z = -1;
        if (f.x == -1 || f.y == -1 || f.z == -1)
        {
            // This face is now invalid, should be removed.
            // For simplicity, let's assume this leads to degenerate check later or filter here.
        }
    }
    // Filter out faces that became invalid due to remapping issues
    faces.erase(std::remove_if(faces.begin(), faces.end(), [](const ivec3 &f)
                               { return f.x == -1 || f.y == -1 || f.z == -1 || f.x == f.y || f.y == f.z || f.x == f.z; }),
                faces.end());

    vertices.swap(new_vertices_list);

    // After cleanup, all adjacencies and PQ are invalid because indices changed.
    s_mesh_structures_valid = false; // Needs rebuild_adjacency_and_face_indices
    s_priority_queue_valid = false;  // Needs init_edges
}

// Main simplification function
void simplify(float targetRatio)
{
    if (!s_mesh_structures_valid)
    { // If SetMesh wasn't called or structures invalidated
        build_adjacency_and_face_indices();
    }
    if (!s_priority_queue_valid)
    { // If PQ is not reflecting current mesh (e.g. after cleanup)
        init_edges();
    }

    int initial_faces = faces.size();
    if (initial_faces == 0)
        return;
    int target_num_faces = int(initial_faces * (1.0f - targetRatio));
    if (target_num_faces < 0)
        target_num_faces = 0;

    int contractions = 0;
    const int max_contractions = initial_faces;

    while (faces.size() > target_num_faces && !edge_queue.empty() && contractions < max_contractions)
    {
        Edge e = edge_queue.top();
        edge_queue.pop();

        int64_t k = edge_key(e.v1, e.v2);
        auto it = edge_map.find(k);

        if (it == edge_map.end() || fabs(it->second.cost - e.cost) > 1e-9)
        {
            continue; // Stale edge from PQ
        }
        // Additional check: ensure vertices and edge still valid in current (possibly modified) graph
        if (e.v1 >= vertices.size() || e.v2 >= vertices.size() ||
            e.v1 < 0 || e.v2 < 0 ||
            e.v1 >= vertex_neighbors.size() || e.v2 >= vertex_neighbors.size() || // vertex_neighbors might shrink
            !vertex_neighbors[e.v1].count(e.v2))
        {
            edge_map.erase(it);
            continue;
        }

        if (contract_edge(e))
        {
            edge_map.erase(k);
            contractions++;
        }
        else
        {
            edge_map.erase(it); // Contraction failed, remove from map
        }
    }
    cleanup(); // Final cleanup
    // After cleanup, structures need rebuild if GetMesh is called or another simplify op.
    // Let GetMesh handles ensure structures are fine, or user calls SetMesh again.
    // For now, `simplify` leaves `s_mesh_structures_valid` potentially false (due to cleanup).
}

// Single step simplification
void simplifyStep()
{
    if (!s_mesh_structures_valid)
    {
        build_adjacency_and_face_indices();
    }
    if (!s_priority_queue_valid)
    {
        init_edges();
    }

    if (edge_queue.empty())
        return;

    Edge e = {-1, -1, 0, vec3(0)};
    bool found_valid_edge = false;

    // Pop until a valid edge is found or queue is empty
    while (!edge_queue.empty())
    {
        e = edge_queue.top();
        edge_queue.pop();
        int64_t k = edge_key(e.v1, e.v2);
        auto it = edge_map.find(k);

        if (it == edge_map.end() || fabs(it->second.cost - e.cost) > 1e-9)
        {
            continue; // Stale
        }
        if (e.v1 >= vertices.size() || e.v2 >= vertices.size() ||
            e.v1 < 0 || e.v2 < 0 ||
            e.v1 >= vertex_neighbors.size() || e.v2 >= vertex_neighbors.size() ||
            !vertex_neighbors[e.v1].count(e.v2))
        {
            edge_map.erase(it);
            continue; // Invalid
        }
        found_valid_edge = true;
        best_edge = e; // Store original indices for GetActiveVertices
        if (contract_edge(e))
        {
            edge_map.erase(k);
        }
        else
        {
            edge_map.erase(it); // Contraction failed for some reason
        }
        break;
    }

    cleanup(); // Cleanup after the step
    // s_mesh_structures_valid and s_priority_queue_valid will be false now.
}

// --- C Interface ---
extern "C"
{
    int SetMesh(vec3 *inV, ivec3 *inF, int nV, int nF)
    {
        vertices.assign(inV, inV + nV);
        faces.assign(inF, inF + nF);

        edge_map.clear();
        while (!edge_queue.empty())
            edge_queue.pop();
        vertex_neighbors.clear();
        vertex_face_indices.clear();
        best_edge = Edge{-1, -1, 0.0, vec3(0.0f)};

        s_mesh_structures_valid = false; // Force rebuild on first operation
        s_priority_queue_valid = false;
        return 0;
    }

    vec3 *GetMeshVertices(int *outV)
    {
        if (!s_mesh_structures_valid)
        { // If cleanup was last op, adjacencies need rebuild for consistency
            build_adjacency_and_face_indices();
        }
        *outV = vertices.size();
        static vector<vec3> buf_v;
        buf_v = vertices;
        return buf_v.data();
    }

    ivec3 *GetMeshFaces(int *outF)
    {
        if (!s_mesh_structures_valid)
        {
            build_adjacency_and_face_indices();
        }
        *outF = faces.size();
        static vector<ivec3> buf_f;
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

    int *GetActiveVertices()
    {
        static vector<int> buf;
        buf.clear();
        buf.push_back(best_edge.v1);
        buf.push_back(best_edge.v2);
        buf.push_back(id_map[best_edge.v1]);
        return buf.data();
    }
}

#endif // MESH_SIMPLIFICATION_CPP