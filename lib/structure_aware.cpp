#ifndef STRUCTURE_AWARE_MESH_SIMPLIFICATION_CPP
#define STRUCTURE_AWARE_MESH_SIMPLIFICATION_CPP

#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>

using namespace std;
using namespace glm;

const float U_PARAMETER = 0.5f;
const int MAX_ITERATIONS = 600000;
const int SHAPE_PRESERVE_MAX_ITERATIONS = 100000;

// Structure-aware edge definition
struct StructureEdge
{
    int v1, v2;
    double cost;
    vec3 optimal_pos;
    vector<int> v1_proxies;
    vector<int> v2_proxies;
    bool is_boundary;
    bool operator>(const StructureEdge &o) const { return cost > o.cost; }
};

// Planar proxy structure
struct PlanarProxy
{
    vec4 plane_eq; // plane equation (a, b, c, d)
    mat4 quadric;
    vector<int> face_indices;
    int proxy_id;
};

// Global data
static vector<vec3> vertices;
static vector<ivec3> faces;
static vector<mat4> quadrics;
static vector<vector<int>> vertex_proxies;
static vector<PlanarProxy> proxies;
static unordered_map<int, mat4> proxy_equations;
static unordered_map<int, int> proxy_face_count;
static priority_queue<StructureEdge, vector<StructureEdge>, greater<StructureEdge>> edge_queue;
static unordered_map<int64_t, StructureEdge> edge_map;
static vector<unordered_set<int>> vertex_neighbors;
static vector<vector<int>> face_proxies;
static StructureEdge best_edge;
static unordered_map<int, int> id_map; // <old_id, new_id>

// Flags for different decimation modes
static bool shape_preserve_flag = false;
static bool shape_preserve_decimation_flag = false;
static bool modify_cost = false;

// Random number generator
static mt19937 rng(random_device{}());

inline int64_t edge_key(int a, int b)
{
    if (a > b)
        std::swap(a, b);
    return ((int64_t)a << 32) | (uint32_t)b;
}

// Compute face area
float compute_face_area(int face_idx)
{
    const ivec3 &f = faces[face_idx];
    vec3 v0 = vertices[f.x];
    vec3 v1 = vertices[f.y];
    vec3 v2 = vertices[f.z];

    vec3 cross_product = cross(v2 - v0, v1 - v0);
    return length(cross_product) * 0.5f;
}

// Compute face quadric
mat4 compute_face_quadric(const vec3 &v0, const vec3 &v1, const vec3 &v2)
{
    vec3 n = normalize(cross(v1 - v0, v2 - v0));
    float d = -dot(n, v0);
    vec4 plane(n, d);

    mat4 Q(0.0f);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            Q[i][j] = plane[i] * plane[j];
    return Q;
}

// Initialize planar proxies using region growing
void initialize_proxies()
{
    int num_faces = faces.size();
    vector<bool> assigned(num_faces, false);
    vector<int> face_to_proxy(num_faces, -1);
    int proxy_id = 0;

    face_proxies.assign(num_faces, vector<int>());
    proxies.clear();
    proxy_face_count.clear();

    // Simple proxy assignment - each face gets its own proxy initially
    // In practice, you'd use region growing based on normal similarity
    for (int i = 0; i < num_faces; ++i)
    {
        if (!assigned[i])
        {
            PlanarProxy proxy;
            proxy.proxy_id = proxy_id;

            const ivec3 &f = faces[i];
            vec3 v0 = vertices[f.x];
            vec3 v1 = vertices[f.y];
            vec3 v2 = vertices[f.z];

            vec3 n = normalize(cross(v1 - v0, v2 - v0));
            float d = -dot(n, v0);
            proxy.plane_eq = vec4(n, d);
            proxy.quadric = compute_face_quadric(v0, v1, v2);
            proxy.face_indices.push_back(i);

            face_to_proxy[i] = proxy_id;
            face_proxies[i].push_back(proxy_id);
            assigned[i] = true;

            proxies.push_back(proxy);
            proxy_face_count[proxy_id] = 1;
            proxy_id++;
        }
    }

    // Build vertex to proxy mapping
    vertex_proxies.assign(vertices.size(), vector<int>());
    for (int i = 0; i < num_faces; ++i)
    {
        const ivec3 &f = faces[i];
        int proxy = face_to_proxy[i];

        auto &v0_proxies = vertex_proxies[f.x];
        auto &v1_proxies = vertex_proxies[f.y];
        auto &v2_proxies = vertex_proxies[f.z];

        if (find(v0_proxies.begin(), v0_proxies.end(), proxy) == v0_proxies.end())
            v0_proxies.push_back(proxy);
        if (find(v1_proxies.begin(), v1_proxies.end(), proxy) == v1_proxies.end())
            v1_proxies.push_back(proxy);
        if (find(v2_proxies.begin(), v2_proxies.end(), proxy) == v2_proxies.end())
            v2_proxies.push_back(proxy);
    }
}

// Calculate planar proxy equations
void calculate_planar_proxy_equations()
{
    proxy_equations.clear();

    for (const auto &proxy : proxies)
    {
        mat4 total_quadric(0.0f);

        for (int face_idx : proxy.face_indices)
        {
            const ivec3 &f = faces[face_idx];
            vec3 v0 = vertices[f.x];
            vec3 v1 = vertices[f.y];
            vec3 v2 = vertices[f.z];

            mat4 face_quadric = compute_face_quadric(v0, v1, v2);
            total_quadric += face_quadric;
        }

        proxy_equations[proxy.proxy_id] = total_quadric;
    }
}

// Initialize vertex quadrics with structure awareness
void init_structure_quadrics()
{
    int nV = vertices.size();
    quadrics.assign(nV, mat4(0.0f));

    for (int i = 0; i < faces.size(); ++i)
    {
        const ivec3 &f = faces[i];
        mat4 Qf = compute_face_quadric(vertices[f.x], vertices[f.y], vertices[f.z]);

        if (modify_cost && !face_proxies[i].empty())
        {
            int proxy_id = face_proxies[i][0];
            if (proxy_equations.count(proxy_id))
            {
                mat4 proxy_quad = proxy_equations[proxy_id];
                Qf = Qf * (1.0f - U_PARAMETER) + proxy_quad * U_PARAMETER;
            }
        }

        quadrics[f.x] += Qf;
        quadrics[f.y] += Qf;
        quadrics[f.z] += Qf;
    }
}

// Build adjacency
void build_adjacency()
{
    int nV = vertices.size();
    vertex_neighbors.assign(nV, unordered_set<int>());

    for (const auto &f : faces)
    {
        vertex_neighbors[f.x].insert(f.y);
        vertex_neighbors[f.x].insert(f.z);
        vertex_neighbors[f.y].insert(f.x);
        vertex_neighbors[f.y].insert(f.z);
        vertex_neighbors[f.z].insert(f.x);
        vertex_neighbors[f.z].insert(f.y);
    }
}

// Compute boundary quadric for structure preservation
mat4 compute_boundary_quadric(int v1, int v2)
{
    mat4 boundary_quadric(0.0f);

    // Find faces adjacent to this edge
    vector<int> adjacent_faces;
    for (int i = 0; i < faces.size(); ++i)
    {
        const ivec3 &f = faces[i];
        if ((f.x == v1 && f.y == v2) || (f.y == v1 && f.z == v2) || (f.z == v1 && f.x == v2) ||
            (f.x == v2 && f.y == v1) || (f.y == v2 && f.z == v1) || (f.z == v2 && f.x == v1))
        {
            adjacent_faces.push_back(i);
        }
    }

    if (adjacent_faces.size() == 2)
    {
        // Check if faces belong to different proxies
        int proxy1 = face_proxies[adjacent_faces[0]].empty() ? -1 : face_proxies[adjacent_faces[0]][0];
        int proxy2 = face_proxies[adjacent_faces[1]].empty() ? -1 : face_proxies[adjacent_faces[1]][0];

        if (proxy1 != proxy2 && proxy1 != -1 && proxy2 != -1)
        {
            float area1 = compute_face_area(adjacent_faces[0]);
            float area2 = compute_face_area(adjacent_faces[1]);

            // Create orthogonal plane quadrics (simplified version)
            vec3 edge_vec = normalize(vertices[v2] - vertices[v1]);
            vec3 n1 = normalize(cross(vertices[faces[adjacent_faces[0]].y] - vertices[faces[adjacent_faces[0]].x],
                                      vertices[faces[adjacent_faces[0]].z] - vertices[faces[adjacent_faces[0]].x]));
            vec3 ortho_normal = normalize(cross(n1, edge_vec));
            float d = -dot(ortho_normal, vertices[v1]);

            vec4 plane(ortho_normal, d);
            mat4 ortho_quad(0.0f);
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    ortho_quad[i][j] = plane[i] * plane[j];

            boundary_quadric = ortho_quad * (area1 + area2);
        }
    }

    return boundary_quadric;
}

// Compute optimal position for structure-aware edge collapse
pair<vec3, double> compute_structure_edge_cost(int v1, int v2)
{
    // Check if vertices have different proxy sets (structure preservation)
    if (shape_preserve_decimation_flag)
    {
        const auto &proxies1 = vertex_proxies[v1];
        const auto &proxies2 = vertex_proxies[v2];

        if (proxies1.size() != proxies2.size())
        {
            return {0.5f * (vertices[v1] + vertices[v2]), 1e6}; // High cost
        }

        // Check if all proxies match
        for (int proxy : proxies1)
        {
            if (find(proxies2.begin(), proxies2.end(), proxy) == proxies2.end())
            {
                return {0.5f * (vertices[v1] + vertices[v2]), 1e6}; // High cost
            }
        }
    }

    mat4 Q = quadrics[v1] + quadrics[v2];

    // Add boundary quadric if preserving structure
    if (shape_preserve_flag)
    {
        Q += compute_boundary_quadric(v1, v2);
    }

    // Prefer vertices with more proxies (corners/features)
    vec3 optP;
    if (vertex_proxies[v1].size() >= 3 && vertex_proxies[v2].size() < 3)
    {
        optP = vertices[v1];
    }
    else if (vertex_proxies[v2].size() >= 3 && vertex_proxies[v1].size() < 3)
    {
        optP = vertices[v2];
    }
    else
    {
        // Standard QEM optimization
        mat4 Qp = Q;
        for (int col = 0; col < 4; ++col)
        {
            Qp[col][3] = (col == 3) ? 1.0f : 0.0f;
        }

        float det = determinant(Qp);
        if (fabs(det) > 1e-6f)
        {
            mat4 invQ = inverse(Qp);
            vec4 res = invQ * vec4(0, 0, 0, 1);
            optP = vec3(res);

            // Sanity check for reasonable position
            vec3 midpoint = 0.5f * (vertices[v1] + vertices[v2]);
            if (isnan(optP.x) || isnan(optP.y) || isnan(optP.z) ||
                length(optP - midpoint) > 10.0f * length(midpoint))
            {
                optP = midpoint;
            }
        }
        else
        {
            optP = 0.5f * (vertices[v1] + vertices[v2]);
        }
    }

    vec4 vH(optP, 1.0f);
    double cost = double(dot(vH, Q * vH));

    // Add structure preservation penalty
    if (shape_preserve_flag)
    {
        // Penalty for proxy mismatch
        const auto &proxies1 = vertex_proxies[v1];
        const auto &proxies2 = vertex_proxies[v2];

        set<int> set1(proxies1.begin(), proxies1.end());
        set<int> set2(proxies2.begin(), proxies2.end());
        set<int> diff1, diff2;

        set_difference(set1.begin(), set1.end(), set2.begin(), set2.end(),
                       inserter(diff1, diff1.begin()));
        set_difference(set2.begin(), set2.end(), set1.begin(), set1.end(),
                       inserter(diff2, diff2.begin()));

        if (!diff1.empty() && !diff2.empty())
        {
            int mismatch = std::min(diff1.size(), diff2.size());
            cost += 10.0 * mismatch; // Penalty for proxy mismatch
        }

        // Area-based penalty for boundary edges
        vector<int> adjacent_faces;
        for (int i = 0; i < faces.size(); ++i)
        {
            const ivec3 &f = faces[i];
            if ((f.x == v1 && (f.y == v2 || f.z == v2)) ||
                (f.y == v1 && (f.x == v2 || f.z == v2)) ||
                (f.z == v1 && (f.x == v2 || f.y == v2)))
            {
                adjacent_faces.push_back(i);
            }
        }

        if (adjacent_faces.size() == 2)
        {
            int proxy1 = face_proxies[adjacent_faces[0]].empty() ? -1 : face_proxies[adjacent_faces[0]][0];
            int proxy2 = face_proxies[adjacent_faces[1]].empty() ? -1 : face_proxies[adjacent_faces[1]][0];

            if (proxy1 != proxy2 && proxy1 != -1 && proxy2 != -1)
            {
                float area1 = compute_face_area(adjacent_faces[0]);
                float area2 = compute_face_area(adjacent_faces[1]);
                cost += area1 + area2; // Area penalty for boundary edges
            }
        }
    }

    return {optP, cost};
}

// Validate edge for collapse with structure awareness
bool validate_structure_edge_collapse(int v1, int v2)
{
    // Standard topological validation
    set<int> v1_neighbors, v2_neighbors, common_neighbors;

    for (int nb : vertex_neighbors[v1])
    {
        if (nb != v2)
            v1_neighbors.insert(nb);
    }
    for (int nb : vertex_neighbors[v2])
    {
        if (nb != v1)
            v2_neighbors.insert(nb);
    }

    set_intersection(v1_neighbors.begin(), v1_neighbors.end(),
                     v2_neighbors.begin(), v2_neighbors.end(),
                     inserter(common_neighbors, common_neighbors.begin()));

    if (common_neighbors.size() != 2)
        return false; // Non-manifold result

    // // Structure preservation checks
    // if (shape_preserve_flag)
    // {
    //     // Check proxy face count constraints
    //     vector<int> adjacent_faces;
    //     for (int i = 0; i < faces.size(); ++i)
    //     {
    //         const ivec3 &f = faces[i];
    //         if ((f.x == v1 && (f.y == v2 || f.z == v2)) ||
    //             (f.y == v1 && (f.x == v2 || f.z == v2)) ||
    //             (f.z == v1 && (f.x == v2 || f.y == v2)))
    //         {
    //             adjacent_faces.push_back(i);
    //         }
    //     }

    //     if (adjacent_faces.size() == 2)
    //     {
    //         int proxy1 = face_proxies[adjacent_faces[0]].empty() ? -1 : face_proxies[adjacent_faces[0]][0];
    //         int proxy2 = face_proxies[adjacent_faces[1]].empty() ? -1 : face_proxies[adjacent_faces[1]][0];

    //         if (proxy1 != -1 && proxy_face_count[proxy1] <= 2)
    //             return false;
    //         if (proxy2 != -1 && proxy_face_count[proxy2] <= 2)
    //             return false;
    //     }
    // }

    return true;
}

// Initialize structure-aware edges
void init_structure_edges()
{
    edge_map.clear();
    while (!edge_queue.empty())
        edge_queue.pop();

    for (int v = 0; v < vertices.size(); ++v)
    {
        for (int nb : vertex_neighbors[v])
        {
            if (v < nb)
            {
                auto [optP, cost] = compute_structure_edge_cost(v, nb);
                StructureEdge e;
                e.v1 = v;
                e.v2 = nb;
                e.cost = cost;
                e.optimal_pos = optP;
                e.v1_proxies = vertex_proxies[v];
                e.v2_proxies = vertex_proxies[nb];
                e.is_boundary = false; // Could be computed based on adjacent faces

                int64_t k = edge_key(v, nb);
                edge_map[k] = e;
                edge_queue.push(e);
            }
        }
    }
}

// Update edges after collapse
void update_structure_edges(int v)
{
    for (int nb : vertex_neighbors[v])
    {
        int64_t k = edge_key(v, nb);
        auto [optP, cost] = compute_structure_edge_cost(v, nb);

        StructureEdge e;
        e.v1 = v;
        e.v2 = nb;
        e.cost = cost;
        e.optimal_pos = optP;
        e.v1_proxies = vertex_proxies[v];
        e.v2_proxies = vertex_proxies[nb];
        e.is_boundary = false;

        edge_map[k] = e;
        edge_queue.push(e);
    }
}

// Contract structure-aware edge
bool contract_structure_edge(const StructureEdge &e)
{
    int v1 = e.v1, v2 = e.v2;

    if (!vertex_neighbors[v1].count(v2))
        return false;
    if (!validate_structure_edge_collapse(v1, v2))
    {
        return false;
    }

    // Update vertex position
    vertices[v1] = e.optimal_pos;

    // Merge quadrics
    quadrics[v1] += quadrics[v2];

    // Merge proxy lists
    for (int proxy : vertex_proxies[v2])
    {
        if (find(vertex_proxies[v1].begin(), vertex_proxies[v1].end(), proxy) == vertex_proxies[v1].end())
        {
            vertex_proxies[v1].push_back(proxy);
        }
    }
    vertex_proxies[v2].clear();

    // Update adjacency
    for (int nb : vertex_neighbors[v2])
    {
        if (nb == v1)
            continue;
        vertex_neighbors[v1].insert(nb);
        vertex_neighbors[nb].erase(v2);
        vertex_neighbors[nb].insert(v1);
    }
    vertex_neighbors[v1].erase(v2);
    vertex_neighbors[v2].clear();

    // Update faces and proxy counts
    vector<ivec3> newF;
    newF.reserve(faces.size());

    for (int i = 0; i < faces.size(); ++i)
    {
        ivec3 nf = faces[i];
        if (nf.x == v2)
            nf.x = v1;
        if (nf.y == v2)
            nf.y = v1;
        if (nf.z == v2)
            nf.z = v1;

        if (nf.x != nf.y && nf.y != nf.z && nf.x != nf.z)
        {
            newF.push_back(nf);
        }
        else
        {
            // Face collapsed, update proxy count
            if (!face_proxies[i].empty())
            {
                int proxy_id = face_proxies[i][0];
                if (proxy_face_count.count(proxy_id))
                {
                    proxy_face_count[proxy_id]--;
                }
            }
        }
    }

    faces.swap(newF);

    // Update face_proxies vector
    vector<vector<int>> new_face_proxies;
    new_face_proxies.reserve(faces.size());
    int old_face_idx = 0;

    for (int i = 0; i < faces.size(); ++i)
    {
        // Find corresponding old face
        while (old_face_idx < face_proxies.size())
        {
            ivec3 old_face = faces[i]; // This is already the new face
            new_face_proxies.push_back(face_proxies[old_face_idx]);
            old_face_idx++;
            break;
        }
    }
    face_proxies = new_face_proxies;

    update_structure_edges(v1);
    return true;
}

// Get k random edges with structure awareness
StructureEdge *get_k_structure_edges(int k)
{
    if (k > (int)vertex_neighbors.size())
        k = (int)vertex_neighbors.size();

    vector<StructureEdge *> candidate_edges;

    for (int i = 0; i < k && i < MAX_ITERATIONS; ++i)
    {
        int v1 = uniform_int_distribution<int>(0, vertices.size() - 1)(rng);
        int v2 = uniform_int_distribution<int>(0, vertices.size() - 1)(rng);

        if (v1 == v2)
        {
            i--;
            continue;
        }

        int64_t key = edge_key(v1, v2);
        if (edge_map.count(key) && validate_structure_edge_collapse(v1, v2))
        {
            candidate_edges.push_back(&edge_map[key]);
        }
    }

    if (candidate_edges.empty())
        return nullptr;

    // Find edge with minimum cost
    StructureEdge *best = candidate_edges[0];
    for (size_t i = 1; i < candidate_edges.size(); ++i)
    {
        if (candidate_edges[i]->cost < best->cost)
        {
            best = candidate_edges[i];
        }
    }

    return best;
}

// Cleanup unused vertices
void cleanup_vertices()
{
    id_map.clear();
    int nV = vertices.size();
    vector<bool> used(nV, false);

    for (const auto &f : faces)
    {
        used[f.x] = used[f.y] = used[f.z] = true;
    }

    vector<vec3> newV;
    vector<vector<int>> new_vertex_proxies;
    vector<int> remap(nV, -1);

    for (int i = 0; i < nV; ++i)
    {
        if (used[i])
        {
            remap[i] = newV.size();
            id_map[i] = newV.size();
            newV.push_back(vertices[i]);
            new_vertex_proxies.push_back(vertex_proxies[i]);
        }
    }

    for (auto &f : faces)
    {
        f.x = remap[f.x];
        f.y = remap[f.y];
        f.z = remap[f.z];
    }

    vertices.swap(newV);
    vertex_proxies.swap(new_vertex_proxies);
}

// Main structure-aware simplification function
void structure_aware_simplify(float targetRatio, bool enable_structure_preservation = true)
{
    shape_preserve_flag = enable_structure_preservation;
    modify_cost = enable_structure_preservation;

    initialize_proxies();
    calculate_planar_proxy_equations();
    init_structure_quadrics();
    build_adjacency();
    init_structure_edges();

    int targetFaces = int(faces.size() * (1.0f - targetRatio));

    while (faces.size() > targetFaces && !edge_queue.empty())
    {
        StructureEdge e = edge_queue.top();
        edge_queue.pop();

        int64_t k = edge_key(e.v1, e.v2);
        if (!edge_map.count(k) || fabs(edge_map[k].cost - e.cost) > 1e-9)
        {
            continue;
        }

        if (contract_structure_edge(e))
        {
            edge_map.erase(k);
        }
    }

    cleanup_vertices();
}

// Multiple choice decimation with structure awareness
void multiple_choice_structure_decimation(int k, int target_collapses)
{
    shape_preserve_flag = false;
    shape_preserve_decimation_flag = false;
    modify_cost = false;

    initialize_proxies();
    init_structure_quadrics();
    build_adjacency();

    for (int i = 0; i < target_collapses; ++i)
    {
        if (faces.size() <= 6)
            break;

        StructureEdge *edge = get_k_structure_edges(k);
        if (!edge)
            break;

        if (contract_structure_edge(*edge))
        {
            int64_t key = edge_key(edge->v1, edge->v2);
            edge_map.erase(key);
        }
    }

    cleanup_vertices();
}

// Shape-preserving decimation
void shape_preserve_structure_decimation(int k, int target_collapses)
{
    shape_preserve_flag = true;
    shape_preserve_decimation_flag = true;
    modify_cost = true;

    initialize_proxies();
    calculate_planar_proxy_equations();
    init_structure_quadrics();
    build_adjacency();

    for (int i = 0; i < target_collapses; ++i)
    {
        if (faces.size() <= 6)
            break;

        StructureEdge *edge = get_k_structure_edges(k);
        if (!edge)
            break;

        if (contract_structure_edge(*edge))
        {
            int64_t key = edge_key(edge->v1, edge->v2);
            edge_map.erase(key);
        }
    }

    cleanup_vertices();
}

// Single step structure-aware simplification
void structure_aware_simplify_step()
{
    shape_preserve_flag = true;
    modify_cost = true;

    initialize_proxies();
    calculate_planar_proxy_equations();
    init_structure_quadrics();
    build_adjacency();
    init_structure_edges();

    int targetFaces = int(faces.size() - 1);

    while (faces.size() > targetFaces && !edge_queue.empty())
    {
        StructureEdge e = edge_queue.top();
        edge_queue.pop();

        int64_t k = edge_key(e.v1, e.v2);
        if (!edge_map.count(k) || fabs(edge_map[k].cost - e.cost) > 1e-9)
        {
            continue;
        }

        if (contract_structure_edge(e))
        {
            edge_map.erase(k);
            best_edge = e;
        }
    }

    cleanup_vertices();
}

// External C interface
extern "C"
{
    int SetMesh(vec3 *inV, ivec3 *inF, int nV, int nF)
    {
        vertices.assign(inV, inV + nV);
        faces.assign(inF, inF + nF);
        return 0;
    }

    vec3 *GetMeshVertices(int *outV)
    {
        *outV = vertices.size();
        static vector<vec3> buf;
        buf = vertices;
        return buf.data();
    }

    ivec3 *GetMeshFaces(int *outF)
    {
        *outF = faces.size();
        static vector<ivec3> buf;
        buf = faces;
        return buf.data();
    }

    void Simplify(float targetRatio)
    {
        structure_aware_simplify(targetRatio);
    }

    void SimplifyStep()
    {
        structure_aware_simplify_step();
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

#endif // __STRUCTURE_AWARE_SIMPLIFICATION_H__