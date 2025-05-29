#ifndef LINDSTROM_TURK_SIMPLIFICATION_CPP
#define LINDSTROM_TURK_SIMPLIFICATION_CPP

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>

using namespace std;
using namespace glm;

// 边结构定义
struct Edge
{
    int v1, v2;
    double cost;
    vec3 optimal_pos;
    bool operator>(const Edge &o) const { return cost > o.cost; }
};

// 全局数据
static vector<vec3> vertices;
static vector<ivec3> faces;
static vector<mat4> quadrics;
static vector<vec3> vertex_normals;
static priority_queue<Edge, vector<Edge>, greater<Edge>> edge_queue;
static unordered_map<int64_t, Edge> edge_map;
static vector<unordered_set<int>> vertex_neighbors;
static vector<unordered_set<int>> vertex_faces;
static Edge best_edge;
static unordered_map<int, int> id_map; // <old_id, new_id>

// Lindstrom-Turk 参数
static float volume_weight = 1000.0f;   // 体积保持权重
static float boundary_weight = 1000.0f; // 边界保持权重

// 计算唯一边键
inline int64_t edge_key(int a, int b)
{
    if (a > b)
        std::swap(a, b);
    return ((int64_t)a << 32) | (uint32_t)b;
}

// 计算面法向量和面积
pair<vec3, float> compute_face_normal_area(const vec3 &v0, const vec3 &v1, const vec3 &v2)
{
    vec3 e1 = v1 - v0;
    vec3 e2 = v2 - v0;
    vec3 normal = cross(e1, e2);
    float area = length(normal) * 0.5f;
    if (area > 1e-10f)
        normal = normalize(normal);
    else
        normal = vec3(0.0f);
    return {normal, area};
}

// 计算面基础二次误差矩阵
mat4 compute_face_quadric(const vec3 &v0, const vec3 &v1, const vec3 &v2)
{
    auto [n, area] = compute_face_normal_area(v0, v1, v2);
    float d = -dot(n, v0);
    vec4 plane(n, d);
    mat4 Q(0.0f);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            Q[i][j] = plane[i] * plane[j] * area;
    return Q;
}

// 检查边是否在边界上
bool is_boundary_edge(int v1, int v2)
{
    int shared_faces = 0;
    for (int face_idx : vertex_faces[v1])
    {
        if (vertex_faces[v2].count(face_idx))
            shared_faces++;
    }
    return shared_faces == 1;
}

// 计算体积保持项
double compute_volume_constraint(int v1, int v2, const vec3 &new_pos)
{
    double volume_change = 0.0;
    vec3 old_pos1 = vertices[v1];
    vec3 old_pos2 = vertices[v2];

    // 遍历包含v1或v2的所有面
    unordered_set<int> affected_faces;
    for (int face_idx : vertex_faces[v1])
        affected_faces.insert(face_idx);
    for (int face_idx : vertex_faces[v2])
        affected_faces.insert(face_idx);

    for (int face_idx : affected_faces)
    {
        const ivec3 &face = faces[face_idx];
        vec3 p0 = vertices[face.x];
        vec3 p1 = vertices[face.y];
        vec3 p2 = vertices[face.z];

        // 计算原始体积贡献
        double old_vol = dot(p0, cross(p1, p2)) / 6.0;

        // 计算收缩后体积贡献
        if (face.x == v1 || face.x == v2)
            p0 = new_pos;
        if (face.y == v1 || face.y == v2)
            p1 = new_pos;
        if (face.z == v1 || face.z == v2)
            p2 = new_pos;

        double new_vol = dot(p0, cross(p1, p2)) / 6.0;
        volume_change += abs(new_vol - old_vol);
    }

    return volume_change;
}

// 计算边界保持项
double compute_boundary_constraint(int v1, int v2, const vec3 &new_pos)
{
    if (!is_boundary_edge(v1, v2))
        return 0.0;

    // 对于边界边，倾向于保持在边界上
    vec3 edge_mid = 0.5f * (vertices[v1] + vertices[v2]);
    return length(new_pos - edge_mid);
}

// 计算Lindstrom-Turk边收缩代价
pair<vec3, double> compute_lindstrom_turk_cost(int v1, int v2)
{
    mat4 Q = quadrics[v1] + quadrics[v2];

    // 尝试多个候选位置
    vector<vec3> candidates;
    candidates.push_back(vertices[v1]);                         // 保持v1
    candidates.push_back(vertices[v2]);                         // 保持v2
    candidates.push_back(0.5f * (vertices[v1] + vertices[v2])); // 中点

    // 尝试最优位置求解
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
        candidates.push_back(vec3(res));
    }

    // 评估所有候选位置
    vec3 best_pos = candidates[0];
    double min_cost = 1e20;

    for (const vec3 &pos : candidates)
    {
        // 基础二次误差
        vec4 vH(pos, 1.0f);
        double quad_cost = double(dot(vH, Q * vH));

        // 体积保持约束
        double vol_cost = compute_volume_constraint(v1, v2, pos);

        // 边界保持约束
        double boundary_cost = compute_boundary_constraint(v1, v2, pos);

        // 总代价
        double total_cost = quad_cost + volume_weight * vol_cost + boundary_weight * boundary_cost;

        if (total_cost < min_cost)
        {
            min_cost = total_cost;
            best_pos = pos;
        }
    }

    return {best_pos, min_cost};
}

// 计算顶点法向量
void compute_vertex_normals()
{
    int nV = vertices.size();
    vertex_normals.assign(nV, vec3(0.0f));

    for (int face_idx = 0; face_idx < faces.size(); ++face_idx)
    {
        const ivec3 &face = faces[face_idx];
        auto [normal, area] = compute_face_normal_area(
            vertices[face.x], vertices[face.y], vertices[face.z]);

        vertex_normals[face.x] += normal * area;
        vertex_normals[face.y] += normal * area;
        vertex_normals[face.z] += normal * area;
    }

    for (int i = 0; i < nV; ++i)
    {
        if (length(vertex_normals[i]) > 1e-10f)
            vertex_normals[i] = normalize(vertex_normals[i]);
    }
}

// 初始化顶点二次误差矩阵
void init_quadrics()
{
    int nV = vertices.size();
    quadrics.assign(nV, mat4(0.0f));

    for (auto &f : faces)
    {
        mat4 Qf = compute_face_quadric(vertices[f.x], vertices[f.y], vertices[f.z]);
        quadrics[f.x] += Qf;
        quadrics[f.y] += Qf;
        quadrics[f.z] += Qf;
    }
}

// 建立邻接表和面-顶点关系
void build_adjacency()
{
    int nV = vertices.size();
    vertex_neighbors.assign(nV, {});
    vertex_faces.assign(nV, {});

    for (int face_idx = 0; face_idx < faces.size(); ++face_idx)
    {
        const ivec3 &f = faces[face_idx];

        // 建立邻接关系
        vertex_neighbors[f.x].insert(f.y);
        vertex_neighbors[f.x].insert(f.z);
        vertex_neighbors[f.y].insert(f.x);
        vertex_neighbors[f.y].insert(f.z);
        vertex_neighbors[f.z].insert(f.x);
        vertex_neighbors[f.z].insert(f.y);

        // 建立面-顶点关系
        vertex_faces[f.x].insert(face_idx);
        vertex_faces[f.y].insert(face_idx);
        vertex_faces[f.z].insert(face_idx);
    }
}

// 初始化边优先队列
void init_edges()
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
                auto [optP, cost] = compute_lindstrom_turk_cost(v, nb);
                Edge e{v, nb, cost, optP};
                int64_t k = edge_key(v, nb);
                edge_map[k] = e;
                edge_queue.push(e);
            }
        }
    }
}

// 更新给定顶点相关边
void update_edges(int v)
{
    for (int nb : vertex_neighbors[v])
    {
        int64_t k = edge_key(v, nb);
        auto [optP, cost] = compute_lindstrom_turk_cost(v, nb);
        Edge e{v, nb, cost, optP};
        edge_map[k] = e;
        edge_queue.push(e);
    }
}

// 收缩边
bool contract_edge(const Edge &e)
{
    int v1 = e.v1, v2 = e.v2;
    if (!vertex_neighbors[v1].count(v2))
        return false;

    // 更新顶点位置
    vertices[v1] = e.optimal_pos;
    quadrics[v1] += quadrics[v2];

    // 更新邻接关系
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

    // 更新面-顶点关系
    for (int face_idx : vertex_faces[v2])
    {
        vertex_faces[v1].insert(face_idx);
    }
    vertex_faces[v2].clear();

    // 更新面
    vector<ivec3> newF;
    newF.reserve(faces.size());
    for (auto &f : faces)
    {
        ivec3 nf = f;
        if (nf.x == v2)
            nf.x = v1;
        if (nf.y == v2)
            nf.y = v1;
        if (nf.z == v2)
            nf.z = v1;

        // 去除退化面
        if (nf.x != nf.y && nf.y != nf.z && nf.x != nf.z)
            newF.push_back(nf);
    }
    faces.swap(newF);

    // 重新建立面-顶点关系
    build_adjacency();

    // 更新相关边
    update_edges(v1);

    return true;
}

// 清理未使用顶点
void cleanup()
{
    id_map.clear();
    int nV = vertices.size();
    vector<bool> used(nV, false);
    for (auto &f : faces)
    {
        used[f.x] = used[f.y] = used[f.z] = true;
    }

    vector<vec3> newV;
    newV.reserve(nV);
    vector<int> remap(nV, -1);
    for (int i = 0; i < nV; ++i)
    {
        if (used[i])
        {
            remap[i] = newV.size();
            id_map[i] = newV.size();
            newV.push_back(vertices[i]);
        }
    }

    for (auto &f : faces)
    {
        f.x = remap[f.x];
        f.y = remap[f.y];
        f.z = remap[f.z];
    }
    vertices.swap(newV);
}

// 主简化函数
void simplify_lindstrom_turk(float targetRatio)
{
    compute_vertex_normals();
    init_quadrics();
    build_adjacency();
    init_edges();

    int targetFaces = int(faces.size() * (1.0f - targetRatio));
    while (faces.size() > targetFaces && !edge_queue.empty())
    {
        Edge e = edge_queue.top();
        edge_queue.pop();
        int64_t k = edge_key(e.v1, e.v2);
        if (!edge_map.count(k) || fabs(edge_map[k].cost - e.cost) > 1e-9)
            continue;
        if (contract_edge(e))
            edge_map.erase(k);
    }
    cleanup();
}

void simplify_step_lindstrom_turk()
{
    compute_vertex_normals();
    init_quadrics();
    build_adjacency();
    init_edges();

    if (!edge_queue.empty())
    {
        Edge e = edge_queue.top();
        edge_queue.pop();
        int64_t k = edge_key(e.v1, e.v2);
        if (!edge_map.count(k) || fabs(edge_map[k].cost - e.cost) > 1e-9)
            return;
        if (contract_edge(e))
        {
            edge_map.erase(k);
            best_edge = e;
        }
    }
    cleanup();
}

// 设置Lindstrom-Turk参数
void set_lindstrom_turk_params(float vol_weight, float bound_weight)
{
    volume_weight = vol_weight;
    boundary_weight = bound_weight;
}

// 接口定义
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

    void Simplify(float targetRatio)
    {
        simplify_lindstrom_turk(targetRatio);
    }

    void SimplifyStep()
    {
        simplify_step_lindstrom_turk();
    }

    void SetLindstromTurkParams(float volumeWeight, float boundaryWeight)
    {
        set_lindstrom_turk_params(volumeWeight, boundaryWeight);
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

    ivec3 *GetMeshFaces(int *outF)
    {
        *outF = faces.size();
        static vector<ivec3> buf;
        buf = faces;
        return buf.data();
    }
}

#endif // LINDSTROM_TURK_SIMPLIFICATION_CPP