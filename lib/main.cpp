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
#include <glm/gtc/matrix_inverse.hpp> // for determinant and inverse

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
static priority_queue<Edge, vector<Edge>, greater<Edge>> edge_queue;
static unordered_map<int64_t, Edge> edge_map; // 使用 int64_t
static vector<unordered_set<int>> vertex_neighbors;
static Edge best_edge;

// 计算唯一边键，返回 64 位键值
inline int64_t edge_key(int a, int b)
{
    if (a > b)
        std::swap(a, b);
    return ((int64_t)a << 32) | (uint32_t)b;
}

// 计算面二次误差矩阵
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

// 计算边收缩最优位置和代价
pair<vec3, double> compute_edge_cost(int v1, int v2)
{
    mat4 Q = quadrics[v1] + quadrics[v2];
    mat4 Qp = Q;
    // 清零最后一行前3列, 置(3,3)=1
    for (int col = 0; col < 4; ++col)
    {
        Qp[col][3] = (col == 3) ? 1.0f : 0.0f;
    }
    vec3 optP;
    float det = determinant(Qp);
    if (fabs(det) > 1e-6f)
    {
        mat4 invQ = inverse(Qp);
        vec4 res = invQ * vec4(0, 0, 0, 1);
        optP = vec3(res);
    }
    else
    {
        // 不可逆则取中点
        optP = 0.5f * (vertices[v1] + vertices[v2]);
    }
    vec4 vH(optP, 1.0f);
    double cost = double(dot(vH, Q * vH));
    return {optP, cost};
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

// 建立邻接表
void build_adjacency()
{
    int nV = vertices.size();
    vertex_neighbors.assign(nV, {});
    for (auto &f : faces)
    {
        vertex_neighbors[f.x].insert(f.y);
        vertex_neighbors[f.x].insert(f.z);
        vertex_neighbors[f.y].insert(f.x);
        vertex_neighbors[f.y].insert(f.z);
        vertex_neighbors[f.z].insert(f.x);
        vertex_neighbors[f.z].insert(f.y);
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
                auto [optP, cost] = compute_edge_cost(v, nb);
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
        auto [optP, cost] = compute_edge_cost(v, nb);
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
    vertices[v1] = e.optimal_pos;
    quadrics[v1] += quadrics[v2];
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
        if (nf.x != nf.y && nf.y != nf.z && nf.x != nf.z)
            newF.push_back(nf);
    }
    faces.swap(newF);
    update_edges(v1);
    return true;
}

// 清理未使用顶点
void cleanup()
{
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
void simplify(float targetRatio)
{
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

void simplifyStep()
{
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
        simplify(targetRatio);
    }

    void SimplifyStep()
    {
        simplifyStep();
    }

    int *GetActiveVertices()
    {
        static vector<int> buf;
        buf.push_back(best_edge.v1);
        buf.push_back(best_edge.v2);
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

#endif // MESH_SIMPLIFICATION_CPP
