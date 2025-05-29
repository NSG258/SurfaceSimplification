#ifndef SPECTRAL_SIMPLIFICATION_CPP
#define SPECTRAL_SIMPLIFICATION_CPP

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
static vector<vector<float>> laplacian_matrix;
static vector<vector<float>> eigenvectors;
static vector<float> eigenvalues;
static vector<float> vertex_importance;
static priority_queue<Edge, vector<Edge>, greater<Edge>> edge_queue;
static unordered_map<int64_t, Edge> edge_map;
static vector<unordered_set<int>> vertex_neighbors;
static vector<float> vertex_areas;
static Edge best_edge;
static unordered_map<int, int> id_map; // <old_id, new_id>

// 谱简化参数
static int num_eigenvectors = 10;
static float spectral_weight = 100.0f;
static float geometric_weight = 1.0f;

// 计算唯一边键
inline int64_t edge_key(int a, int b)
{
    if (a > b)
        std::swap(a, b);
    return ((int64_t)a << 32) | (uint32_t)b;
}

// 计算顶点面积权重（使用混合面积）
void compute_vertex_areas()
{
    int nV = vertices.size();
    vertex_areas.assign(nV, 0.0f);

    for (const auto &face : faces)
    {
        vec3 v0 = vertices[face.x];
        vec3 v1 = vertices[face.y];
        vec3 v2 = vertices[face.z];

        vec3 e1 = v1 - v0;
        vec3 e2 = v2 - v0;
        float area = length(cross(e1, e2)) * 0.5f;

        // 分配1/3面积给每个顶点
        vertex_areas[face.x] += area / 3.0f;
        vertex_areas[face.y] += area / 3.0f;
        vertex_areas[face.z] += area / 3.0f;
    }
}

// 构建拉普拉斯矩阵（余切权重）
void build_laplacian_matrix()
{
    int nV = vertices.size();
    laplacian_matrix.assign(nV, vector<float>(nV, 0.0f));

    // 计算余切权重
    for (const auto &face : faces)
    {
        int i = face.x, j = face.y, k = face.z;
        vec3 vi = vertices[i], vj = vertices[j], vk = vertices[k];

        // 计算三个角的余切值
        vec3 eij = vj - vi, eik = vk - vi;
        vec3 eji = vi - vj, ejk = vk - vj;
        vec3 eki = vi - vk, ekj = vj - vk;

        float cot_k = dot(eik, ekj) / length(cross(eik, ekj));
        float cot_i = dot(eij, eik) / length(cross(eij, eik));
        float cot_j = dot(eji, ejk) / length(cross(eji, ejk));

        // 添加权重到拉普拉斯矩阵
        laplacian_matrix[i][j] += cot_k * 0.5f;
        laplacian_matrix[j][i] += cot_k * 0.5f;
        laplacian_matrix[j][k] += cot_i * 0.5f;
        laplacian_matrix[k][j] += cot_i * 0.5f;
        laplacian_matrix[k][i] += cot_j * 0.5f;
        laplacian_matrix[i][k] += cot_j * 0.5f;
    }

    // 设置对角线元素
    for (int i = 0; i < nV; ++i)
    {
        float sum = 0.0f;
        for (int j = 0; j < nV; ++j)
        {
            if (i != j)
                sum += laplacian_matrix[i][j];
        }
        laplacian_matrix[i][i] = -sum;
    }
}

// 幂迭代法计算特征向量
vector<float> power_iteration(const vector<vector<float>> &matrix, int max_iter = 100)
{
    int n = matrix.size();
    vector<float> v(n, 1.0f);

    for (int iter = 0; iter < max_iter; ++iter)
    {
        vector<float> Av(n, 0.0f);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                Av[i] += matrix[i][j] * v[j];
            }
        }

        // 归一化
        float norm = 0.0f;
        for (float x : Av)
            norm += x * x;
        norm = sqrt(norm);

        if (norm > 1e-10f)
        {
            for (int i = 0; i < n; ++i)
                v[i] = Av[i] / norm;
        }
    }

    return v;
}

// 简化的特征分解（仅计算几个主要特征向量）
void compute_eigenvectors()
{
    int nV = vertices.size();
    eigenvectors.clear();
    eigenvalues.clear();

    // 使用修正的拉普拉斯矩阵进行特征分解
    vector<vector<float>> L = laplacian_matrix;

    // 计算前几个特征向量
    for (int k = 0; k < std::min(num_eigenvectors, nV); ++k)
    {
        vector<float> eigenvec = power_iteration(L);
        eigenvectors.push_back(eigenvec);

        // 计算特征值
        float eigenval = 0.0f;
        for (int i = 0; i < nV; ++i)
        {
            float Lv_i = 0.0f;
            for (int j = 0; j < nV; ++j)
            {
                Lv_i += L[i][j] * eigenvec[j];
            }
            eigenval += eigenvec[i] * Lv_i;
        }
        eigenvalues.push_back(eigenval);

        // 从矩阵中移除这个特征向量的贡献
        for (int i = 0; i < nV; ++i)
        {
            for (int j = 0; j < nV; ++j)
            {
                L[i][j] -= eigenval * eigenvec[i] * eigenvec[j];
            }
        }
    }
}

// 计算顶点的谱重要性
void compute_vertex_importance()
{
    int nV = vertices.size();
    vertex_importance.assign(nV, 0.0f);

    // 基于低频特征向量计算重要性
    for (int v = 0; v < nV; ++v)
    {
        float importance = 0.0f;
        for (int k = 0; k < eigenvectors.size(); ++k)
        {
            // 权重与特征值成反比（低频更重要）
            float weight = 1.0f / (abs(eigenvalues[k]) + 1e-6f);
            importance += weight * abs(eigenvectors[k][v]);
        }
        vertex_importance[v] = importance;
    }

    // 归一化重要性值
    float max_importance = *max_element(vertex_importance.begin(), vertex_importance.end());
    if (max_importance > 1e-10f)
    {
        for (float &imp : vertex_importance)
            imp /= max_importance;
    }
}

// 计算谱距离
float compute_spectral_distance(int v1, int v2)
{
    float dist = 0.0f;
    for (int k = 0; k < eigenvectors.size(); ++k)
    {
        float diff = eigenvectors[k][v1] - eigenvectors[k][v2];
        dist += diff * diff;
    }
    return sqrt(dist);
}

// 计算谱边收缩代价
pair<vec3, double> compute_spectral_cost(int v1, int v2)
{
    // 几何代价（基础二次误差）
    vec3 pos1 = vertices[v1];
    vec3 pos2 = vertices[v2];
    vec3 mid_pos = 0.5f * (pos1 + pos2);

    float geometric_cost = length(pos1 - pos2);

    // 谱代价
    float spectral_dist = compute_spectral_distance(v1, v2);
    float importance_penalty = vertex_importance[v1] + vertex_importance[v2];
    float spectral_cost = spectral_dist * importance_penalty;

    // 组合代价
    double total_cost = geometric_weight * geometric_cost +
                        spectral_weight * spectral_cost;

    return {mid_pos, total_cost};
}

// 建立邻接表
void build_adjacency()
{
    int nV = vertices.size();
    vertex_neighbors.assign(nV, {});

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
                auto [optP, cost] = compute_spectral_cost(v, nb);
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
        auto [optP, cost] = compute_spectral_cost(v, nb);
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

    // 更新顶点重要性（简单平均）
    vertex_importance[v1] = (vertex_importance[v1] + vertex_importance[v2]) * 0.5f;

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

    // 重新计算谱信息（简化版本 - 仅在大量简化后重新计算）
    static int contract_count = 0;
    if (++contract_count % 100 == 0)
    {
        compute_vertex_areas();
        build_laplacian_matrix();
        compute_eigenvectors();
        compute_vertex_importance();
    }

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
    vector<float> newImp;
    newV.reserve(nV);
    newImp.reserve(nV);
    vector<int> remap(nV, -1);

    for (int i = 0; i < nV; ++i)
    {
        if (used[i])
        {
            remap[i] = newV.size();
            id_map[i] = newV.size();
            newV.push_back(vertices[i]);
            newImp.push_back(vertex_importance[i]);
        }
    }

    for (auto &f : faces)
    {
        f.x = remap[f.x];
        f.y = remap[f.y];
        f.z = remap[f.z];
    }

    vertices.swap(newV);
    vertex_importance.swap(newImp);
}

// 主谱简化函数
void simplify_spectral(float targetRatio)
{
    compute_vertex_areas();
    build_laplacian_matrix();
    compute_eigenvectors();
    compute_vertex_importance();
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

void simplify_step_spectral()
{
    compute_vertex_areas();
    build_laplacian_matrix();
    compute_eigenvectors();
    compute_vertex_importance();
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

// 设置谱简化参数
void set_spectral_params(int num_eigen, float spec_weight, float geom_weight)
{
    num_eigenvectors = num_eigen;
    spectral_weight = spec_weight;
    geometric_weight = geom_weight;
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
        simplify_spectral(targetRatio);
    }

    void SimplifyStep()
    {
        simplify_step_spectral();
    }

    void SetSpectralParams(int numEigenvectors, float spectralWeight, float geometricWeight)
    {
        set_spectral_params(numEigenvectors, spectralWeight, geometricWeight);
    }

    float *GetVertexImportance(int *outSize)
    {
        *outSize = vertex_importance.size();
        static vector<float> buf;
        buf = vertex_importance;
        return buf.data();
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

#endif // SPECTRAL_SIMPLIFICATION_CPP