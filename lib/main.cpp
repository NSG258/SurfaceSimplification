#include <iostream>
#include <cstring>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

glm::vec3* mesh_vertices = nullptr;
glm::ivec3* mesh_faces = nullptr;

int mesh_vertex_count = 0;
int mesh_face_count = 0;

extern "C" {

int SetMesh(glm::vec3* vertices, glm::ivec3* faces, int num_v, int num_f) {
    if (mesh_vertices) {
        delete[] mesh_vertices;
    }
    if (mesh_faces) {
        delete[] mesh_faces;
    }

    mesh_vertex_count = num_v;
    mesh_face_count = num_f;

    mesh_vertices = new glm::vec3[num_v];
    mesh_faces = new glm::ivec3[num_f];

    std::memcpy(mesh_vertices, vertices, sizeof(glm::vec3) * num_v);
    std::memcpy(mesh_faces, faces, sizeof(glm::ivec3) * num_f);


    return 0;
}

glm::vec3* GetMeshVertices(int* out_num_v) {
    *out_num_v = mesh_vertex_count;
    return mesh_vertices;
}

glm::ivec3* GetMeshFaces(int* out_num_f) {
    *out_num_f = mesh_face_count;
    return mesh_faces;
}

}
