import numpy as np
from lib.interface import *
import trimesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches
import textwrap

def visualize_mesh_comparison(original_mesh, simplified_mesh):
    """
    在同一界面可视化原始网格和简化后的网格
    """
    fig = plt.figure(figsize=(16, 8))

    # 原始网格
    ax1 = fig.add_subplot(121, projection='3d')
    vertices_orig = original_mesh.vertices
    faces_orig = original_mesh.faces

    # 绘制原始网格的线框
    for face in faces_orig:
        triangle = vertices_orig[face]
        for i in range(3):
            start = triangle[i]
            end = triangle[(i + 1) % 3]
            ax1.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], 
                     'b-', alpha=0.6, linewidth=0.5)

    ax1.set_title(f'Original Mesh\nVertices: {len(vertices_orig)}\nFaces: {len(faces_orig)}',
                  fontsize=11, fontweight='bold')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')

    # 设置相同的视角和范围
    ax1.set_xlim([vertices_orig[:, 0].min(), vertices_orig[:, 0].max()])
    ax1.set_ylim([vertices_orig[:, 1].min(), vertices_orig[:, 1].max()])
    ax1.set_zlim([vertices_orig[:, 2].min(), vertices_orig[:, 2].max()])

    # 简化后的网格
    ax2 = fig.add_subplot(122, projection='3d')
    vertices_simp = simplified_mesh.vertices
    faces_simp = simplified_mesh.faces

    for face in faces_simp:
        triangle = vertices_simp[face]
        for i in range(3):
            start = triangle[i]
            end = triangle[(i + 1) % 3]
            ax2.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], 
                     'r-', alpha=0.6, linewidth=0.5)

    vertex_reduction = (1 - len(vertices_simp) / len(vertices_orig)) * 100
    face_reduction = (1 - len(faces_simp) / len(faces_orig)) * 100

    ax2.set_title(f'Simplified Mesh\nVertices: {len(vertices_simp)} (-{vertex_reduction:.1f}%)\n'
                  f'Faces: {len(faces_simp)} (-{face_reduction:.1f}%)',
                  fontsize=11, fontweight='bold')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')

    ax2.set_xlim([vertices_orig[:, 0].min(), vertices_orig[:, 0].max()])
    ax2.set_ylim([vertices_orig[:, 1].min(), vertices_orig[:, 1].max()])
    ax2.set_zlim([vertices_orig[:, 2].min(), vertices_orig[:, 2].max()])
    ax2.view_init(elev=ax1.elev, azim=ax1.azim)

    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    plt.suptitle('Mesh Simplification Comparison - Wireframe View', fontsize=16, fontweight='bold')

    # 添加统计信息并避免文字重叠
    stats_text = f"""
Simplification Statistics:
• Vertices Reduced: {len(vertices_orig) - len(vertices_simp)} ({vertex_reduction:.1f}%)
• Faces Reduced: {len(faces_orig) - len(faces_simp)} ({face_reduction:.1f}%)
• Compression Ratio: {len(faces_simp)/len(faces_orig):.3f}
    """
    wrapped_text = "\n".join(textwrap.wrap(stats_text, width=70))
    fig.text(0.02, 0.02, wrapped_text, fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8))

    plt.show()
    plt.savefig('mesh_comparison_wireframe.png')

def visualize_mesh_surface_comparison(original_mesh, simplified_mesh):
    """
    表面渲染对比
    """
    fig = plt.figure(figsize=(16, 8))

    # 原始网格
    ax1 = fig.add_subplot(121, projection='3d')
    vertices_orig = original_mesh.vertices
    faces_orig = original_mesh.faces

    ax1.plot_trisurf(vertices_orig[:, 0], vertices_orig[:, 1], vertices_orig[:, 2],
                     triangles=faces_orig, alpha=0.8, cmap='viridis', shade=True)

    ax1.set_title(f'Original Mesh (Surface Rendering)\nVertices: {len(vertices_orig)}\nFaces: {len(faces_orig)}',
                  fontsize=11, fontweight='bold')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')

    # 简化网格
    ax2 = fig.add_subplot(122, projection='3d')
    vertices_simp = simplified_mesh.vertices
    faces_simp = simplified_mesh.faces

    ax2.plot_trisurf(vertices_simp[:, 0], vertices_simp[:, 1], vertices_simp[:, 2],
                     triangles=faces_simp, alpha=0.8, cmap='plasma', shade=True)

    vertex_reduction = (1 - len(vertices_simp) / len(vertices_orig)) * 100
    face_reduction = (1 - len(faces_simp) / len(faces_orig)) * 100

    ax2.set_title(f'Simplified Mesh (Surface Rendering)\nVertices: {len(vertices_simp)} (-{vertex_reduction:.1f}%)\n'
                  f'Faces: {len(faces_simp)} (-{face_reduction:.1f}%)',
                  fontsize=11, fontweight='bold')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')

    # 设置统一坐标范围
    all_vertices = np.vstack([vertices_orig, vertices_simp])
    for ax in [ax1, ax2]:
        ax.set_xlim([all_vertices[:, 0].min(), all_vertices[:, 0].max()])
        ax.set_ylim([all_vertices[:, 1].min(), all_vertices[:, 1].max()])
        ax.set_zlim([all_vertices[:, 2].min(), all_vertices[:, 2].max()])

    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    plt.suptitle('Mesh Simplification Comparison - Surface View', fontsize=16, fontweight='bold')
    plt.show()
    plt.savefig('mesh_comparison_surface.png')

def main():
    print("=== 网格简化程序 - Quadric Error Metrics ===")
    
    # 加载原始网格
    print("正在加载原始网格...")
    mesh_original = trimesh.load('objs/00b2c8c60d2f45a893ee73fd1f107e27.obj', process=False)
    V_original = mesh_original.vertices.astype(np.float32)
    F_original = mesh_original.faces.astype(np.int32)
    
    print(f"原始网格 - 顶点数: {len(V_original)}, 面数: {len(F_original)}")
    
    # 执行网格简化
    print("正在执行网格简化...")
    SetMesh(V_original, F_original)
    V_simplified, F_simplified = GetMesh()
    
    print(f"简化后网格 - 顶点数: {len(V_simplified)}, 面数: {len(F_simplified)}")
    
    # 创建简化后的网格对象
    mesh_simplified = trimesh.Trimesh(vertices=V_simplified, faces=F_simplified)
    
    # 保存简化后的网格
    mesh_simplified.export('new_mesh.obj')
    print("简化后的网格已保存为 'new_mesh.obj'")
    
    # 可视化对比
    print("正在生成可视化...")
    
    # 线框模式对比
    print("显示线框模式对比...")
    visualize_mesh_comparison(mesh_original, mesh_simplified)
    
    # 表面渲染模式对比
    print("显示表面渲染模式对比...")
    visualize_mesh_surface_comparison(mesh_original, mesh_simplified)
    
    # 打印详细统计信息
    vertex_reduction = (1 - len(V_simplified) / len(V_original)) * 100
    face_reduction = (1 - len(F_simplified) / len(F_original)) * 100
    
    print("\n=== 简化统计 ===")
    print(f"原始顶点数: {len(V_original)}")
    print(f"简化后顶点数: {len(V_simplified)}")
    print(f"顶点减少: {len(V_original) - len(V_simplified)} ({vertex_reduction:.2f}%)")
    print(f"原始面数: {len(F_original)}")
    print(f"简化后面数: {len(F_simplified)}")
    print(f"面减少: {len(F_original) - len(F_simplified)} ({face_reduction:.2f}%)")
    print(f"压缩比: {len(F_simplified)/len(F_original):.4f}")

if __name__ == "__main__":
    main()