import os
import numpy as np

# 第一步：加载原始非流形 obj 文件
input_path = os.path.join("non_manifold", "crate.obj")
output_path = os.path.join("objs", "crate.obj")

epsilon = 1e-7

# 读取 .obj 文件
vertices = []
faces = []
with open(input_path, 'r') as f:
    for line in f:
        if line.startswith('v '):
            parts = list(map(float, line.strip().split()[1:]))
            vertices.append(parts)
        elif line.startswith('f '):
            parts = line.strip().split()[1:]
            face = [int(p.split('/')[0]) for p in parts]  # 支持 v/vt/vn 格式
            faces.append(face)

vertices = np.array(vertices)
faces = np.array(faces)

# 建立空间位置哈希（合并重复顶点）
rounded = np.round(vertices / epsilon).astype(np.int64)
pos_to_new_index = {}
new_vertices = []
old_to_new = {}

for i, pos in enumerate(rounded):
    key = tuple(pos)
    if key not in pos_to_new_index:
        pos_to_new_index[key] = len(new_vertices)
        new_vertices.append(vertices[i])
    old_to_new[i] = pos_to_new_index[key]

# 重建 face 索引
new_faces = []
for face in faces:
    new_face = [old_to_new[idx - 1] + 1 for idx in face]  # obj 索引从 1 开始
    if len(set(new_face)) == len(new_face):  # 避免塌陷面（如两个点合并）
        new_faces.append(new_face)

# 写出新的 obj 文件
with open(output_path, 'w') as f:
    for v in new_vertices:
        f.write(f'v {v[0]} {v[1]} {v[2]}\n')
    for face in new_faces:
        f.write('f ' + ' '.join(map(str, face)) + '\n')

print(f'✅ 顶点合并完成，输出文件: {output_path}')
print(f'原始顶点数: {len(vertices)}, 合并后: {len(new_vertices)}')
print(f'原始面数: {len(faces)}, 合并后: {len(new_faces)}')
