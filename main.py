import numpy as np
from lib.interface import *
import trimesh

mesh = trimesh.load('objs/00d1be4411f848efaeb72b936d4d1692.obj', process=False)
V = mesh.vertices.astype(np.float32)
F = mesh.faces.astype(np.int32)


SetMesh(V, F)
V_new, F_new = GetMesh()

mesh_new = trimesh.Trimesh(vertices=V_new, faces=F_new)
mesh_new.export('new_mesh.obj')