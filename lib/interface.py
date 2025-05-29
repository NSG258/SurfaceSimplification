from ctypes import *
import numpy as np
import os

libpath = os.path.dirname(os.path.abspath(__file__))
Processor = cdll.LoadLibrary(os.path.join(libpath, 'main.dll'))
# Processor = cdll.LoadLibrary(os.path.join(libpath, 'main.so'))



def SetMesh(V, F):
  print(f"Vertices count: {V.shape[0]}\nFaces count: {F.shape[0]}")
  handle = Processor.SetMesh(c_void_p(V.ctypes.data), c_void_p(F.ctypes.data), V.shape[0], F.shape[0]);
  return handle

def Simplify(targetRatio=0.5):
  '''
    一步到位，直接获取简化结果
  '''
  Processor.Simplify.argtypes = [c_float]
  Processor.Simplify(targetRatio)
  vertices, faces = GetMesh()
  return vertices, faces
  
def SimplifyStep():
  '''
    SimplifyStep(): 每次调用的时候，只收缩一次
    GetActiveVertices(): 返回收缩的那两个顶点
  '''
  Processor.SimplifyStep()
  
  vertices, faces = GetMesh()
  
  Processor.GetActiveVertices.restype = POINTER(c_int)
  # Processor.GetActiveVertices.argtypes = [POINTER(c_int)]
  # num_active_v = c_int()
  ptr_active_v = Processor.GetActiveVertices()
  arr_active_v = np.ctypeslib.as_array(ptr_active_v, shape=(3,))

  return vertices, faces, arr_active_v

def GetMesh():
  Processor.GetMeshVertices.restype = POINTER(c_float)
  Processor.GetMeshVertices.argtypes = [POINTER(c_int)]
  Processor.GetMeshFaces.restype = POINTER(c_int)
  Processor.GetMeshFaces.argtypes = [POINTER(c_int)]
  
  num_v = c_int()
  ptr_v = Processor.GetMeshVertices(byref(num_v))
  arr_v = np.ctypeslib.as_array(ptr_v, shape=(num_v.value * 3,))
  arr_v = arr_v.reshape((num_v.value, 3))
  
  num_f = c_int()
  ptr_f = Processor.GetMeshFaces(byref(num_f))
  arr_f = np.ctypeslib.as_array(ptr_f, shape=(num_f.value * 3,))
  arr_f = arr_f.reshape((num_f.value, 3))
  return arr_v.copy(), arr_f.copy()
