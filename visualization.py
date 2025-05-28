import sys
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtGui import QFont
import vtk
import os
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import trimesh
import numpy as np
from lib.interface import *

def set_face_colors(polydata, face_indices, color):
    """
    polydata: vtkPolyData对象
    face_indices: 需要变色的面索引列表（比如 [3, 7, 15]）
    color: RGB三元组，取值0~255 (比如红色(255, 0, 0))
    """
    # 创建颜色数组（每个面4个字节：RGBA）
    num_faces = polydata.GetNumberOfCells()
    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(4)  # RGBA
    colors.SetNumberOfTuples(num_faces)

    # 默认颜色（白色不透明）
    default_color = [224, 224, 255, 255]
    for i in range(num_faces):
        colors.SetTuple(i, default_color)

    # 给指定面的颜色设为color，颜色也需加上Alpha通道255
    for face_id in face_indices:
        colors.SetTuple(face_id, list(color) + [255])

    # 赋值给polydata的CellData
    polydata.GetCellData().SetScalars(colors)

class VTKWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.test_id = 0
        self.mesh1_path = None
        self.mesh2_path = None  # 右侧原始mesh路径

        self.frame = QtWidgets.QFrame()
        self.layout = QtWidgets.QVBoxLayout()

        self.uid_label = QtWidgets.QLabel("请输入物体的id")
        font = QFont("Times new roman", 14)
        self.uid_label.setFont(font)
        # 输入框和按钮
        self.uid_input = QtWidgets.QLineEdit()
        self.uid_input.setPlaceholderText("输入 Object UID，例如 bunny")
        self.uid_input.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.uid_input.setFixedHeight(60)

        # 按钮
        self.load_left_btn = QtWidgets.QPushButton("加载网格")
        self.load_left_btn.clicked.connect(self.on_load_left_clicked)
        self.load_right_btn = QtWidgets.QPushButton("执行简化")
        self.load_right_btn.setEnabled(False)
        self.load_right_btn.clicked.connect(self.load_right_mesh)
        
        self.load_left_btn.setFixedHeight(80)
        self.load_left_btn.setStyleSheet("background-color: lightgreen; font-weight: bold;")
        self.load_right_btn.setFixedHeight(80)
        self.load_right_btn.setStyleSheet("background-color: lightblue; font-weight: bold;")

        self.reset_camera_btn = QtWidgets.QPushButton("重置相机")
        self.reset_camera_btn.clicked.connect(self.reset_camera)

        
        # 状态显示标签
        self.origin_label = QtWidgets.QLabel("")
        font = QFont("Comic Sans MS", 14)
        self.origin_label.setFont(font)

        self.simplified_label = QtWidgets.QLabel("")
        font = QFont("Comic Sans MS", 14)
        self.simplified_label.setFont(font)

        # VTK渲染窗口
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)

        # 布局顺序：按钮、标签、渲染窗口
        # 主水平布局（分左侧和右侧）
        main_layout = QtWidgets.QHBoxLayout()
        left_layout = QtWidgets.QVBoxLayout()
        # 第一行：输入框
        left_layout.addWidget(self.uid_label)
        left_layout.addWidget(self.uid_input)
        # 第二行：加载按钮 + 简化按钮（水平排列）
        button_row = QtWidgets.QHBoxLayout()
        button_row.addWidget(self.load_left_btn)
        button_row.addWidget(self.load_right_btn)
        left_layout.addLayout(button_row)
        # 第三行：状态标签
        left_layout.addWidget(self.origin_label)
        left_layout.addWidget(self.simplified_label)
        left_layout.addWidget(self.reset_camera_btn)
        # 添加 stretch 拉伸填充，避免控件顶到底部
        left_layout.addStretch()
        # 右侧：VTK 渲染窗口
        vtk_layout = QtWidgets.QVBoxLayout()
        vtk_layout.addWidget(self.vtkWidget)
        # 将左右布局加到主布局中
        main_layout.addLayout(left_layout, 1)
        main_layout.addLayout(vtk_layout, 4)
        # 设置主布局到 frame
        self.frame.setLayout(main_layout)
        self.setCentralWidget(self.frame)


        self.renderWindow = self.vtkWidget.GetRenderWindow()

        # 创建左右两个Renderer
        self.renderer_left = vtk.vtkRenderer()
        self.renderer_right = vtk.vtkRenderer()

        self.renderer_left.SetViewport(0.0, 0.0, 0.5, 1.0)
        self.renderer_right.SetViewport(0.5, 0.0, 1.0, 1.0)

        # 共享一个相机
        shared_camera = vtk.vtkCamera()
        self.renderer_left.SetActiveCamera(shared_camera)
        self.renderer_right.SetActiveCamera(shared_camera)

        self.renderWindow.AddRenderer(self.renderer_left)
        self.renderWindow.AddRenderer(self.renderer_right)

        self.renderer_left.SetBackground(1,1,1)
        self.renderer_right.SetBackground(1,1,1)

        # 重置相机，自动适应左侧mesh大小
        self.renderer_left.ResetCamera()

        self.vtkWidget.Initialize()
        self.vtkWidget.Start()
        self.show()

    def on_load_left_clicked(self):
        uid = self.uid_input.text().strip()
        
        # 清空右侧视图
        self.renderer_left.RemoveAllViewProps()
        self.renderer_right.RemoveAllViewProps()
        if not uid:
            self.origin_label.setText("请输入有效的 Object UID")
            return

        self.mesh1_path = os.path.join("objs", uid + ".obj")
        self.mesh2_path = os.path.join("results", uid + ".obj")

        if not os.path.exists(self.mesh1_path):
            self.origin_label.setText(f"找不到文件: {self.mesh1_path}")
            return

        mesh = trimesh.load(self.mesh1_path, process=False)
        V = mesh.vertices.astype(np.float32)
        F = mesh.faces.astype(np.int32)
        self.origin_label.setText(f"原始网格，顶点数: {len(V)}, 面数: {len(F)}")
        self.load_left_mesh(self.mesh1_path)


        self.renderWindow.Render()

        self.load_right_btn.setEnabled(True)

    def load_left_mesh(self, mesh_path):
        self.renderer_left.RemoveAllViewProps()

        reader = vtk.vtkOBJReader()
        reader.SetFileName(mesh_path)
        reader.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(reader.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.8, 0.9, 0.8)
        actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetEdgeColor(0.1, 0.1, 0.1)
        actor.GetProperty().SetLineWidth(1.0)

        self.renderer_left.AddActor(actor)
        self.renderer_left.ResetCamera()
        self.renderWindow.Render()

    def numpy_to_vtk_polydata(self, vertices, faces):
        # 转换numpy数组到vtkPolyData
        points = vtk.vtkPoints()
        for v in vertices:
            points.InsertNextPoint(v[0], v[1], v[2])

        polys = vtk.vtkCellArray()
        for f in faces:
            polys.InsertNextCell(3, f)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetPolys(polys)
        return polydata

    def load_right_mesh(self):
        # 禁用按钮防止重复点击
        self.load_right_btn.setEnabled(False)

        self.simplified_label.setText("正在执行右侧网格简化，请稍候...")
        self.simplified_label.repaint()  # 强制刷新界面显示

        # 加载原始mesh
        mesh = trimesh.load(self.mesh1_path, process=False)
        V = mesh.vertices.astype(np.float32)
        F = mesh.faces.astype(np.int32)

        # 执行简化
        SetMesh(V, F)
        SimplifyStep()
        V_simplified, F_simplified = GetMesh()
        mesh_simplified = trimesh.Trimesh(vertices=V_simplified, faces=F_simplified)
    
        # 保存简化后的网格
        mesh_simplified.export(self.mesh2_path)

        # 更新状态标签
        self.simplified_label.setText(f"简化完成，顶点数: {len(V_simplified)}, 面数: {len(F_simplified)}\n\
        顶点减少: {len(V) - len(V_simplified)} 面减少: {len(F) - len(F_simplified)}\n \
        压缩比: {len(F_simplified)/len(F):.4f}")

        # 清空右侧renderer所有actor
        self.renderer_right.RemoveAllViewProps()

        # 转换简化结果到vtkPolyData
        polydata = self.numpy_to_vtk_polydata(V_simplified, F_simplified)
        

        mapper2 = vtk.vtkPolyDataMapper()
        mapper2.SetInputData(polydata)
        actor2 = vtk.vtkActor()
        actor2.SetMapper(mapper2)

        actor2.GetProperty().SetColor(0.8, 0.8, 0.9)  # 浅蓝色
        actor2.GetProperty().EdgeVisibilityOn()
        actor2.GetProperty().SetEdgeColor(0.1, 0.1, 0.1)  # 深灰色
        actor2.GetProperty().SetLineWidth(1.0)

        self.renderer_right.AddActor(actor2)

        # 渲染更新
        self.renderWindow.Render()

        # 恢复按钮
        self.load_right_btn.setEnabled(True)

    def animate_camera_to(self, target_position, target_focal_point, duration=1.0, steps=60, on_finish=None):
        camera = self.renderer_left.GetActiveCamera()
        self._anim_camera = camera  # 缓存，供 update 调用

        self._anim_on_finish = on_finish  # 保存回调

        # 起始位置
        start_pos = np.array(camera.GetPosition())
        start_fp = np.array(camera.GetFocalPoint())

        self.camera_interp_step = 0
        self.camera_interp_steps = steps
        self.camera_pos_delta = (np.array(target_position) - start_pos) / steps
        self.camera_fp_delta = (np.array(target_focal_point) - start_fp) / steps
        self.camera_start_pos = start_pos
        self.camera_start_fp = start_fp

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.update_camera_interpolation)
        self.timer.start(int(duration * 1000 / steps))
        
    def update_camera_interpolation(self):
        step = self.camera_interp_step
        camera = self._anim_camera

        if step >= self.camera_interp_steps:
            self.timer.stop()
            if self._anim_on_finish:
                self._anim_on_finish()  # 动画结束回调
            return

        new_pos = self.camera_start_pos + step * self.camera_pos_delta
        new_fp = self.camera_start_fp + step * self.camera_fp_delta

        camera.SetPosition(*new_pos)
        camera.SetFocalPoint(*new_fp)
        self.renderWindow.Render()

        self.camera_interp_step += 1

    def reset_camera(self):
        # 获取右侧的 renderer（也可以选择左侧）
        renderer = self.renderer_right if self.renderer_right.GetActors().GetNumberOfItems() > 0 else self.renderer_left
        camera = renderer.GetActiveCamera()
        actors = renderer.GetActors()
        actors.InitTraversal()
        actor = actors.GetNextActor()
        if not actor:
            return
        polydata = actor.GetMapper().GetInput()
        if not polydata:
            return
        face_id = self.test_id
        
        cell = polydata.GetCell(face_id)
        p0 = np.array(polydata.GetPoint(cell.GetPointId(0)))
        p1 = np.array(polydata.GetPoint(cell.GetPointId(1)))
        p2 = np.array(polydata.GetPoint(cell.GetPointId(2)))
        face_center = (p0 + p1 + p2) / 3

        # 获取法向量（用于相机方向）
        normal = np.cross(p1 - p0, p2 - p0)
        normal = normal / np.linalg.norm(normal)

        # 设置相机的位置：在法线方向上方离面中心一定距离
        distance = 1.5 * np.linalg.norm(polydata.GetLength())  # 可调
        camera_position = face_center + normal * distance
        set_face_colors(polydata, [self.test_id], (255, 128, 128))
        self.test_id += 1

        # 设置相机参数
        def on_camera_done():
            camera.SetViewUp(0, 1, 0)
            renderer.ResetCameraClippingRange()
            self.renderWindow.Render()

        self.animate_camera_to(camera_position, face_center, on_finish=on_camera_done)

if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    window = VTKWindow()
    sys.exit(app.exec_())
