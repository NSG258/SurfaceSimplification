import sys
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtGui import QFont
import vtk
import os
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import trimesh
import numpy as np
from lib.interface import *
import time

def get_active_faces(faces, uids):
    has_v1 = np.any(faces == uids[0], axis=1)
    has_v2 = np.any(faces == uids[1], axis=1)
    
    both = np.where(np.logical_and(has_v1, has_v2))[0]
    only_v1 = np.where(np.logical_and(has_v1, ~has_v2))[0]
    only_v2 = np.where(np.logical_and(~has_v1, has_v2))[0]

    return both, only_v1, only_v2

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
    def __init__(self, lib_suffix, parent=None):
        super().__init__(parent)

        self.origin_mesh_path = None
        self.simplified_mesh_path = None  # 右侧原始mesh路径
        self.display_mode = "result"
        self.current_v = None
        self.current_f = None
        self.uid = None
        self.method = "QEM"
        self.simplify_ratio = 0.5
        self.lib_suffix = lib_suffix

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
        
        # 创建滑动条
        self.ratio_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.ratio_slider.setMinimum(0)
        self.ratio_slider.setMaximum(100)
        self.ratio_slider.setValue(int(self.simplify_ratio * 100))
        self.ratio_slider.setTickInterval(10)
        self.ratio_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)

        # 标签显示数值
        self.ratio_label = QtWidgets.QLabel(f"物体简化比例 {self.simplify_ratio:.2f}")
        font = QFont("Times New Roman", 10)
        self.ratio_label.setFont(font)
        self.ratio_slider.valueChanged.connect(self.update_ratio_live)

        # 按钮
        self.load_left_btn = QtWidgets.QPushButton("加载网格/重置")
        self.load_left_btn.clicked.connect(self.load_left_mesh)
        self.load_right_btn = QtWidgets.QPushButton("执行简化")
        self.load_right_btn.setEnabled(False)
        self.load_right_btn.clicked.connect(self.load_right_mesh)
        
        self.load_left_btn.setFixedHeight(80)
        self.load_left_btn.setFixedWidth(320)
        self.load_left_btn.setStyleSheet("""
            QPushButton {
                background-color: lightgreen;      /* 背景颜色 */
                color: black;                   /* 字体颜色 */
                border: none;                   /* 无边框线 */
                border-radius: 10px;            /* 圆角半径 */
                padding: 6px 12px;              /* 内边距 */
                font-weight: bold;              /* 字体加粗 */
            }
            QPushButton:hover {
                background-color: green;      /* 悬停时颜色 */
            }
            QPushButton:disabled {
                color: #E0E0E0;             /* 禁用状态字体颜色（灰色） */
            }
        """)
        self.load_right_btn.setFixedHeight(80)
        self.load_right_btn.setFixedWidth(320)
        self.load_right_btn.setStyleSheet("""
            QPushButton {
                background-color: lightblue;      /* 背景颜色 */
                color: black;                   /* 字体颜色 */
                border: none;                   /* 无边框线 */
                border-radius: 10px;            /* 圆角半径 */
                padding: 6px 12px;              /* 内边距 */
                font-weight: bold;              /* 字体加粗 */
            }
            QPushButton:hover {
                background-color: blue;      /* 悬停时颜色 */
            }
            QPushButton:disabled {
                color: #E0E0E0;             /* 禁用状态字体颜色（灰色） */
            }
        """)

        self.simplfy_step_btn = QtWidgets.QPushButton("查看单步简化")
        self.simplfy_step_btn.clicked.connect(self.simplfy_step)
        self.simplfy_step_btn.setEnabled(False)
        self.reset_camera_btn = QtWidgets.QPushButton("下一步")
        self.reset_camera_btn.clicked.connect(self.simplify_single_step)
        self.reset_camera_btn.setEnabled(False)
        
        self.simplfy_step_btn.setFixedHeight(80)
        self.simplfy_step_btn.setFixedWidth(320)
        self.simplfy_step_btn.setStyleSheet("""
            QPushButton {
                background-color: orange;      /* 背景颜色 */
                color: black;                   /* 字体颜色 */
                border: none;                   /* 无边框线 */
                border-radius: 10px;            /* 圆角半径 */
                padding: 6px 12px;              /* 内边距 */
                font-weight: bold;              /* 字体加粗 */
            }
            QPushButton:hover {
                background-color: #FF8C00;      /* 悬停时颜色 */
            }
            QPushButton:disabled {
                color: #E0E0E0;             /* 禁用状态字体颜色（灰色） */
            }
        """)
        self.reset_camera_btn.setFixedHeight(80)
        self.reset_camera_btn.setFixedWidth(320)
        self.reset_camera_btn.setStyleSheet("""
            QPushButton {
                background-color: yellow;      /* 背景颜色 */
                color: black;                   /* 字体颜色 */
                border: none;                   /* 无边框线 */
                border-radius: 10px;            /* 圆角半径 */
                padding: 6px 12px;              /* 内边距 */
                font-weight: bold;              /* 字体加粗 */
            }
            QPushButton:hover {
                background-color: #EAB540;      /* 悬停时颜色 */
            }
            QPushButton:disabled {
                color: #E0E0E0;             /* 禁用状态字体颜色（灰色） */
            }
        """)

        self.origin_mesh_label = QtWidgets.QLabel("简化前的网格")
        font = QFont("黑体", 20)
        self.origin_mesh_label.setFont(font)
        
        self.simplified_mesh_label = QtWidgets.QLabel("简化后的网格")
        font = QFont("黑体", 20)
        self.simplified_mesh_label.setFont(font)
        
        # 状态显示标签
        self.origin_label = QtWidgets.QLabel("")
        font = QFont("Times New Roman", 14)
        self.origin_label.setFont(font)

        self.simplified_label = QtWidgets.QLabel("")
        font = QFont("Times New Roman", 14)
        self.simplified_label.setFont(font)
        
        self.export_btn = QtWidgets.QPushButton("导出图像")
        self.export_btn.clicked.connect(self.export_images)
        self.export_btn.setFixedHeight(60)
        self.export_btn.setFixedWidth(320)
        self.export_btn.setStyleSheet("""
            QPushButton {
                background-color: pink;      /* 背景颜色 */
                color: black;                   /* 字体颜色 */
                border: none;                   /* 无边框线 */
                border-radius: 10px;            /* 圆角半径 */
                padding: 6px 12px;              /* 内边距 */
                font-weight: bold;              /* 字体加粗 */
            }
            QPushButton:hover {
                background-color: violet;      /* 悬停时颜色 */
            }
            QPushButton:disabled {
                color: #E0E0E0;             /* 禁用状态字体颜色（灰色） */
            }
        """)
        
        # 添加简化方法选择按钮
        
        self.method_label = QtWidgets.QLabel("")
        font = QFont("Times New Roman", 14)
        self.method_label.setFont(font)
        self.method_label.setText("当前方法: 经典QEM")
        
        self.qem_button = QtWidgets.QPushButton("QEM")
        self.qem_button.setFixedSize(120, 120)
        self.lt_button = QtWidgets.QPushButton("Lindstrom\nTurk")
        self.lt_button.setFixedSize(120, 120)
        self.spectral_button = QtWidgets.QPushButton("Spectral")
        self.spectral_button.setFixedSize(120, 120)
        self.structure_aware_button = QtWidgets.QPushButton("Structure\nAware")
        self.structure_aware_button.setFixedSize(120, 120)
        
        self.other_method_style_sheet = """
            QPushButton {
                background-color: white;
                color: black;
                border: none;
                border-radius: 60px;  /* 一半的大小，形成圆形 */
                font-weight: bold;
            }
        """
        self.method_style_sheet = """
            QPushButton {
                background-color: black;
                color: white;
                border: none;
                border-radius: 60px;  /* 一半的大小，形成圆形 */
                font-weight: bold;
            }
        """
        self.qem_button.setStyleSheet(self.method_style_sheet)
        self.lt_button.setStyleSheet(self.other_method_style_sheet)
        self.spectral_button.setStyleSheet(self.other_method_style_sheet)
        self.structure_aware_button.setStyleSheet(self.other_method_style_sheet)

        self.qem_button.clicked.connect(lambda: self.set_simplify_method("QEM"))
        self.lt_button.clicked.connect(lambda: self.set_simplify_method("LindstromTurk"))
        self.spectral_button.clicked.connect(lambda: self.set_simplify_method("Spectral"))
        self.structure_aware_button.clicked.connect(lambda: self.set_simplify_method("StructureAware"))


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
        button_row1 = QtWidgets.QHBoxLayout()
        button_row1.addWidget(self.load_left_btn)
        button_row1.addWidget(self.load_right_btn)
        left_layout.addLayout(button_row1)
        slider_layout = QtWidgets.QHBoxLayout()
        slider_layout.addWidget(self.ratio_label)
        slider_layout.addWidget(self.ratio_slider)
        left_layout.addLayout(slider_layout)  # 你自己的主布局，比如QVBoxLayout
        # 第三行：状态标签
        left_layout.addWidget(self.origin_label)
        left_layout.addWidget(self.simplified_label)
        
        button_row2 = QtWidgets.QHBoxLayout()
        button_row2.addWidget(self.simplfy_step_btn)
        button_row2.addWidget(self.reset_camera_btn)
        left_layout.addLayout(button_row2)
        
        left_layout.addWidget(self.export_btn)
        
        left_layout.addWidget(self.method_label)

        button_row3  = QtWidgets.QHBoxLayout()
        button_row3.addWidget(self.qem_button)
        button_row3.addWidget(self.lt_button)
        button_row3.addWidget(self.spectral_button)
        button_row3.addWidget(self.structure_aware_button)
        left_layout.addLayout(button_row3)
        # 添加 stretch 拉伸填充，避免控件顶到底部
        left_layout.addStretch()
        # 右侧：VTK 渲染窗口
        vtk_layout = QtWidgets.QVBoxLayout()
        label_row = QtWidgets.QHBoxLayout()
        label_row.addStretch(1)
        label_row.addWidget(self.origin_mesh_label)
        label_row.addStretch(1)
        label_row.addWidget(self.simplified_mesh_label)
        label_row.addStretch(1)
        vtk_layout.addLayout(label_row)
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
    def update_ratio_live(self, value):
        self.simplify_ratio = value / 100.0
        self.ratio_label.setText(f"物体简化比例 {self.simplify_ratio:.2f}")
        
        
    def set_simplify_method(self, method):
        self.method = method
        self.origin_mesh_path = None
        self.simplified_mesh_path = None
        self.display_mode = "result"
        self.current_v = []
        self.current_f = []
        self.uid = None
        self.simplified_label.clear()
        self.origin_label.clear()
        self.load_right_btn.setEnabled(False)
        self.simplfy_step_btn.setEnabled(False)
        self.reset_camera_btn.setEnabled(False)
                
        self.render_mesh(self.renderer_left, self.current_v, self.current_f, (0.8, 0.9, 0.8), face_colors=None)
        self.render_mesh(self.renderer_right, self.current_v, self.current_f, (0.8, 0.8, 0.9), face_colors=None)
        self.renderWindow.Render()
        
        self.qem_button.setStyleSheet(self.other_method_style_sheet)
        self.lt_button.setStyleSheet(self.other_method_style_sheet)
        self.spectral_button.setStyleSheet(self.other_method_style_sheet)
        self.structure_aware_button.setStyleSheet(self.other_method_style_sheet)

        if self.method == "QEM":
            init_processor(f"main{self.lib_suffix}")
            self.qem_button.setStyleSheet(self.method_style_sheet)
            self.method_label.setText("当前方法: 经典QEM")
            self.method_label.repaint()
        elif self.method == "LindstromTurk":
            init_processor(f"lindstrom_turk{self.lib_suffix}")
            self.lt_button.setStyleSheet(self.method_style_sheet)
            self.method_label.setText("当前方法: Lindstrom Turk")
            self.method_label.repaint()
        elif self.method == "Spectral":
            init_processor(f"spectral{self.lib_suffix}")
            self.spectral_button.setStyleSheet(self.method_style_sheet)
            self.method_label.setText("当前方法: Spectral")
            self.method_label.repaint()
        else:
            init_processor(f"structure_aware{self.lib_suffix}")
            self.structure_aware_button.setStyleSheet(self.method_style_sheet)
            self.method_label.setText("当前方法: Structure Aware")
            self.method_label.repaint()
            
        print(f"简化方法已设置为：{method}")
    def export_images(self):
        # 获取当前窗口尺寸，用于截图
        window = self.renderWindow
        w2if = vtk.vtkWindowToImageFilter()
        w2if.SetInput(window)
        w2if.Update()

        # 保存整个窗口截图（可选）
        writer = vtk.vtkPNGWriter()
        writer.SetFileName("full_view.png")
        writer.SetInputConnection(w2if.GetOutputPort())
        writer.Write()

        # 如果你想分开保存左半和右半，可以截取部分区域：
        def save_viewport_image(renderer, filename):
            window = self.renderWindow

            # 截图当前窗口
            w2if = vtk.vtkWindowToImageFilter()
            w2if.SetInput(window)
            w2if.SetInputBufferTypeToRGBA()
            w2if.ReadFrontBufferOff()  # 捕获back buffer，更清晰
            w2if.Update()

            # 获取 renderer 的 viewport 区域
            x0, y0 = int(renderer.GetViewport()[0] * window.GetSize()[0]), int(renderer.GetViewport()[1] * window.GetSize()[1])
            x1, y1 = int(renderer.GetViewport()[2] * window.GetSize()[0]), int(renderer.GetViewport()[3] * window.GetSize()[1])

            # 提取对应区域
            extract = vtk.vtkExtractVOI()
            extract.SetInputConnection(w2if.GetOutputPort())
            extract.SetVOI(x0, x1 - 1, y0, y1 - 1, 0, 0)
            extract.Update()

            writer = vtk.vtkPNGWriter()
            writer.SetFileName(filename)
            writer.SetInputConnection(extract.GetOutputPort())
            writer.Write()

        save_viewport_image(self.renderer_left, os.path.join("images", f"{self.uid}_left.png"))
        save_viewport_image(self.renderer_right, os.path.join("images", f"{self.uid}_{self.method}.png"))

        QtWidgets.QMessageBox.information(self, "导出完成", "图像已保存为 left_view.png 和 right_view.png")
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
    
    def load_left_mesh(self):
        uid = self.uid_input.text().strip()
        # 清空右侧视图
        self.renderer_left.RemoveAllViewProps()
        self.renderer_right.RemoveAllViewProps()
        self.uid = uid
        if not uid:
            self.origin_label.setText("请输入有效的 Object UID")
            return
        self.load_right_btn.setEnabled(False)
        self.reset_camera_btn.setEnabled(False)
        self.simplfy_step_btn.setEnabled(False)
        self.origin_mesh_path = os.path.join("objs", uid + ".obj")
        self.simplified_mesh_path = os.path.join("results", uid + ".obj")
        self.display_mode = "result"
        if not os.path.exists(self.origin_mesh_path):
            self.origin_label.setText(f"找不到文件: {self.origin_mesh_path}")
            return
        
        mesh = trimesh.load(self.origin_mesh_path, process=False)
        
        V = mesh.vertices.astype(np.float32)
        if uid == "gun" or uid == "sword":
            V = V @ np.array([[0, 0, 1.0],
                              [0, 1.0, 0],
                              [1.0, 0, 0]], dtype=V.dtype)
        self.current_v = V
        F = mesh.faces.astype(np.int32)
        self.current_f = F
        self.origin_label.setText(f"原始网格，顶点数: {len(V)}, 面数: {len(F)}")
        
        # 清空右侧renderer所有actor
        self.renderer_left.RemoveAllViewProps()

        # 转换简化结果到vtkPolyData
        polydata = self.numpy_to_vtk_polydata(self.current_v, self.current_f)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.8, 0.9, 0.8)
        actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetEdgeColor(0.1, 0.1, 0.1)
        actor.GetProperty().SetLineWidth(1.0)

        self.renderer_left.AddActor(actor)
        self.renderer_left.ResetCamera()
        self.renderWindow.Render()
        self.load_right_btn.setEnabled(True)
        # self.reset_camera_btn.setEnabled(True)
        self.simplfy_step_btn.setEnabled(True)
        
    def render_mesh(self, renderer, vertices, faces, base_color, face_colors=None):
        renderer.RemoveAllViewProps()

        polydata = self.numpy_to_vtk_polydata(vertices, faces)

        if face_colors is not None:
            # 设置颜色（assume face_colors: dict{face_id: (r, g, b)})
            colors = vtk.vtkUnsignedCharArray()
            colors.SetNumberOfComponents(3)
            colors.SetName("Colors")

            for i in range(len(faces)):
                if i in face_colors:
                    colors.InsertNextTypedTuple(face_colors[i])
                else:
                    base_color_int = tuple(int(x * 255) for x in base_color)
                    colors.InsertNextTypedTuple(base_color_int)  # 默认灰色
            polydata.GetCellData().SetScalars(colors)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        actor.GetProperty().SetColor(*base_color)
        actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetEdgeColor(0.1, 0.1, 0.1)
        actor.GetProperty().SetLineWidth(1.0)

        renderer.AddActor(actor)


    def load_right_mesh(self):
        # 禁用按钮防止重复点击
        if not self.display_mode == "result":
            self.simplified_label.setText("处于单步调试状态，无法一步到位...")
            self.simplified_label.repaint()
            return
        
        self.load_right_btn.setEnabled(False)
        self.simplfy_step_btn.setEnabled(False)
        self.reset_camera_btn.setEnabled(False)
        self.simplified_label.setText("正在执行右侧网格简化，请稍候...")
        self.simplified_label.repaint()  # 强制刷新界面显示
        V = self.current_v
        F = self.current_f
        # 执行简化
        SetMesh(V, F)
        Simplify(self.simplify_ratio)
        V_simplified, F_simplified = GetMesh()
        mesh_simplified = trimesh.Trimesh(vertices=V_simplified, faces=F_simplified)
        # 保存简化后的网格
        mesh_simplified.export(self.simplified_mesh_path)
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
    def animate_vertex_move(self, renderer, source_id, target_id, steps=50, delay=0.01):
        # 获取 actor
        # renderer = renderer.GetRenderer()
        actors = renderer.GetActors()
        actors.InitTraversal()
        actor = actors.GetNextActor()
        if not actor:
            print("No actor found.")
            return

        polydata = actor.GetMapper().GetInput()
        points = polydata.GetPoints()

        # 当前顶点坐标
        start_pos = np.array(points.GetPoint(source_id))
        target_pos = np.array(points.GetPoint(target_id))
        
        # 生成插值路径
        for i in range(steps + 1):
            t = i / steps
            interp_pos = (1 - t) * start_pos + t * target_pos
            points.SetPoint(source_id, *interp_pos)
            points.Modified()  # 通知 VTK 点已更新
            renderer.GetRenderWindow().Render()
            time.sleep(delay)

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
        
    
    def simplify_single_step(self):
        if not self.display_mode == "step":
            self.simplified_label.setText("无法进行单步调试...")
            self.simplified_label.repaint()
            return
        
        new_v, new_f, active_v = SimplifyStep()
        active_faces, only_v1, only_v2 = get_active_faces(self.current_f, active_v)
        mask = np.any(new_f == active_v[2], axis=1)
        only_v1_new = np.where(mask)[0]

        # 为 active faces 上色：face_id -> (r,g,b)
        color_dict_left = {}
        color_dict_left.update({fid: (255, 128, 128) for fid in active_faces})
        color_dict_left.update({fid: (255, 192, 128) for fid in only_v1})
        color_dict_left.update({fid: (255, 192, 128) for fid in only_v2})
        
        color_dict_right = {}
        color_dict_right.update({fid: (255, 255, 128) for fid in only_v1_new})

        # 渲染左侧（当前 mesh，染色）
        self.render_mesh(self.renderer_left, self.current_v, self.current_f, (0.8, 0.9, 0.8), face_colors=color_dict_left)

        # 渲染右侧（简化后 mesh，无染色）
        self.render_mesh(self.renderer_right, self.current_v, self.current_f, (0.8, 0.8, 0.9), face_colors=color_dict_left)
        
        self.current_v = new_v
        self.current_f = new_f

        # 获取右侧的 renderer（也可以选择左侧）
        renderer = self.renderer_left
        camera = renderer.GetActiveCamera()
        actors = renderer.GetActors()
        actors.InitTraversal()
        actor = actors.GetNextActor()
        if not actor:
            return
        polydata = actor.GetMapper().GetInput()
        if not polydata:
            return
        face_id = active_faces[0]
        print(len(self.current_v), active_v, active_faces)
        
        cell = polydata.GetCell(face_id)
        p0 = np.array(polydata.GetPoint(cell.GetPointId(0)))
        p1 = np.array(polydata.GetPoint(cell.GetPointId(1)))
        p2 = np.array(polydata.GetPoint(cell.GetPointId(2)))
        face_center = (p0 + p1 + p2) / 3.0
        # 获取法向量（用于相机方向）
        normal = np.cross(p1 - p0, p2 - p0)
        normal = normal / np.linalg.norm(normal)
        # 设置相机的位置：在法线方向上方离面中心一定距离
        distance = 0.5 * np.linalg.norm(polydata.GetLength())  # 可调
        camera_position = face_center + normal * distance


        # 设置相机参数
        def on_camera_done():
            camera.SetViewUp(0, 1, 0)
            renderer.ResetCameraClippingRange()
            self.renderWindow.Render()
            time.sleep(1)
            self.animate_vertex_move(self.renderer_right, active_v[1], active_v[0])
            self.simplified_label.setText(f"简化完成，当前顶点数: {len(self.current_v)}, 面数: {len(self.current_f)}\n")

        self.animate_camera_to(camera_position, face_center, on_finish=on_camera_done)
        
    def simplfy_step(self):
        self.display_mode = "step"
        self.load_right_btn.setEnabled(False)
        self.reset_camera_btn.setEnabled(True)
        V = self.current_v
        F = self.current_f
        # 执行简化
        self.simplified_label.setText("已切换为单步简化查看模式")
        SetMesh(V, F)

if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    lib_suffix = ".dll"
    window = VTKWindow(lib_suffix=lib_suffix)
    sys.exit(app.exec_())
