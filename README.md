# SurfaceSimplification

## 项目结构说明

- 前端可视化为`visualization.py`
- 后端代码均位于 `lib` 目录下。
- 前后端接口分别在 `interface.py` 和 `main.cpp` 中实现。
- 依赖库 `glm`（G-Lib Math）需要从以下[仓库](https://github.com/g-truc/glm.git)下载
- 依赖库 `eigen` 需要从以下[仓库](https://github.com/PX4/eigen.git)下载

下载后将其放置于 `lib` 目录下。

## 安装方法

- 运行项目前，请执行 `compile.sh` 脚本进行编译
- Windows系统则使用 `compile_win.bat` 进行编译

python 依赖

    numpy==2.2.1
    vtk==9.4.2
    trimesh==4.6.10
    PyQt5

直接使用pip安装即可

## 运行方法

![可视化界面](display_images/ui.png "可视化界面")

安装好python依赖后，在命令行中输入：

    python visualization.py


- 操作过程可见于 [视频](https://cloud.tsinghua.edu.cn/f/828e5ce1eb1c4fdd8cfa/)
- 一些测试样例存在 `/obj` 文件夹下，在搜索框中输入其名称即可进行化简，不需要带文件后缀
- 需要查看单步简化，首先点击重置按钮，并点击单步简化按钮，之后可查看单步简化
- 若想要使用自己寻找的网格，则可以用fix.py进行水密化处理，否则就会见到非水密网格所取得的抽象结果

## 部分结果展示

本次实验共在Objaverse数据集上选取了550个物体进行了量化实验，下图是部分可视化的结果，实验物体包括斯坦福公开3D模型和Objaverse中的数据

![结果展示](display_images/result_00.png "结果展示")

量化结果展示如下

| **Metrics**      | **QEM**   | **Lindstrom-Turk** | **Spectral** | **Structure-Aware** |
|------------------|-----------|--------------------|--------------|----------------------|
| CD ↓             | 0.0160    | **0.0042**         | 0.0129       | 0.0094               |
| EMD ↓            | 0.0347    | **0.0226**         | 0.0314       | 0.0290               |
| SSIM ↑           | 0.9477    | **0.9674**         | 0.9319       | 0.9579               |
