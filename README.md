零偏斜时钟树构建的DME算法
---
### TODO
1. 由于本文档涉及很多专业概念，请老师务必仔细阅读本 readme 文档
2. 配置网表生成程序运行环境: 配置方法见本文档“运行环境”部分
3. 打开 main.c, 修改 benchmark_file_path 的绝对路径 (将 path-to-ZST-DME 修改为老师存储本作业代码的路径) 以及 netlist_file_path 的绝对路径(同理)
4. 运行 main.c (推荐使用CLion, 注意配置好CMakeLists.txt)，生成网表文件
5. 配置 embedding 程序运行环境：配置方法见本文档“运行环境”部分(建议复制draw.cpp的代码到这一步新建的项目文件中, 作为 embedding 程序)
6. 修改 embedding 程序的 main 函数中 netlist.txt 的绝对路径，同上
7. 运行 embedding 程序 (这一步无法使用CLion, 推荐按照运行环境的要求使用 Visual Studio, 注意配置好项目属性和 graphics.h), 画出布线图像。

### 问题的提出
在二维曼哈顿平面上给定一组引脚点(sink)(s0, s1, s2,...)，和一个信号源点，仅通过横线和竖线将这些引脚点与信号源点连通，但是需要满足两个条件：
- 从信号源点到每一个引脚点的线长都相等
- 这个线长要尽可能的小

之所以这么做是因为：
- 根据线性延时模型( Linear Delay Model )，线长线性地决定了延时大小
- 从信号源点到每一个引脚点的线长都相等可以保证芯片中各个功能块的运行是同步的
- 这个线长要尽可能的小：保证芯片的运行周期尽可能的小(频率或者说运行速度尽可能的快)

### 主要概念
- sink: 引脚点
- 合并：点A和点B合并到点C，即为把点C作为点A和点B的父节点，点A和点B作为点C的左右子节点
- merging point: 如果我们将点A和点B合并到点C, 那么点C就是一个merging point, 根据算法步骤中的介绍，这个merging point 实际上就是其两个子节点的中点
- Triple: 点C的 Triple 表示点C及其左右子节点
- 曼哈顿距离：[参考链接: 百度百科-曼哈顿距离](https://baike.baidu.com/item/%E6%9B%BC%E5%93%88%E9%A1%BF%E8%B7%9D%E7%A6%BB/743092?fr=aladdin)
- 曼哈顿平面：使用曼哈顿距离度量两点间距离的平面
- 时钟树：节点包含坐标信息的满二叉树
- 网表文件：存储时钟树所有节点的 Triple 的文件
- embedding: 根据网表文件进行实际布线的工艺，本大作业通过画图进行模拟
- 其他概念：如果以上没有您需要查找的概念，那么可以在算法步骤或者其他部分查找

### 算法步骤
这个算法称为 DME 算法，详见附件 DME.pdf 的介绍
1. 网表生成程序(主函数在 main.c )<br> 
生成二叉树并存储为文件(网表文件)，这个二叉树的叶子节点是引脚点(s0, s1, s2...), 其他节点是起辅助作用的点(merging segment), 这里取 merging segment 的中点即 merging point 代替 merging segment.<br>
(1) 读取 benchmark 文件(文件中是所有引脚点的坐标)，将这些坐标存储在一个数组中，称为递归数组;<br>
(2) 另外定义一个路径数组(元素是有序的)，用于之后构建二叉树;<br>
(3) 计算递归数组所有引脚点对中曼哈顿距离最小的点对;<br>
(4) 计算这个点对的中点;(可以证明这个中点也是 merging segment 的中点，即 merging point)<br>
(5) 将点对和中点(一个Triple的**值**)按顺序添加到路径数组当中;<br>
(6) 将这个点对从递归数组中删去，再将它们的中点添加到递归数组当中;<br>
(7) 针对新的递归数组，重复(3)-(7)步骤，直到递归数组中只剩下一个点为止(即进行(8));<br>
(8) 将路径数组中的所有点都变换为树节点结构体数据类型;<br>
(9) 按照路径数组的顺序，将第1个点第2个点作为第3个点的左右子节点，构建Triple；将第4个点和第5个点作为第6个点的左右子节点，构建Triple......依次类推，直到整个路径数组遍历完成;<br>
(10) 返回最后的子树的根节点 root，这个根节点便代表了整颗二叉树。<br>
(11)从树的根部开始，自顶向下地对树进行**中序遍历**，将遍历到的节点的值(即横纵坐标)以Triple为单位写到一个文件中，这个文件类似于实际芯片布线中的网表文件，布线机器通过读取该文件计算布线时的走线方向和位置等;<br>
2. embedding 程序(主函数在 embedding/draw.cpp) <br>
连线画图，用以模拟在实际芯片时钟树中的布线<br>
(1) 按照网表文件的顺序，将第1个点和第2个点连接起来，将第2个点和第3个点连接起来(一个Triple)，然后将第4个点和第5个点连接起来，将第5个点和第6个点连接起来(注意顺序，每三个点Triple为一组，即"左子节点-父节点-右子节点"为一个操作单位)......依次类推，直到文件读取操作完毕;<br>
(2) 结束画图，这个图像便是最终芯片时钟树布线的示意。<br>
*注意，在连接两点时，不能用直线相连，而只能用横纵折线相连 (例如: 点A(2, 7)和点B(3, 4)相连时候，需要按照(2, 7)-(3, 7)-(3, 4)进行连接)<br>
*对于信号源点，实际芯片布线中，它是在曼哈顿平面上与 root 节点相连的，这个在代码中没有也不必体现

### 主要文件说明
- main.c 网表生成主控 
- benchmark.txt 输入的 benchmark 文件
- reader.h benchmark 文件读取解析
- util.h 二叉树基本结构和函数定义
- topoparser.h 二叉树生成部分的代码
- /embedding embedding 程序，需要更换运行环境; draw.cpp 为程序入口，netlist.txt 为运行网表生成程序后生成的网表文件
- netlist.h 网表文件生成
- DME.pdf 本大作业所使用的算法的介绍
- CMakeLists.txt 项目部署配置文件(CLion)


### 架构设计
- 上述已经基本说明
- 简单讲就是: <br>
``
benchmark (输入) + util.h(基础) -> reader.h (读取) -> topoparser.h (树生成) -> netlist.h (网表生成) ->结束网表生成程序
``      
``
draw.cpp (画图) ->结束 embedding 程序
``

### 运行环境
- 网表生成程序: 使用 CLion 作为 IDE, 语言规范为 C++ 11, 部署工具为 CMake ( 在 CMakeLists.txt 里写明了部署所需信息，特别注意下面这一行：)  
``add_executable(ZST_DME main.c)``
- embedding 程序：使用 Visual Studio 作为 IDE (含 Visual C++ 2015-2019), 之后在 [EasyX 链接](http://www.easyx.cn/downloads/)上下载安装 graphics.h 头文件(只能在 Visual C++ 环境下使用)，打开 Visual Studio 新建项目"C++ 控制台应用程序" (注意是C++不是C, 即必须是.cpp文件)，打开“项目—属性—C/C++—预处理器—预处理器定义”， 删除之前内容，改为以下内容:<br>
``_CRT_SECURE_NO_DEPRECATE``  
``_CRT_NONSTDC_NO_DEPRECATE``  
然后建议将 embedding/draw.cpp文件的代码复制到这个新建的项目文件中运行，否则，直接运行 draw.cpp 不能保证运行环境是完全配套的。

### 实验结果
![embedding picture](http://a1.qpic.cn/psc?/05f296b7-f920-4499-af25-b1090ac6d0d1/4KbNA3H1osI2VUAtoM9GOuAux2Vp4zBPd.qyPSH3*NMcExc1n.Ag0Yc6bD*2GbKS.AD8LkynUKcx*F20L*CGbQ!!/c&ek=1&kp=1&pt=0&bo=vgO4A74DuAMRADc!&tl=1&tm=1593594000&sce=0-12-12&rf=0-18)
这个是画出的时钟树布线图像。

### 声明
- 关于需要手动设置两个文件的路径：因为我们发现不同的 IDE 对于相对路径的处理不同，例如：有的 IDE 把./xxx 视为与 .exe 文件同级的路径，有的 IDE 把./xxx 视为与 main.c 文件同级的路径；此外，如果是与.exe文件同级的情况，不同的 IDE 生成的.exe文件的路径又不同。因此推荐老师手动设置绝对路径，以避免不必要的麻烦。
- 关于语言标准: 本大作业主要采用 C++ 11, 如果运行遇到一些问题，可能是是项目运行环境中的语言标准不一(例如我们发现 C++ 11 可以将 char* 赋给 const char* 类型的参数，但是其他语言标准不可以)
- 本大作业由中法工程师学院2018级学生王泽坤、高翊恒、杨冠宇、陈安琪、芮一鸿独立完成
- 本大作业参考资料仅限芯片时钟树布线的算法基础介绍和曼哈顿距离的数学基础介绍，编程部分只参考了 graphics.h 的 API 文档
- 为了提高大作业的小组合作开发效率，代码在[小组成员王泽坤的 Github](https://github.com/ZenMoore/ZST-DME) 上进行了开源，不属于抄袭
- 本大作业对 DME.pdf 中介绍的算法进行了一点简化，即取 merging segment 的中点即 merging point 代替 merging segment, 在这个变化的基础上，无需 TRR 等的计算，详见“算法步骤”部分。