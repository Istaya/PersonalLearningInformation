/*  Transformation(变形)
    2D Transformation(2D变形)
    为了将模型空间，视口(摄像机)空间这两种独立的体系统一到相同的世界空间中，需要人为的将模型与摄像机的位置移动到世界中心的原点上。
    首先需要学习的是Modeling(模型)变化,在学习3D空间的变形前，先学习2D空间下模型的变形规则。

    #Scale(比例变化)
    首先考虑的是2D的比例变化问题，对于已经处于世界中心原点的物体，可以通过直接对该模型各点坐标的X，Y数值乘上一个系数S来进行大小的缩放。
    X' = S * X;
    Y' = S = Y;
    写成矩阵的形式。
    Eigen::MatrixXf A, A', S;
    A << X,
         Y;
    A' << X',
          Y';
    S << s, 0,
         0, s;
    A' = S * A;
    当X轴与Y轴缩放比例不相同时，只需要将S矩阵中s值进行替换就可以实现需要的效果。
    S << sx, 0,
         0, sy;
    A' = S * A;
    关于X轴或Y轴对称只需要保证一条轴上的数值不变，另一条轴上的数值取反；中心对称则两轴数值都取反即可。
    Eigen::MatrixXf Sx, Sy, So;
    Sx << 1, 0,
          0, -1;
    Sy << -1, 0,
          0, 1;
    So << -1, 0,
          0, -1;
    A' = Sx & Sy & So * A;

    #Shear Matrix(剪切矩阵)
    剪切矩阵可以使模型沿着X轴或Y轴方向倾斜，该过程也被称之为斜切。
    矩阵表示为
    S << 1, ax,
         ay, 1;
    A' = S * A;
    其中ax为在x轴方向上斜切的距离，ay为在y轴方向上斜切的距离。