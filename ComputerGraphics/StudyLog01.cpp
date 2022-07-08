//  #include <*/Eigen>
//  #include <optional>
//  #include <algorithm>
//  #include <math.h>

//  向量及线性代数部分复习
/*  Vectors(向量部分)
    {
        #Definition(定义)
        Eigen::Vector2f A, B;
        A << 0, 1;
        B << 2, 3;
        Eigen::Vector2f AB;
        AB << B - A;

        #Normalization(归一化)
        AB.normalized == AB / ||AB||.

        #Addition(加法)
        float x, y;
        Eigen::MatrixXf A(2, 1), AT(1, 2);
        A << x,
             y;
        AT << x, y;
        ||A|| == sqrt(pow(x, 2) + pow(y, 2));

        #Multiplication(乘法)

        #Dot Product(点乘)
        float k, xa, xb, xc, ya, yb, yc, za, zb, zc;
        Eigen::Vector3f A, B, C, Bp;
        A · B == ||A|| * ||B|| * cosθ;
        cosθ == A · B / (||A|| * ||B||) == A.normalized · B.normalized;
        A · B == B · A;
        A · (B + C) == (A · B) + (A · C);
        (kA) · B == A · (kB) = k (A · B);

        #In 2D(二维计算)
        A · B == xa * xb + ya * yb;

        #In 3D(三维计算)
        A · B == xa * xb + ya * yb + za * zb;

        #Projection(投影)
        //  B投影到A上
        Bp == k * A.normalized;
        k == ||B|| * cosθ;
        0° <= θ <= 90°;

        #Cross Product(叉乘)
        float k, Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz;
        Eigen::Vector3f X, Y, Z;
        //空间内任意一组坐标轴可两两叉乘得到第三条轴
        X.cross(Y) == +Z;
        X.cross(Y) == -Y.cross(X);
        X.cross(X) == 0(0向量);
        X.corss((Y + Z)) == X.cross(Y) + X.cross(Z);
        X.cross((k * Y)) == k * (X.cross(Y));
        X << Xx,
             Yx,
             Zx;
        Y << Xy,
             Yy,
             Zy;
        Z << Xz,
             Yz,
             Zz;
        X.cross(Y) == Yx * Zy - Yy * Zx;
                 Zx * Xy - Xx * Zy;
                 Xx * Yy - Yx * Xy;
        Eigen::MatrixXf A(3, 3), b(3, 1);
        A << 0, -Zx, Yx,
             Zx, 0, -Xx,
             -Yx, Xx, 0;
        b << Xy,
             Yy,
             Zy;
        X.cross(Y) == A * b;
    }
