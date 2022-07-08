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
        A.dot(B) == ||A|| * ||B|| * cosθ;
        cosθ == A.dot(B) / (||A|| * ||B||) == A.normalized.dot(B.normalized);
        A.dot(B) == B.dot(A);
        A.dot(B + C) == A.dot(B) + A.cot(C);
        (kA).dot(B) == A.dot(kB) = k * A.dot(B);

        #In 2D(二维计算)
        A.dot(B) == xa * xb + ya * yb;

        #In 3D(三维计算)
        A.dot(B) == xa * xb + ya * yb + za * zb;

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

        #Orthonormal
        //Any set of 3 vectors (in 3D) that
        Eigen::Vector3f U, V, W, P;
        U.dot(V) == V.dot(W) == U.dot(W) == 0;
        W == U.dot(V);
        P == P.dot(U).cross(U) + P.dot(V).cross(V) + P.dot(W).cross(W);
    }
    Matrix(矩阵部分)
    {
        #Definition(定义)
        //一个 m * n 的矩阵；
        float m, n;
        Eigen::MatrixXf N(m, n);
        
        #Multiplication(乘法)
        float M1, M2, M3, M4, M5, M6;
        float N1, N2, N3, N4, N5, N6 ,N7, N8;
        Eigen::MatrixXf M(2, 3), N(4, 2), O;
        M << M1, M2,
             M3, M4,
             M5, M6;
        N << N1, N2, N3, N4,
             N5, N6, N7, N8;
        O << M1 * N1 + M2 * N5, M1 * N2 + M2 * N6, M1 * N3 + M2 * N7, M1 * N4 + M2 * N8,
             M3 * N1 + M4 * N5, M3 * N2 + M4 * N6, M3 * N3 + M4 * N7, M3 * N4 + M4 * N8,
             M5 * N1 + M6 * N5, M5 * N2 + M6 * N6, M5 * N3 + M6 * N7, M5 * N4 + M6 * N8;
        O == M * N; 
        Eigen::MatrixXf A, B, C;
        A * B != B * A;
        (A * B) * C = A * (B * C);
        A * (B + C) = A * B + A * C;
        (A + B) * C = A * C + B * C;
        Eigen::Matrix2f i, m;
        i << -1, 0,
             0, 1;
        m << x,
             y;
        n << -x;
             y;
        n == i * m;

        #Transpose(转置)
        Eigen::MatrixXf A, B;
        A << 1, 2,
             3, 4,
             5, 6;
        B << 1, 3, 5,
             2, 4, 6;
        A.transpose == B;
        (A * B).transpose == B.transpose * A.transpose;
        //I的概念：从左上到右下对角线为1其他位置为0的矩阵；
        A * A.inverse == A.inverse * A == I;
        (A * B).inverse == B.inverse * A.inverse;

        #Vector Multiplication(矢量乘法)
        Eigen::MatrixXf a, b, A;
        a.dot(b) == a.transpose * b;
        a.cross(b) = A * b;
    }
*/