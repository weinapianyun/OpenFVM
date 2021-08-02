/***************************************************************************
 *   Copyright (C) 2004-2008 by OpenFVM team                               *
 *   http://sourceforge.net/projects/openfvm/                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <time.h>

#include "ttime.h"

#include "globals.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "variables.h"
#include "gradient.h"
#include "geocalc.h"
#include "setup.h"
#include "ttime.h"

#include "msolver.h"

void
PrintVector (Vector * x) // 在命令行中打印向量X
{
    int i;
    int nr;
    nr = V_GetDim (x);

    for (i = 0; i < nr; i++)
    {
        printf ("%.16e\n", V_GetCmp (x, i + 1));
    }
}

void
WriteVector (Vector * x) // 将向量X 写入文件vector.m中
{
    FILE *fp; // 文件指针
    int i;
    int nr;

    nr = V_GetDim (x); // 获取向量 x 的规格
    fp = fopen ("vector.m", "w"); // 创建并以只写格式打开 vector.m文件

    fprintf (fp, "Vec_s = [\n"); // 将 "vec_s ="输入到 fp(vector.m)中

    for (i = 0; i < nr; i++)
    {
        fprintf (fp, "%.16e\n", V_GetCmp (x, i + 1)); // 得到 X的值并写入fp中
    }

    fprintf (fp, "];\n"); // 以 “];\n”作为文件vector.m的结束
    fclose (fp); // 关闭文件
}

void
WriteMatrix (QMatrix * A, int symmetric) // 将矩阵A 写入文件Matrix.m中
{
    FILE *fp;

    int i, j;
    int nr; // 存储行数
    int nj; // 存储每一行的元素个数
    int nonzeros; // 存储 矩阵全部元素的个数

    nr = Q_GetDim (A); // 矩阵的维度(行数)
    fp = fopen ("matrix.m", "w"); // 创建并以只写格式打开 matrix.m文件

    nonzeros = 0;
    for (i = 0; i < nr; i++)
    {
        nj = Q_GetLen (A, i + 1); // 获得矩阵A每一行的长度
        nonzeros += nj; // 矩阵全部元素的个数
    }

    fprintf (fp, "%% Size = %d %d\n", nr, nr);
    fprintf (fp, "%% Nonzeros = %d\n", nonzeros);
    fprintf (fp, "zzz = zeros(%d,3);\n", nonzeros);
    fprintf (fp, "zzz = [\n");

    for (i = 0; i < nr; i++)
    {
        nj = Q_GetLen (A, i + 1); // 获取矩阵A行的长度
        for (j = 0; j < nj; j++)
        {
            fprintf (fp, "%d %zu  %.16e\n", i + 1, Q_GetPos (A, i + 1, j),
                     Q_GetVal (A, i + 1, j)); // 打印 序号、位置、元素值
        }
    }

    fprintf (fp, "];\n"); // 打印 "["并换行
    fprintf (fp, "Mat_s = spconvert(zzz);\n"); // 打印 "Mat_s = zzz矩阵转换形式"

    if (symmetric == LOGICAL_TRUE) // 若矩阵A 为对称矩阵
        fprintf (fp, "Mat_s = triu(Mat_s, 1) + tril(Mat_s', 0);\n");
}

int
CheckIfDiagonalMatrix (QMatrix * A) // 检查矩阵A是否对角占优
{ /*  严格对角占优矩阵： 每一行中对角元素的值的模 > 其余元素值的模之和
   *    弱对角占优矩阵： 每一行中对角元素的值的模 >= 其余元素值的模之和  */

    int i, j;
    int nr; // 矩阵行数
    int nj; // 矩阵每一行元素数
    int cond1, cond2; // 两个占优条件 判断参数

    double sum_ap_parcial, sum_an_parcial; // 存储每一行的对角元素、其他元素之和
    double sum_ap_total, sum_an_total; // 存储整个矩阵的对角元素之和、其他元素之和

    sum_ap_total = 0.0;
    sum_an_total = 0.0;

    cond1 = LOGICAL_FALSE;
    cond2 = LOGICAL_FALSE;
    nr = Q_GetDim (A); // 获取矩阵A的行数

    sum_ap_parcial = 0.0;

    for (i = 0; i < nr; i++) // 遍历矩阵的每一行
    {
        nj = Q_GetLen (A, i + 1); // 获取每一行的元素数
        if (nj > 0)
            sum_ap_parcial = LABS (Q_GetVal (A, i + 1, i)); // 矩阵A 第 i+1 行第 i+1 个元素值 的绝对值

        sum_an_parcial = 0.0;

        for (j = 0; j < nj; j++)
        {
            if(j != i) // 非对角线元素
                sum_an_parcial += LABS (Q_GetVal (A, i + 1, j)); // 计算矩阵A i+1行的其他元素值的绝对值之和
        }

        if (sum_ap_parcial >= sum_an_parcial) // 通过比较代数值来判断 对角元素 的占优性
            cond1 = LOGICAL_TRUE;
        sum_ap_total += sum_ap_parcial; // 对角元素绝对值之和
        sum_an_total += sum_an_parcial; // 其他元素绝对值之和
    }

    if (sum_ap_total >= sum_an_total) // 通过比较对角元素之和 与 其他元素之和 来判断占优性
        cond2 = LOGICAL_TRUE;

    if (!cond1 && !cond2) // 当两个占优条件有一个不满足，矩阵A 均不满足对角占优
        return LOGICAL_FALSE;
    else
        return LOGICAL_TRUE;
}

void
SolveMatrix (QMatrix * A, Vector * x, Vector * b, int *iter, double *res,
             double *time, int msolver, int mprecond, int miter) // 代数矩阵迭代求解
{
    PrecondProcType PrecondProc; // PrecondProc 为定义的函数指针
    double start, end; // 存储起止时间
    start = ttime (); // 求解开始时刻

    // (0-Null, 1-Jacobi, 2-SOR, 3-ILU) 预处理方法

    PrecondProc = 0; // 函数指针初始化
    switch (mprecond) // 判断设定的预处理方法
    {
        case pNONE:
            PrecondProc = NULL; // 无预处理
            break;

        case pJACOBI:
            PrecondProc = JacobiPrecond; // Jacobi预处理
            break;

        case pSOR:
            PrecondProc = SSORPrecond; // 超松弛预处理
            break;

        case pILU:
            PrecondProc = ILUPrecond; // 不完全LU分解预处理
            break;

        case pASM:
            PrecondProc = ILUPrecond; // 不完全LU分解预处理
            break;
    }

    // (0-Jacobi, 1-SOR, 2-QMR, 3-GMRES, 4-CG, 5-CGN, 6-CGS, 7-BiCG, 8-BiCGStab)

    switch (msolver) // 判断设定的迭代方法
    {
        case sJACOBI:
            // Jacobi
            JacobiIter (A, x, b, miter, PrecondProc, 1.0); // 使用 jacibo迭代法求解矩阵

            *iter = GetLastNoIter (); // 获取最后一次迭代后的总迭代次数
            *res = GetLastAccuracy (); // 获取在迭代最后达到的精度
            break;

        case sSOR:
            // SOR 超松弛迭代法
            SSORIter (A, x, b, miter, PrecondProc, 1.0);

            *iter = GetLastNoIter ();
            *res = GetLastAccuracy ();
            break;

        case sQMR:
            // QMR 拟最小残差法
            QMRIter (A, x, b, miter, PrecondProc, 1.0);

            *iter = GetLastNoIter ();
            *res = GetLastAccuracy ();
            break;

        case sGMRES:
            // GMRES 广义最小残差法
            GMRESIter (A, x, b, miter, PrecondProc, 1.0);

            *iter = GetLastNoIter ();
            *res = GetLastAccuracy ();
            break;

        case sCG:
            // CG 共轭梯度法
            CGIter (A, x, b, miter, PrecondProc, 1.0);

            *iter = GetLastNoIter ();
            *res = GetLastAccuracy ();
            break;

        case sCGN:
            // CGN  normal(自然)方程的共轭梯度法
            CGNIter (A, x, b, miter, PrecondProc, 1.0);

            *iter = GetLastNoIter ();
            *res = GetLastAccuracy ();
            break;

        case sCGS:
            // CGS  平方共轭梯度法
            CGSIter (A, x, b, miter, PrecondProc, 1.0);

            *iter = GetLastNoIter ();
            *res = GetLastAccuracy ();
            break;

        case sBICG:
            // BICG 双共轭梯度法
            BiCGIter (A, x, b, miter, PrecondProc, 1.0);

            *iter = GetLastNoIter ();
            *res = GetLastAccuracy ();
            break;

        case sBICGS:
            // BICGS 稳定双共轭梯度法
            BiCGSTABIter (A, x, b, miter, PrecondProc, 1.0);

            *iter = GetLastNoIter ();
            *res = GetLastAccuracy ();
            break;

        default: // 默认打印语句 “Solver method (%d) not implemented”
            printf ("\nError: Solver method (%d) not implemented\n", msolver);
            exit (LOGICAL_ERROR);
            break;
    }

    end = ttime (); // 求解结束时刻
    *time = end - start; // 求解时间
}
