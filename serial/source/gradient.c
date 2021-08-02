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

#include "mesh.h"
#include "bcond.h"
#include "geocalc.h"
#include "globals.h"
#include "variables.h"

#include "gradient.h"

msh_vector
Gradient(Vector *phi, Vector *phif, int bound, int element) // 梯度计算的 CG 单元基法
{
    // Cell based
    int j;
    int neighbor, face, pair;
    double phij; // 界面面心的变量值

    //double dNf, dPf;
    double lambda; // 插值因子

    double v1, v2;
    double v1min, v1max, v2min, v2max;
    double factor; // 梯度计算的中间变量

    double fv;
    msh_vector rv; // 存储梯度变量值

    rv.x = 0.0; // 梯度变量初始化
    rv.y = 0.0;
    rv.z = 0.0;

    for (j = 0; j < elements[element].nbfaces; j++) // 遍历单元的所有界面
    {
        face = elements[element].face[j]; // 单元界面编号
        pair = faces[face].pair;          // 共面界面的编号

        if (pair != -1) // 非边界单元
        {
            neighbor = faces[pair].element; // 获取界面相邻单元的编号
            /*
            dNf =
              GeoMagVector (GeoSubVectorVector
                    (elements[neighbor].celement, faces[face].cface));
            dPf =
              GeoMagVector (GeoSubVectorVector
                    (elements[element].celement, faces[face].cface));
            lambda = dPf / (dPf + dNf); // 插值因子的计算
            */
            lambda = 0.5;

            // Element face variable , 插值算法 Uf = UNf λ + UP (1 − λ)
            phij = V_GetCmp(phi, neighbor + 1) * lambda +
                   V_GetCmp(phi, element + 1) * (1.0 - lambda);

            factor = phij / elements[element].Vp; // 中间变量 Uf / Vp

            // Element center gradient , grad(U) = 求和（Uf * Af ）/ Vp
            rv.x += factor * faces[face].A.x;
            rv.y += factor * faces[face].A.y;
            rv.z += factor * faces[face].A.z;
        }
        else
        {
            // Element face variable
            if (bound == LOGICAL_TRUE)
                phij = V_GetCmp(phif, face + 1); // 边界界面直接取边界条件的值
            else
                phij = V_GetCmp(phi, element + 1); // 直接取网格单元中心的值作为面上的值

            factor = phij / elements[element].Vp;

            // Element center gradient 网格单元中心的梯度值
            rv.x += factor * faces[face].A.x;
            rv.y += factor * faces[face].A.y;
            rv.z += factor * faces[face].A.z;
        }
    }

    // Eliminate local extrema , 消除局部极值
    v1min = +GREAT;
    v1max = -GREAT;

    v2min = +GREAT;
    v2max = -GREAT;

    for (j = 0; j < elements[element].nbfaces; j++)
    {
        face = elements[element].face[j];
        v1 = V_GetCmp(phi, element + 1) + rv.x * (faces[face].cface.x - elements[element].celement.x) + rv.y * (faces[face].cface.y - elements[element].celement.y) + rv.z * (faces[face].cface.z - elements[element].celement.z);
        v1min = LMIN(v1min, v1);
        v1max = LMAX(v1max, v1);

        pair = faces[face].pair;

        if (pair != -1)
        {
            neighbor = faces[pair].element;

            v2 = (V_GetCmp(phi, element + 1) + V_GetCmp(phi, neighbor + 1)) * 0.5;
            v2min = LMIN(v2min, v2);
            v2max = LMAX(v2max, v2);
        }
    }

    fv = 1.0;

    if (LABS(v1max - v1min) > 2 * LABS(v2max - v2min))
        fv = LABS(v2max - v2min) / LABS(v1max - v1min);

    rv.x *= fv;
    rv.y *= fv;
    rv.z *= fv;

    return rv; // 单元中心梯度
}

msh_vector
GradientN(Vector *phin, Vector *phif, int bound, int element) // 梯度计算的 CG 顶点基法
{
    // Node based
    int j, k;
    int node, face, pair; // 节点编号，界面编号，共面界面编号
    double phij;          // 界面面心的变量值
    double factor;        // 中间变量

    msh_vector rv; // 存储变量值

    rv.x = 0.0;
    rv.y = 0.0;
    rv.z = 0.0;

    for (j = 0; j < elements[element].nbfaces; j++) // 遍历单元的所有界面
    {
        face = elements[element].face[j];
        pair = faces[face].pair;

        if (pair != -1) // 非边界界面
        {
            phij = 0.0;
            for (k = 0; k < faces[face].nbnodes; k++) // 遍历界面包含的节点
            {
                node = faces[face].node[k];
                phij += V_GetCmp(phin, node + 1); // phin为节点的变量值
            }

            // Element face variable
            if (faces[face].nbnodes > 0)
                phij /= faces[face].nbnodes; // Uf = 求和(Ufn) / Nf
            factor = phij / elements[element].Vp;

            // Element center gradient
            rv.x += factor * faces[face].A.x;
            rv.y += factor * faces[face].A.y;
            rv.z += factor * faces[face].A.z;
        }
        else
        {
            // Element face variable
            if (bound == LOGICAL_TRUE) // 边界界面直接取边界条件的值
            {
                phij = V_GetCmp(phif, face + 1);
            }
            else // 由节点变量值计算界面变量值
            {
                phij = 0.0;
                for (k = 0; k < faces[face].nbnodes; k++)
                {
                    node = faces[face].node[k];
                    phij += V_GetCmp(phin, node + 1);
                }

                // Element face variable
                if (faces[face].nbnodes > 0)
                    phij /= faces[face].nbnodes;
            }
            factor = phij / elements[element].Vp;

            // Element center gradient
            rv.x += factor * faces[face].A.x;
            rv.y += factor * faces[face].A.y;
            rv.z += factor * faces[face].A.z;
        }
    }
    return rv;
}


msh_vector
GradientX(Vector *phi, Vector *phif, int element, double wf) // 最小二乘法 梯度计算方法
{
    int i, j, k;              // 计数
    int neighbor, face, pair; // 节点编号，界面编号，共面的编号
    double df[3] = {0.0};     // 界面两侧单元中心的位移矢量

    msh_vector rv; // 变量的梯度值

    rv.x = 0.0; // 初始化
    rv.y = 0.0;
    rv.z = 0.0;

    // 创建 最小二乘法(Gg = h) , 所需的系数矩阵G , 右端向量h , 梯度解向量g
    V_Constr(&g, "Gradient", 3, Normal, True);
    V_Constr(&h, "Right element", 3, Normal, True);
    Q_Constr(&G, "Left element", 3, False, Rowws, Normal, True);

    V_SetAllCmp(&g, 0.0); // 将向量 g、h 初始化为 0.0
    V_SetAllCmp(&h, 0.0);

    for (j = 0; j < 3; j++) // 将矩阵 G 初始化为 0.0
    {
        Q_SetLen(&G, j + 1, 3);
        for (k = 0; k < 3; k++)
            Q_SetEntry(&G, j + 1, k, k + 1, 0.0);
    }

    for (i = 0; i < elements[element].nbfaces; i++) // 遍历单元的所有界面
    {
        face = elements[element].face[i]; // 获取当前界面编号
        pair = faces[face].pair;          // 获取与当前界面共面的界面编号

        // 用数组 df[3] 存储界面两侧单元 P、N 辅助节点间的位矢
        df[0] = faces[face].d.x;
        df[1] = faces[face].d.y;
        df[2] = faces[face].d.z;

        if (pair != -1) // 非边界界面
        {
            neighbor = faces[pair].element; // 获取界面相邻单元的编号

            for (j = 0; j < 3; j++) // 给右端向量h 系数矩阵 G 赋值
            {
                // h(i) = wf*wf*(phi_P-phi_N)*df(i), 其中 i 为 x y z
                V_SetCmp(&h, j + 1, V_GetCmp(&h, j + 1) + wf * wf * (V_GetCmp(phi, element + 1) - V_GetCmp(phi, neighbor + 1)) * df[j]);
                for (k = 0; k < 3; k++)
                    // G(i,j) = wf*wf*df(i)*df(j), 其中 i、j 为 x y z
                    Q_SetEntry(&G, j + 1, k, k + 1, Q_GetVal(&G, j + 1, k) + wf * wf * df[j] * df[k]);
            }
        }
        else // 边界界面
        {
            for (j = 0; j < 3; j++) // 边界界面不存在相邻单元N, 故用界面的值 phi_f 代替 phi_N
            {
                // h(i) = wf*wf*(phi_P-phi_N)*df(i), 其中 i 为 x y z
                V_SetCmp(&h, j + 1, V_GetCmp(&h, j + 1) + wf * wf * (V_GetCmp(phi, element + 1) - V_GetCmp(phif, face + 1)) * df[j]);
                for (k = 0; k < 3; k++)
                    // G(i,j) = wf*wf*df(i)*df(j), 其中 i、j 为 x y z
                    Q_SetEntry(&G, j + 1, k, k + 1, Q_GetVal(&G, j + 1, k) + wf * wf * df[j] * df[k]);
            }
        }
    }
    // 求解变量的梯度值
    mc = MulInv_QV(&G, &h); // g = h/G , 即 g = G-1 * h
    Asgn_VV(&g, mc);

    // 将Vector类型的向量 g 转换为msh_vector类型的 rv
    rv.x = V_GetCmp(&g, 1);
    rv.y = V_GetCmp(&g, 2);
    rv.z = V_GetCmp(&g, 3);

    return rv; // 返回 rv
}
