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

#include "globals.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "variables.h"
#include "gradient.h"
#include "geocalc.h"
#include "setup.h"
#include "msolver.h"

#include "temperature.h"

void
CorrectFaceT () // 修正界面的温度值
{
    unsigned int i, j;
    register unsigned int face, pair; // 界面的编号、共面状态
    register unsigned int element, neighbor; // 当前单元、邻接单元的编号

    //double dNf, dPf;
    double lambda; // 界面插值因子
    double Tpl; // 单元P辅助节点 p'的温度值

    msh_vector gradTp; // 单元中心的温度梯度

    for (i = 0; i < nbfaces; i++) // 遍历所有界面
    {
        face = i;
        element = faces[face].element; // 界面所在单元号
        pair = faces[face].pair; // 界面共面状态

        if (pair != -1) // 非边界界面
        {
            neighbor = faces[pair].element; // 共面单元的单元编号

            /*
            dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
            dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

            lambda = dPf / (dPf + dNf);
            */

            lambda = 0.5; // 界面插值因子
            // 插值计算界面上的温度值  Tf = lambda * Tn + (1-lambda) * Tp
            V_SetCmp (&xTf, face + 1, V_GetCmp (&xT, neighbor + 1) * lambda +
                                      V_GetCmp (&xT, element + 1) * (1.0 - lambda));
        }
        else // 边界界面
        {
            if (faces[face].bc == ADIABATICWALL || faces[face].bc == OUTLET // 绝热墙、出口、空条件
                || faces[face].bc == EMPTY)                                     // 在其他类型边界条件中温度可设为定值
            {
                // zero gradient 温度梯度为 0
                Tpl = V_GetCmp (&xT, element + 1); // 给单元P辅助节点 p'的值Tp'赋初值

                if (parameter.orthof != 0.0) // 非正交修正
                {
                    gradTp = Gradient (&xT, &xTf, LOGICAL_TRUE, element); // 计算单元P的温度梯度
                    // 计算辅助节点p'的温度值 Tp' = Tp + grad(Tp)*(rp'-rp)
                    Tpl += parameter.orthof * GeoDotVectorVector (gradTp,
                                                                  GeoSubVectorVector (faces[face].rpl,
                                                                                      elements[element].celement));
                }
                V_SetCmp (&xTf, face + 1, Tpl); // 用 Tp'更新界面的温度值 Tpf
            }
        }
    }
}

void
BuildEnergyMatrix (double dt, double schemefactor) // 组装离散化能量守恒方程的代数矩阵 AX=b
{
    // schemefactor 为方案(延时修正)因子
    int i, j, n; // n为单元的非边界界面数目

    register unsigned int face, pair;
    register unsigned int element, neighbor;

    double aep; // 半离散能量守恒方程的 aep系数
    double aen[MAXFACES]; // 半离散能量守恒方程的 aen系数
    unsigned int ani[MAXFACES]; // 存储每个单元各界面相邻单元的编号

    double bep; // 半离散能量守恒方程的 bep系数

    //double dNf, dPf;
    double lambda; // 界面插值因子
    double xsi; // 计算流场变量时，界面的差分因子

    msh_vector gradTp; // 单元 p的温度梯度
    msh_vector gradTn; // 单元 n的温度梯度

    double densp; // 计算中使用的密度
    double spheatp; // 计算中使用的比热容
    double thcondj; // 计算中使用的热导率

    for (i = 0; i < nbelements; i++) // 遍历网格单元
    {
        element = i;

        aep = 0.0; // 各系数初始化
        bep = 0.0;
        n = 0;

        if (parameter.orthof != 0.0) // 非正交网格下 计算单元中心温度梯度
            gradTp = Gradient (&xT, &xTf, LOGICAL_TRUE, element);

        densp = V_GetCmp (&dens, element + 1); // 获取当前单元的密度值
        spheatp = V_GetCmp (&spheat, element + 1); // 获取当前单元的比热容

        for (j = 0; j < elements[element].nbfaces; j++) // 遍历每个单元的所有界面
        {
            // 获取当前界面的编号 共面界面的编号
            face = elements[element].face[j];
            pair = faces[face].pair;

            if (pair != -1) // 非边界界面
            {
                neighbor = faces[pair].element;
                /*
            dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
            dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

            lambda = dPf / (dPf + dNf);
            */

                lambda = 0.5;
                // 计算界面上的的热导率  kf = (1-lambda) * kp + lamnda * kn
                thcondj = V_GetCmp (&thcond, element + 1) * (1.0 - lambda) +
                          V_GetCmp (&thcond, neighbor + 1) * lambda;

                // Conduction  热传导
                // 若 sf(schemefactor)!=1, 存在延时修正， aep(传导的修正项) = sf * kf * Af / |df| / Vp
                aep += schemefactor * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;
                // aen(传导项) = -k * Af / |df| / Vp
                aen[n] = -schemefactor * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;
                // 若 sf !=1, 存在延时修正，热传导项的延时项 bep_Cond = -(1-sf)*kf*(Tn0 - Tp0)/ |df| * Af/Vp
                bep += -(1.0 - schemefactor) * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj)
                       / elements[element].Vp * V_GetCmp (&xT0, element + 1);
                bep += +(1.0 - schemefactor) * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj)
                       / elements[element].Vp * V_GetCmp (&xT0, neighbor + 1);

                // Convection  热对流
                if (parameter.scheme[iT] == UDS)
                {
                    // UDS  上风差分形式
                    if (V_GetCmp (&uf, face + 1) > 0.0)
                        xsi = 0.0;
                    else
                        xsi = 1.0;
                }
                else
                {
                    //CDS 中心差分
                    xsi = lambda; // 默认为 0.5
                }

                // Convection   schemefactor的处理与热传导项一致
                // aep(对流项) = sf * (1-xsi) * dens * Cp * Af / Vp * Uf
                aep += schemefactor * (1.0 - xsi) * densp * spheatp * V_GetCmp (&uf, face + 1) *
                       faces[face].Aj / elements[element].Vp;
                // aen(对流项) = sf * xsi * dens * Cp * Af / Vp * Uf
                aen[n] += schemefactor * xsi * densp * spheatp * V_GetCmp (&uf, face + 1) *
                          faces[face].Aj / elements[element].Vp;
                ani[n] = neighbor; // 存储相邻单元的编号
                n++; // 统计单元的非边界界面数

                // 对流延时项, bep_Conv = -(1-sf)*[dens*Cp*Af/Vp*Uf*Tp0* ((1-xsi)Tp0 + xsi*Tn0)]
                bep += -(1.0 - schemefactor) * (1.0 - xsi) * densp * spheatp *
                        V_GetCmp (&uf,face + 1) * faces[face].Aj / elements[element].Vp
                        * V_GetCmp (&xT0, element + 1);
                bep += -(1.0 - schemefactor) * xsi * densp * spheatp *
                        V_GetCmp (&uf,face + 1) * faces[face].Aj / elements[element].Vp
                        * V_GetCmp (&xT0, neighbor + 1);

                if (parameter.orthof != 0.0) // 非正交修正
                {
                    gradTn = Gradient (&xT, &xTf, LOGICAL_TRUE, neighbor); // 计算单元N的温度梯度

                    // 扩散项中 温度梯度的非正交修正项  bep_northcor = kf*Aj/Vp*[grad(Tn)*(rn'-rn)-grad(Tp)*(rp'-rp)]/|df|
                    bep += parameter.orthof * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
                           elements[element].Vp * (GeoDotVectorVector (gradTn, GeoSubVectorVector (faces[face].rnl,elements[neighbor].celement))
                           - GeoDotVectorVector (gradTp,GeoSubVectorVector (faces[face].rpl, elements[element].celement)));
                }
            }
            else // 边界界面
            {
                thcondj = material.bthcond; // 边界界面上的热导率

                if (faces[face].bc != EMPTY && faces[face].bc != ADIABATICWALL) // 温度梯度不为 0
                {
                    // Conduction  Tp的系数 , aep_cond = sf*kf*Af/ |df|/Vp
                    aep += schemefactor * thcondj * faces[face].Aj / (faces[face].dj +
                                                                      faces[face].kj) / elements[element].Vp;
                    // 边界界面 扩散项中的Tn使用界面的Tf , bep_cond = sf*kf*Af/Vp * Tf /|df|
                    bep += schemefactor * thcondj * faces[face].Aj /
                            (faces[face].dj +faces[face].kj) / elements[element].Vp * V_GetCmp (&xTf, face + 1);

                    // Convection
                    if (parameter.scheme[iT] == UDS)
                    {
                        // UDS 上风差分形式
                        if (V_GetCmp (&uf, face + 1) > 0.0) // 界面温度值Tf 采用单元中心温度值Tp
                        {
                            // aep_conv = sf*dens*Cp*Uf*Af/Vp
                            aep += schemefactor * densp * spheatp * V_GetCmp (&uf, face + 1) *
                                   faces[face].Aj / elements[element].Vp;
                        }
                        else // 界面温度值Tf 采用已知界面值(初值) , 并将对流项归入源项
                        {
                            // bep_conv = -sf*dens*Cp*Uf*Af/Vp *Tf
                            bep += -schemefactor * densp * spheatp * V_GetCmp (&uf, face +  1) *
                                   faces[face].Aj / elements[element].Vp * V_GetCmp (&xTf, face + 1);
                        }
                    }
                    else
                    {
                        // CDS  直接取边界界面值 Tf, 归入源项  bep_conv = -sf*dens*Cp*Uf*Af/Vp*Tf
                        bep += -schemefactor * densp * spheatp * V_GetCmp (&uf, face + 1) *
                               faces[face].Aj / elements[element].Vp * V_GetCmp (&xTf, face + 1);
                    }

                    if (parameter.orthof != 0.0) // 非正交修正
                    {
                        // 扩散项中 温度梯度的非正交修正项  bep_northcor = kf*Aj/Vp* grad(Tp)* (rp'-rp)/|df|
                        bep += parameter.orthof * thcondj * faces[face].Aj / (faces[face].dj +faces[face].kj)
                                / elements[element].Vp * GeoDotVectorVector (gradTp,
                                                                             GeoSubVectorVector (faces[face].rpl, elements[element].celement));
                    }
                }
            }
        }

        if (dt > 0) // 非定常项的处理
        {
            // Unsteady term  在系数中增加非稳定项
            aep += densp * spheatp / dt; //  + dens*Cp/dt
            bep += densp * spheatp / dt * V_GetCmp (&xT0, element + 1); //  + dens*Cp/dt*Tp0

        }

        if (aep == 0.0 || aep != aep) // 若 app类型改变 或 取值为 0
        {
            printf ("\nError: Problem setting up energy matrix\n");
            exit (LOGICAL_ERROR);
        }

        Q_SetLen (&Ae, element + 1, n + 1); // Ae 矩阵每行的系数数目
        Q_SetEntry (&Ae, element + 1, 0, element + 1, aep); // 将每个单元aep的值放入矩阵Ae的对角位置

        for (j = 0; j < n; j++) // j 为矩阵每行系数元素的数目
        {
            Q_SetEntry (&Ae, element + 1, j + 1, ani[j] + 1, aen[j]); // 将aen[j]的值放入Ae的对应位置
        }
        V_SetCmp (&bT, element + 1, bep); // 设置 bT向量对应位置的值
    }
}

void
SolveEnergyExplicit (double dt) // 显式求解 离散化能量守恒方程
{
    /* 显式求解能量守恒方程 , 与BuildEnergyMatrix函数相比
     * 该函数显式处理半离散化方程中不含待求解变量Tp的其他项 , 最后直接求解Tp */
    unsigned int i, j, n; // 计数

    register unsigned int face, pair;
    register unsigned int element, neighbor;

    double aep; // 半离散能量守恒方程的 aep系数
    double aen[MAXFACES]; // 半离散能量守恒方程的 aen系数
    unsigned int ani[MAXFACES]; // 存储每个单元各界面相邻单元的编号

    double bep; // 半离散能量守恒方程的 bep系数

    //double dNf, dPf;
    double lambda; // 界面插值因子
    double xsi; // 计算流场变量时，界面的差分因子

    msh_vector gradTp; // 单元 p的温度梯度
    msh_vector gradTn; // 单元 n的温度梯度

    double densp; // 计算中使用的密度
    double spheatp; // 计算中使用的比热容
    double thcondj;  // 计算中使用的热导率

    double sumT; // 显式计算各界面相邻单元的影响 , 求和[aeN*Tn0]

    for (i = 0; i < nbelements; i++) // 遍历网格单元
    {
        element = i;

        aep = 0.0; // 各系数初始化
        bep = 0.0;
        n = 0;

        if (parameter.orthof != 0.0) // 非正交网格下 计算单元中心温度梯度
            gradTp = Gradient (&xT, &xTf, LOGICAL_TRUE, element);

        densp = V_GetCmp (&dens, element + 1); // 获取当前单元的密度值
        spheatp = V_GetCmp (&spheat, element + 1); // 获取当前单元的比热容

        for (j = 0; j < elements[element].nbfaces; j++) // 遍历当前单元的所有界面
        {
            face = elements[element].face[j];
            pair = faces[face].pair;

            if (pair != -1) // 非边界界面
            {
                neighbor = faces[pair].element; // 获取当前界面相邻的单元编号

                /*
                dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
                dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

                lambda = dPf / (dPf + dNf);
                */

                lambda = 0.5;
                // 计算界面上的的热导率
                thcondj = V_GetCmp (&thcond, element + 1) * (1.0 - lambda) +
                          V_GetCmp (&thcond, neighbor + 1) * lambda;

                // Conduction 热扩散项系数的计算
                aep += thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;
                aen[n] = -thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

                // Convection  热对流项 界面插值因子的确定
                if (parameter.scheme[iT] == UDS)
                {
                    // UDS
                    if (V_GetCmp (&uf, face + 1) > 0.0)
                        xsi = 0.0;
                    else
                        xsi = 1.0;
                }
                else
                {
                    //CDS
                    xsi = lambda;

                }

                // Convection  热对流项系数的计算
                aep += (1.0 - xsi) * densp * spheatp * V_GetCmp (&uf, face + 1) *
                       faces[face].Aj / elements[element].Vp;
                aen[n] += xsi * densp * spheatp * V_GetCmp (&uf, face + 1) *
                          faces[face].Aj / elements[element].Vp;
                ani[n] = neighbor;
                n++;

                if (parameter.orthof != 0.0) // 热扩散项温度梯度的非正交修正
                {
                    gradTn = Gradient (&xT, &xTf, LOGICAL_TRUE, neighbor);
                    // Non-orthogonal correction term
                    bep += parameter.orthof * thcondj * faces[face].Aj / (faces[face].dj +faces[face].kj) /
                            elements[element].Vp *(GeoDotVectorVector (gradTn,GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement))
                            - GeoDotVectorVector (gradTp, GeoSubVectorVector (faces[face].rpl,elements[element].celement)));
                }

            }
            else // 边界界面
            {
                thcondj = material.bthcond;
                if (faces[face].bc != EMPTY && faces[face].bc != ADIABATICWALL) // 温度梯度不为 0的边界
                {
                    // Conduction  热扩散项系数的计算
                    aep += thcondj * faces[face].Aj / (faces[face].dj +
                                                       faces[face].kj) / elements[element].Vp;
                    bep += thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
                           elements[element].Vp * V_GetCmp (&xTf, face + 1);

                    // Convection  热对流项系数的计算
                    if (parameter.scheme[iT] == UDS)
                    {
                        // UDS
                        if (V_GetCmp (&uf, face + 1) > 0.0)
                        {
                            aep += densp * spheatp * V_GetCmp (&uf, face + 1) *
                                   faces[face].Aj / elements[element].Vp;
                        }
                        else
                        {
                            bep += -densp * spheatp * V_GetCmp (&uf, face + 1) *
                                   faces[face].Aj / elements[element].Vp * V_GetCmp (&xTf, face + 1);
                        }
                    }
                    else
                    {
                        // CDS
                        bep += -densp * spheatp * V_GetCmp (&uf, face + 1) * faces[face].Aj /
                               elements[element].Vp * V_GetCmp (&xTf, face + 1);
                    }

                    if (parameter.orthof != 0.0) // 热扩散项温度梯度的非正交修正
                    {
                        // Non-orthogonal correction term
                        bep += parameter.orthof * thcondj * faces[face].Aj / (faces[face].dj + faces[face].kj) /
                                elements[element].Vp * GeoDotVectorVector (gradTp,GeoSubVectorVector (faces[face].rpl,
                                                                                                      elements[element].celement));
                    }
                }
            }
        }

        if (dt > 0) // 非定常项的处理
        {
            // Unsteady term - Euler
            aep += densp * spheatp / dt;
            bep += densp * spheatp / dt * V_GetCmp (&xT0, element + 1);
        }

        sumT = 0.0;
        for (j = 0; j < n; j++) // 显式处理当前单元 各界面相邻单元的影响
        {
            sumT += aen[j] * V_GetCmp (&xT0, ani[j] + 1); // sumT = 求和[ anT * Tn0 ]
        }
        V_SetCmp (&xT, element + 1, (bep - sumT) / aep); // 显式求解单元的温度值 , Tp = (bep-sumT) / aep
    }
}

void
CalculateTemperature (char *var, int *fiter, double dt, double maxCp,
                      int verbose, int pchecks)  // 求解流场的温度场分布
{
    unsigned int i;
    double mres; // 矩阵求解的残差
    int miter; // 矩阵求解的迭代次数
    double mtime; // 矩阵求解的计算时间

    double tempc; // 温度收敛参数——温度场两次迭代值间的差异

    if (parameter.calc[iT] == LOGICAL_FALSE)  // 若温度场不需要求解，则返回
        return;
    // 创建 Vector类型的向量 Tp
    V_Constr (&xTp, "Temperature at cell center - previous iteration",
              nbelements, Normal, True);

    // Store previous time step values
    Asgn_VV (&xT0, &xT); // 存储上一时间步的变量值

    fiter[iT]++; // 温度场(T) 求解更新的次数

    for (i = 0; i <= parameter.northocor; i++) // 依据 非正交修正次数 多次计算
    {
        if (parameter.timemethod[iT] == IMPLICITEULER  // 当求解采用隐式格式、延时修正格式时
            || parameter.timemethod[iT] == CRANKNICOLSON)
        {
            // 创建 QMatrix类型的矩阵 Ae
            Q_Constr (&Ae, "Energy matrix", nbelements, False,
                      Rowws, Normal, True);
            // 创建 Vector类型的向量 bT
            V_Constr (&bT, "Energy source", nbelements, Normal, True);
        }

        // Store previous iteration values
        if (parameter.northocor > 0) // 当需要多次迭代求解,且本次求解不是最后一次
        {
            Asgn_VV (&xTp, &xT); // 存储上一次迭代求解的值
        }

        if (parameter.timemethod[iT] == EXPLICITEULER) // 若能量守恒方程选择 显式求解
        {
            SolveEnergyExplicit (dt); // 直接显式求解
        }

        if (parameter.timemethod[iT] == IMPLICITEULER) // 若能量守恒方程选择 隐式求解
        {
            BuildEnergyMatrix (dt, 1.0); // sf因子取 1.0 , 全隐式求解
        }

        if (parameter.timemethod[iT] == CRANKNICOLSON) // 若能量守恒方程选择 延时修正求解
        {
            BuildEnergyMatrix (dt, 0.5); // sf因子取 0.5 , 延时修正
        }

        if (parameter.timemethod[iT] == IMPLICITEULER  // 当求解采用隐式格式、延时修正格式时
            || parameter.timemethod[iT] == CRANKNICOLSON)
        {
            if (pchecks == LOGICAL_TRUE) // 开启对计算参数的检查
            {
                if (!CheckIfDiagonalMatrix (&Ae)) // 检查 Ae矩阵是否具有 对角占优性
                {
                    printf("\nWarning: Energy matrix is not diagonal dominant\n");
                    WriteMatrix (&Ae, LOGICAL_FALSE); // 若不是，则将矩阵 Ae,向量 bT写入文件中
                    WriteVector (&bT);
                    //exit (LOGICAL_ERROR);
                }
            }

            // Set matrix solution accuracy
            SetRTCAccuracy (parameter.mtol[iT]); // 设置求解的终止残差

            // Solve matrix to get temperature T  求解能量守恒方程得到 温度场 T
            SolveMatrix (&Ae, &xT, &bT, &miter, &mres, &mtime,parameter.msolver[iT],
                         parameter.mprecond[iT],parameter.miter[iT]);

            if (verbose == LOGICAL_TRUE) // 在命令行中打印求解信息(每求解一次能量守恒方程)
                printf("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
                       var[iT], miter, mres);

            // 矩阵求解结束的残差大于设定的终止残差 或 求解达到最大迭代次数, 保证本次迭代求解是有效的
            if ((mres > parameter.mtol[iT] && miter == parameter.miter[iT])
                || LASResult () != LASOK)  // 线性代数运算的结果标志
            {
                printf ("\nError: Problem solving matrix %c\n", var[iT]); // 打印 "矩阵求解出现问题"
                exit (LOGICAL_ERROR); // 程序退出，并返回 LOGICAL_ERROR
            }
        }

        tempc = 0.0; // 温度收敛参数 初始化

        // Calculate temperature convergence  计算温度场的收敛参数
        if (parameter.northocor > 0) // 当需要多次迭代求解,且本次求解不是最后一次
        {
            tempc = l2Norm_V (Sub_VV (&xTp, &xT)); // 使用l2范数 度量两次温度场迭代向量间的差异
        }

        if (verbose == LOGICAL_TRUE) // 将tempc作为第 i 次修正后的非正交误差, 在命令行中打印出来
            printf ("\nNon-orthogonality error %d (energy): %+E\n", i, tempc);

        CorrectFaceT (&xT, &xTf); // 修正界面上的温度值

        if (parameter.timemethod[iT] == IMPLICITEULER // 当求解采用隐式格式、延时修正格式时
            || parameter.timemethod[iT] == CRANKNICOLSON)
        {
            Q_Destr (&Ae); // 计算完毕，释放迭代计算所创建矩阵、向量
            V_Destr (&bT);
        }

        if (tempc < parameter.mtol[iT]) // 若温度场的相对收敛偏差 小于 设置的迭代终止残差，跳出循环
            break;
    }
    V_Destr (&xTp); // 释放存储上一次迭代温度场的向量 xTp
}
