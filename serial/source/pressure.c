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

#include "pressure.h"

void
CorrectFaceP () // 修正单元的界面上的压力值
{

  unsigned int i, k;

  unsigned int node; // 节点编号

  register unsigned int face, pair; // 界面的编号、共面状态
  register unsigned int element, neighbor; // 当前单元、邻接单元的编号

  //double dNf, dPf;
  double lambda; // 界面插值因子

  double ppl; // 单元近似节点的压力值 (非正交修正)

  double phij; // 界面面心的变量值(压力值)

  msh_vector gradpp; // 单元中心的梯度

  double apj; // 边界网格单元的 ap参数值

  double ghf; // 压力计算的中间变量——体积力的影响

  msh_vector g; // 体积力 f

  g.x = parameter.g[0]; // 流场体积力(重力加速度)赋值
  g.y = parameter.g[1];
  g.z = parameter.g[2];

  for (i = 0; i < nbfaces; i++)  // 遍历所有界面
  {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      if (parameter.orthof != 0.0) // 非正交网格单元
          gradpp = Gradient (&xp, &xpf, LOGICAL_TRUE, element); // CG单元基法计算压力梯度

      if (pair != -1) // 非界面单元
	{

	  neighbor = faces[pair].element;

	  // Cell-based linear interpolation

	  /*
	  dNf = GeoMagVector (GeoSubVectorVector(elements[neighbor].celement, faces[face].cface));
	  dPf = GeoMagVector (GeoSubVectorVector(elements[element].celement, faces[face].cface));

	  lambda = dPf / (dPf + dNf);
	  */

	  lambda = 0.5;

	  phij = V_GetCmp (&xp, neighbor + 1) * lambda +
	    V_GetCmp (&xp, element + 1) * (1.0 - lambda); // 界面面心的压力值

	  V_SetCmp (&xpf, face + 1, phij); // 将界面的压力值放入解向量 xpf中

	}
      else // 边界单元
	{

	  ppl = V_GetCmp (&xp, element + 1); // 获取边界网格单元的压力值 Pp

	  apj = V_GetCmp (&ap, element + 1); // 获取边界网格单元的 ap参数值

	  if (parameter.orthof != 0.0) // 非正交网格
	  {
	      //单元近似节点 p'的压力值  Pp'= Pp + grad(Pp)*(rp'-rp)
	      ppl += parameter.orthof * GeoDotVectorVector
	              (gradpp,GeoSubVectorVector
	              (faces[face].rpl,elements[element].celement));
	  }

	  // dens * dot(g,df)  体积力的影响
	  ghf = V_GetCmp (&dens, element + 1) * GeoDotVectorVector (g, faces[face].d);

	  ppl += ghf; // pf = Pp' + dens * g

	  // 根据边界的特点对相应界面的压力值进行一定修正
	  if (faces[face].bc == PERMEABLE) // 渗流界面
	  {
	      // pf = pf0 * (1-F) + Pp' * F , F为相函数
	      V_SetCmp (&xpf, face + 1, V_GetCmp (&xpf, face + 1) *
	            (1.0 - V_GetCmp (&xs, element + 1)) + ppl * V_GetCmp (&xs, element + 1));

	  }

	  if (faces[face].bc == INLET) // 入口界面
	  {
	      // pf = Pp' - Uf * ap * |df| , 入口存在压力梯度，近似grad(P) = Uf * ap
	      V_SetCmp (&xpf, face + 1, ppl - V_GetCmp (&uf, face + 1) * apj *
	                (faces[face].dj + faces[face].kj));
	  }

	  if (faces[face].bc == WALL ||
	      faces[face].bc == MOVINGWALL ||
	      faces[face].bc == ADIABATICWALL ||
	      faces[face].bc == SLIP || faces[face].bc == SURFACE)
	  {
	      // pf = Pp'
	      V_SetCmp (&xpf, face + 1, ppl);

	  }

	}
  }

}

void
BuildContinuityMatrix (double dt) // 创建连续性方程的矩阵
{

  unsigned int i, j, n, nj; // n为单元的非边界界面数目

  register unsigned int face, pair;
  register unsigned int element, neighbor;

  double acp; // 半离散连续性方程的 acp系数
  double acn[MAXFACES]; // 半离散连续性方程的 acn系数
  unsigned int ani[MAXFACES]; // 存储界面相邻单元的编号
  double bcp; // 半离散连续性方程的 系数bp

  double apj; // 暂存 边界网格单元的 ap参数值

  double Huj, Hvj, Hwj; // 存储界面上的 Huf Hvf Hwf值
  double Hf; // 存储界面上的 Hf值

  //double dNf, dPf;
  double lambda; // 插值因子

  msh_vector gradpp; // 单元 p的压力梯度
  msh_vector gradpn; // 单元 n的压力梯度

  // Equation: div(U) = 0  不可压流体的连续性方程
  // 界面速度 Uf = uf - grad(P)/amf , 其中 uf 为无压力影响下的速度，uf = Huf / amf
  // 将速度带入连续性方程得：  求和[ grad(P)*Af / amf ] = 求和[ uf * Af ]

  for (i = 0; i < nbelements; i++) // 遍历网格单元
  {

      element = i;

      acp = 0.0;
      bcp = 0.0;
      n = 0;

      if (parameter.orthof != 0.0) // 非正交网格单元
          gradpp = Gradient (&xp, &xpf, LOGICAL_TRUE, element); // CG单元基法计算单元P的压力梯度

      for (j = 0; j < elements[element].nbfaces; j++) // 遍历每个网格单元的所有界面
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  if (pair != -1) // 非边界界面
	  {

	      neighbor = faces[pair].element; // 共界面单元的编号

              /*
	      dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
	      dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement,  faces[face].cface));

	      lambda = dPf / (dPf + dNf);
	      */

	      lambda = 0.5; // 界面插值因子
	      // 插值计算界面上的 amf系数 ;  amf = lambda * ap + (1 - lambda) * an
	      apj = V_GetCmp (&ap, neighbor + 1) * lambda + V_GetCmp (&ap, element + 1) * (1.0 - lambda);

	      // 插值计算界面上的 Huf Hvf Hwf 值 ; Hf = lambda * Hp + (1 - lambda) * Hn
	      Huj = V_GetCmp (&hu, neighbor + 1) * lambda + V_GetCmp (&hu, element + 1) * (1.0 - lambda);
	      Hvj = V_GetCmp (&hv, neighbor + 1) * lambda + V_GetCmp (&hv, element + 1) * (1.0 - lambda);
	      Hwj = V_GetCmp (&hw, neighbor + 1) * lambda + V_GetCmp (&hw, element + 1) * (1.0 - lambda);

	      // 计算界面上的 Hf值
	      Hf = Huj * faces[face].n.x + Hvj * faces[face].n.y + Hwj * faces[face].n.z;
	      // acp = 求和[ -Af / (amf * |df|) ]
	      acp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj;
	      // acn = - Af / (amf * |df|)
	      acn[n] = 1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj;

	      ani[n] = neighbor; // 存储共面单元的编号
	      n++; // 统计单元的非边界界面数目

	      V_SetCmp (&uf, face + 1, Hf / apj); // 计算界面上的速度值 uf = Hf / amp

	      bcp += V_GetCmp (&uf, face + 1) * faces[face].Aj; // bcp = 求和[ uf * Af ]

	      // Non-orthogonal correction term  非正交修正项
	      if (parameter.orthof != 0.0)
	          gradpn = Gradient (&xp, &xpf, LOGICAL_TRUE, neighbor); // CG单元基法计算单元N 的压力梯度

	      if (parameter.orthof != 0.0) // 非正交网格单元时 对压力梯度的修正
	      {
	          // [grad(Pn)*(rn'-rn)-grad(Pp)*(rp'-rp)] * Af / amf , 将修正项归入bcp中
	          bcp += -1.0 * parameter.orthof / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj *
	                  (GeoDotVectorVector (gradpn, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement)) -
	                  GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement)));
	      }

	  }
	  else // 边界界面
	  {
	      apj = V_GetCmp (&ap, element + 1); // 取界面所在单元中心的 ap值作为 amp

	      // 取界面所在单元中心的Hu Hv Hw 值作为界面上的 Huf Hvf Hwf值
	      Huj = V_GetCmp (&hu, element + 1);
	      Hvj = V_GetCmp (&hv, element + 1);
	      Hwj = V_GetCmp (&hw, element + 1);
	      Hf = Huj * faces[face].n.x + Hvj * faces[face].n.y + Hwj * faces[face].n.z; // 计算界面上的 Hf值

	      if (faces[face].bc == PERMEABLE) // 渗透边界
	      {
	          // acp = 求和[ -Af / (amf * |df|) * (1-F)] , F 为相函数
              acp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj *
                      (1.0 - V_GetCmp (&xs, element + 1));
              // bcp = 求和[ -Af / (amf * |df|) * pf * (1-F) ]
              bcp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj *
                      V_GetCmp (&xpf, face + 1) * (1.0 - V_GetCmp (&xs, element + 1));
              // uf = Hf / amp * (1 - F)
              V_SetCmp (&uf, face + 1, Hf / apj * (1.0 - V_GetCmp (&xs, element + 1)));
              // bcp = bcp + 求和[ uf * Af ]
              bcp += V_GetCmp (&uf, face + 1) * faces[face].Aj;

              // Non-orthogonal correction term
              if (parameter.orthof != 0.0)
              {
                  // 修正项： grad(Pp)*(rp'-rp)*Af/(amf*|df|)*(1-F)
                  bcp += 1.0 * parameter.orthof / (apj * (faces[face].dj + faces[face].kj)) *
                          (GeoDotVectorVector (gradpp,  GeoSubVectorVector (faces[face].rpl, elements[element].celement))) *
                          faces[face].Aj * (1.0 - V_GetCmp (&xs, element + 1));
              }
	      }

	      if (faces[face].bc == OUTLET) // 出口边界
	      {
              // velocity gradient = 0
              // specified pressure

              // acp = 求和[ -Af / (amf * |df|) ]
              acp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj;
              // bcp = 求和[-Af / (amf * |df|) * pf]
              bcp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj * V_GetCmp (&xpf, face + 1);

              V_SetCmp (&uf, face + 1, Hf / apj); // uf = Hf / amp
              // bcp = 求和[ -Af/(amf * |df|) * pf + uf * Af ]
              bcp += V_GetCmp (&uf, face + 1) * faces[face].Aj;

              // Non-orthogonal correction term
              if (parameter.orthof != 0.0)
              {
                  // 修正项： grad(Pp)*(rp'-rp)*Af/(amf*|df|)
                  bcp += 1.0 * parameter.orthof / (apj * (faces[face].dj + faces[face].kj)) *
                          (GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement)))
                          * faces[face].Aj;
              }
	      }

	      if (faces[face].bc == PRESSUREINLET) // 压力边界
	      {
              // specified pressure
              // velocity gradient = 0

              // acp = 求和[ -Af / (amf * |df|) ]
              acp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj;
              // bcp = 求和[-Af / (amf * |df|) * pf]
              bcp += -1.0 / (apj * (faces[face].dj + faces[face].kj)) * faces[face].Aj * V_GetCmp (&xpf, face + 1);

              V_SetCmp (&uf, face + 1, Hf / apj); // uf = Hf / amp
              // bcp = 求和[ -Af/(amf * |df|) * pf + uf * Af ]
              bcp += V_GetCmp (&uf, face + 1) * faces[face].Aj;

              // Non-orthogonal correction term
              if (parameter.orthof != 0.0)
              {
                  bcp += 1.0 * parameter.orthof / (apj * (faces[face].dj + faces[face].kj)) *
                          (GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement)))
                          * faces[face].Aj;
              }
	      }

	      if (faces[face].bc == INLET ||    // 其他边界
		  faces[face].bc == MOVINGWALL ||
		  faces[face].bc == WALL ||
		  faces[face].bc == ADIABATICWALL ||
		  faces[face].bc == SURFACE)
	      {
              // pressure gradient = 0
              // specified velocity

              V_SetCmp (&uf, face + 1,V_GetCmp (&xuf, face + 1) *
                        faces[face].n.x + V_GetCmp (&xvf, face + 1) * faces[face].n.y +
                        V_GetCmp (&xwf, face + 1) * faces[face].n.z);
              bcp += V_GetCmp (&uf, face + 1) * faces[face].Aj; // bcp = 求和[uf * Af]
	      }
	  }
	}

      if (acp == 0.0 || acp != acp) // 若 aCp类型改变 或 取值为0
      {
          printf ("\nError: Problem setting up continuity matrix\n");
          exit (LOGICAL_ERROR);
      }

      nj = 0;

      for (j = 0; j < n; j++)
      {
          if (ani[j] > element) // 矩阵每行中，位于当前网格单元右边的系数数目 nj
              nj++;
      }

      Q_SetLen (&Ac, element + 1, nj + 1); // Am矩阵每行里对角线右侧 系数的数目

      Q_SetEntry (&Ac, element + 1, 0, element + 1, acp);  // 将单元的acp值放入矩阵Ac的对角位置

      nj = 0;

      for (j = 0; j < n; j++)
      {
          if (ani[j] > element) // 矩阵每行中，位于当前网格单元右边的系数数目
          {
              Q_SetEntry (&Ac, element + 1, nj + 1, ani[j] + 1, acn[j]); // 将acn[j]的值放入Am的对应位置
              nj++;
          }
      } // 最终会构成一个上三角矩阵Ac

      V_SetCmp (&bp, element + 1, bcp); // 设置bp向量对应单元处的值

  }

}

void
CalculatePressure (char *var, int *fiter, double dt,
                   double maxCp, int verbose, int pchecks) // 计算压力场
{
  unsigned int i, j;

  double mres; // 矩阵求解的残差
  int miter; // 矩阵求解的迭代次数
  double mtime; // 矩阵求解的计算时间

  double presc; // 压力收敛参数——压力场两次迭代值间的差异

  if (parameter.calc[ip] == LOGICAL_FALSE) // 若压力场不需要计算
      return;
  // 创建 Vector类型的向量，存储上一次迭代的压力场
  V_Constr (&xpp, "Pressure at cell center - previous iteration", nbelements, Normal, True);

  // Store previous time step values
  Asgn_VV (&xp0, &xp); // 更新 xp0 中上一时间步的压力场

  fiter[ip]++; // 压力场的求解更新次数(时间步) +1

  for (i = 0; i <= parameter.northocor; i++) // 依据 非正交修正次数 进行循环计算
  {
      // 创建 QMatrix 类型的矩阵Ac , Vector类型的向量 bp
      Q_Constr (&Ac, "Continuity matrix", nbelements, True, Rowws, Normal, True);
      V_Constr (&bp, "Continuity source", nbelements, Normal, True);

      // Store previous iteration values
      if (parameter.northocor > 0)
      {
        Asgn_VV (&xpp, &xp); // 存储上一次迭代的压力场
      }

      // Build the continuity matrix (mass conservation)
      BuildContinuityMatrix (dt); // 创建连续性方程的矩阵

      if (pchecks == LOGICAL_TRUE)
      {
          if (!CheckIfDiagonalMatrix (&Ac)) // 检查 Ac 是否为 对角占优矩阵
          {
              printf("\nWarning: Continuity matrix is not diagonal dominant\n");
              // 若不是，则将矩阵 Ac,向量 bp写入文件中
              WriteMatrix (&Ac, LOGICAL_TRUE);
              WriteVector (&bp);
              //exit (LOGICAL_ERROR);
          }
      }

      // Set matrix solution accuracy  设置求解的终止残差
      SetRTCAccuracy (parameter.mtol[ip]);

      // Solve matrix to get pressure p  求解连续性方程得到 压力场P
      SolveMatrix (&Ac, &xp, &bp, &miter, &mres, &mtime,parameter.msolver[ip],
                   parameter.mprecond[ip],parameter.miter[ip]);

      if (verbose == LOGICAL_TRUE) // 在命令行中打印求解信息(每求解一次连续性方程)
          printf("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
                 var[ip], miter, mres, mtime);

      // 矩阵求解结束的残差大于设定的终止残差 或 求解达到最大迭代次数, 保证本次迭代求解是有效的
      if ((mres > parameter.mtol[ip] && miter == parameter.miter[ip])
      || LASResult () != LASOK) // 线性代数运算的结果标志
      {
          printf ("\nError: Problem solving matrix %c\n", var[ip]); // 打印 "矩阵求解出现问题"
          exit (LOGICAL_ERROR); // 程序退出，并返回 LOGICAL_ERROR
      }

      presc = 0.0; // 压力收敛参数 初始化

      // Calculate pressure convergence  计算压力场的收敛参数
      if (parameter.northocor > 0)
      {
        presc = l2Norm_V (Sub_VV (&xpp, &xp)); // 使用l2范数 度量两次压力场迭代向量间的差异
      }

      if (verbose == LOGICAL_TRUE) // 将 presc 作为第 i 次修正后的非正交误差
          printf ("\nNon-orthogonality error %d (continuity): %+E\n", i, presc);

      CorrectFaceP (); // 修正界面上的压力值

      Q_Destr (&Ac); // 计算完毕，释放迭代计算所创建矩阵、向量
      V_Destr (&bp);

      if (presc < parameter.mtol[ip]) // 若压力场的相对收敛偏差 小于 设置的迭代终止残差，跳出循环
          break;
  }

  V_Destr (&xpp); // 释放存储上一次迭代压力场的向量xpp

}
