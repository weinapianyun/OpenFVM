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

#include "velocity.h"

void
CorrectVelocityField () // 修正网格单元中心的速度值
{

  unsigned int i, j;

  // 从寄存器中读取 编号数据
  register unsigned int face;
  register unsigned int element;

  msh_vector gradp; // 临时存储压力梯度值

  msh_vector sum1; // 存储 界面法向向量各分量的 ABS值之和
  msh_vector sum2; // 存储 每个面心速度值与各面矢量分量的积 Uf * A.x 之和

  double u, v, w; // 临时存储速度的修正值

  if (parameter.calc[is] == LOGICAL_TRUE) // 若相函数需要计算
  {

      for (i = 0; i < nbelements; i++) // 遍历网格单元
	{

	  element = i;

	  sum1.x = 0.0; // 初始化
	  sum1.y = 0.0;
	  sum1.z = 0.0;

	  sum2.x = 0.0;
	  sum2.y = 0.0;
	  sum2.z = 0.0;

	  for (j = 0; j < elements[element].nbfaces; j++)
	    {

	      face = elements[element].face[j];

	      // 得到单元所有界面法向矢量在 X Y Z 方向投影的绝对值
	      sum1.x += LABS (faces[face].A.x);
	      sum1.y += LABS (faces[face].A.y);
	      sum1.z += LABS (faces[face].A.z);

	      // 得到界面面心速度与各界面法向投影的积 , Uf * A.x = Ufx * | A | , A = |A|*n
	      sum2.x += V_GetCmp (&uf, face + 1) * faces[face].A.x;
	      sum2.y += V_GetCmp (&uf, face + 1) * faces[face].A.y;
	      sum2.z += V_GetCmp (&uf, face + 1) * faces[face].A.z;

	    }

	  // 修正网格单元中心的速度 Ux = 求和(Ufx * |A|) / 求和|A.x|
	  u = sum2.x / sum1.x;
	  v = sum2.y / sum1.y;
	  w = sum2.z / sum1.z;

	  //printf("u: %f\n", u);
	  //printf("v: %f\n", v);
	  //printf("w: %f\n", w);

	  // 用修正后的速度 U 更新解向量 X
	  V_SetCmp (&xu, element + 1, u);
	  V_SetCmp (&xv, element + 1, v);
	  V_SetCmp (&xw, element + 1, w);

	}

  }
  else // 若此时  无需计算相函数
  {

      for (i = 0; i < nbelements; i++) // 遍历网格单元
	{

	  element = i;

	  gradp = Gradient (&xp, &xpf, LOGICAL_TRUE, element); // 计算单元压力梯度值

	  // SIMPLE法-使用动量方程变换的离散形式修正速度场 —— U = (Hu - gradP) / ap
	  V_SetCmp (&xu, element + 1, (V_GetCmp (&hu, element + 1) - gradp.x) / V_GetCmp (&ap, element + 1));
	  V_SetCmp (&xv, element + 1, (V_GetCmp (&hv, element + 1) - gradp.y) / V_GetCmp (&ap, element + 1));
	  V_SetCmp (&xw, element + 1, (V_GetCmp (&hw, element + 1) - gradp.z) / V_GetCmp (&ap, element + 1));

	}

  }

}

void
CorrectFaceUVW () // 修正单元界面面心的速度值
{

  int i;

  unsigned int face, pair; // 网格单元的面编号、共面状态

  register unsigned int element; // 网格单元的单元编号
  unsigned int neighbor; // 相邻单元的编号

  double apj; // 存储 界面面心的 ap 系数

  //double dNf, dPf;
  double lambda; // 界面变量插值因子

  double ppl; // 存储单元 P 的上一轮迭代压力值
  double pnl; // 存储单元 N 的上一轮迭代压力值

  msh_vector gradpp; // 单元 P 的压力梯度
  msh_vector gradpn; // 相邻单元 N 的压力梯度

  for (i = 0; i < nbfaces; i++) // 遍历所有界面
    {

      face = i;

      element = faces[face].element;

      pair = faces[face].pair;

      if (parameter.orthof != 0.0) // 网格非正交
          gradpp = Gradient (&xp, &xpf, LOGICAL_TRUE, element); // 计算单元 P 的压力梯度值

      if (pair != -1) // 非边界界面
	{

	  neighbor = faces[pair].element;

	  /*   // 计算 lambda 的值
	  dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
	  dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement, faces[face].cface));

	  lambda = dPf / (dPf + dNf);
	  */

	  lambda = 0.5;

	  apj = V_GetCmp (&ap, neighbor + 1) * lambda +
	          V_GetCmp (&ap, element + 1) * (1.0 - lambda); // 插值得到界面上的 ap系数值

	  if (parameter.orthof != 0.0)
	    gradpn = Gradient (&xp, &xpf, LOGICAL_TRUE, neighbor); // 计算单元 N 的压力梯度值

	  ppl = V_GetCmp (&xp, element + 1); // 先将单元 P 上一轮压力迭代值Pp 放入 ppl中

	  if (parameter.orthof != 0.0)
	    {
	      // 计算单元 P 的辅助节点P'的压力值 ppl
	      ppl += parameter.orthof * GeoDotVectorVector (gradpp,
                                                        GeoSubVectorVector (faces[face].rpl,
                                                                            elements[element].celement));
	    }

	  pnl = V_GetCmp (&xp, neighbor + 1); // 将单元 N 上一轮压力迭代值Pn 放入 pnl中

	  if (parameter.orthof != 0.0)
	    {
	      // 计算单元 N 的辅助节点N'的压力值 pnl
	      pnl += parameter.orthof * GeoDotVectorVector (gradpn,
                                                        GeoSubVectorVector (faces[face].rnl,
                                                                            elements[neighbor].celement));
	    }

	  // 根据修正的压力值计算界面面心速度 Uf = Ufl - 1/ap*(Pnl-Ppl)/ |df|
	  V_SetCmp (&uf, face + 1, V_GetCmp (&uf, face + 1) -
              1.0 / apj * (pnl - ppl) /
              (faces[face].dj + faces[face].kj)); // df的表达？？

	  // 根据界面的单位法向矢量，由Uf得到 X Y Z 方向的速度分量
	  V_SetCmp (&xuf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.x);
	  V_SetCmp (&xvf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.y);
	  V_SetCmp (&xwf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.z);

	}
      else // 边界界面
	{

	  apj = V_GetCmp (&ap, element + 1); // 将界面所在单元上轮迭代的的ap系数作为界面当前的ap系数

	  ppl = V_GetCmp (&xp, element + 1); // 将单元 P 上一轮压力迭代值Pp 放入 ppl中

	  if (parameter.orthof != 0.0) // 网格非正交
	    {
	      // 计算单元 P 的辅助节点P'的压力值 ppl
	      ppl += parameter.orthof * GeoDotVectorVector (gradpp, GeoSubVectorVector (faces[face].rpl, elements[element].celement));
	    }

	  // 根据边界面的性质来修正界面的速度场
	  if (faces[face].bc == PERMEABLE)
	    {
	      //  Uf = Ufl - 1/(apj * |df|) * (Ppf - ppl)
	      V_SetCmp (&uf, face + 1, (V_GetCmp (&uf, face + 1) - 1.0 /
	        (apj * (faces[face].dj + faces[face].kj)) * (V_GetCmp (&xpf, face + 1) - ppl)));

	      // Wall
	      /*
	      V_SetCmp (&xuf, face + 1, V_GetCmp (&uf,  face + 1) * faces[face].n.x * (1.0 - V_GetCmp(&xs, element + 1)));
	      V_SetCmp (&xvf, face + 1, V_GetCmp (&uf,  face + 1) * faces[face].n.y * (1.0 - V_GetCmp(&xs, element + 1)));
	      V_SetCmp (&xwf, face + 1, V_GetCmp (&uf,  face + 1) * faces[face].n.z * (1.0 - V_GetCmp(&xs, element + 1)));
	      */
		
	      // Slip  Uf = U * (1 - gamma)
	      V_SetCmp (&xuf, face + 1, V_GetCmp (&xu, element + 1) * (1.0 - V_GetCmp(&xs, element + 1)));
	      V_SetCmp (&xvf, face + 1, V_GetCmp (&xv, element + 1) * (1.0 - V_GetCmp(&xs, element + 1)));
	      V_SetCmp (&xwf, face + 1, V_GetCmp (&xw, element + 1) * (1.0 - V_GetCmp(&xs, element + 1)));

	    }

	  if (faces[face].bc == OUTLET) // 出口边界
	    {

	      // velocity gradient = 0
	      // specified pressure

	      V_SetCmp (&uf, face + 1, V_GetCmp (&uf, face + 1) - 1.0 / (apj * (faces[face].dj + faces[face].kj)) * (V_GetCmp (&xpf, face + 1) - ppl));

	      V_SetCmp (&xuf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.x); 
	      V_SetCmp (&xvf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.y);
	      V_SetCmp (&xwf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.z);

	    }

	  if (faces[face].bc == PRESSUREINLET) // 恒压入口
	    {

	      // velocity gradient = 0
	      // specified pressure

	      V_SetCmp (&uf, face + 1, V_GetCmp (&uf, face + 1) - 1.0 / (apj * (faces[face].dj + faces[face].kj)) * (V_GetCmp (&xpf, face + 1) - ppl));

	      V_SetCmp (&xuf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.x);
	      V_SetCmp (&xvf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.y);
	      V_SetCmp (&xwf, face + 1, V_GetCmp (&uf, face + 1) * faces[face].n.z);

	    }

	  if (faces[face].bc == INLET ||  // 恒流入口
	      faces[face].bc == MOVINGWALL ||
	      faces[face].bc == WALL ||
	      faces[face].bc == ADIABATICWALL || faces[face].bc == SURFACE)
	    {

	      // pressure gradient = 0
	      // specified velocity

	      V_SetCmp (&uf, face + 1, V_GetCmp (&xuf, face + 1) * faces[face].n.x + V_GetCmp (&xvf, face + 1) * faces[face].n.y + V_GetCmp (&xwf, face + 1) * faces[face].n.z);

	      V_SetCmp (&xpf, face + 1, ppl);

	    }

	  if (faces[face].bc == SLIP) // 滑移面
	    {

	      // pressure gradient = 0
	      // velocity gradient = 0

	      V_SetCmp (&xuf, face + 1, V_GetCmp (&xu, element + 1));
	      V_SetCmp (&xvf, face + 1, V_GetCmp (&xv, element + 1));
	      V_SetCmp (&xwf, face + 1, V_GetCmp (&xw, element + 1));

	      V_SetCmp (&uf, face + 1, 0.0);

	    }

	}

    }

}

void
CalculateCorrectionFactors () // 计算速度求解的系数 Hu  Hv Hw —— 近似代表无压力项时的速度分布
{
  // BuildMomentumMatrix 函数中已对 Hu 进行初始化 , Hu = Su
  // Hu = Hu - 求和[aN * UN] = -求和[aN * UN] + Su
  AddAsgn_VV (&hu, Mul_QV (Sub_QQ (Diag_Q (&Am), &Am), &xu));
  AddAsgn_VV (&hv, Mul_QV (Sub_QQ (Diag_Q (&Am), &Am), &xv));
  AddAsgn_VV (&hw, Mul_QV (Sub_QQ (Diag_Q (&Am), &Am), &xw));

  //PrintVector(&xu);
  //PrintVector(&hu);

}

void
BuildMomentumMatrix (double dt)
{

  int i, j, n; // 计数

  unsigned int face, pair;
 
  register unsigned int element;
  unsigned int neighbor;

  double app; // 半离散方程的系数ap

  double apn[MAXFACES]; // 半离散方程的系数aN
  unsigned int ani[MAXFACES]; // 存储界面相邻单元的编号

  double bpu, bpv, bpw; // 半离散方程的系数bup

  double densp; // 单元密度
  double viscj; // 单元各界面的粘度

  //msh_vector gradup, gradvp, gradwp;
  //msh_vector gradun, gradvn, gradwn;
  //msh_vector gradvisc;

  msh_vector gradp;

  //double dNf, dPf;
  double lambda; // 插值调节因子
  double xsi; // 差分调节因子

  msh_vector g; // 体积力

  g.x = parameter.g[0]; // 流场体积力(重力加速度)赋值
  g.y = parameter.g[1];
  g.z = parameter.g[2];

  // Equation: dU/dt + div(rho*U*U) - div(mi*grad(U)) = qU

  for (i = 0; i < nbelements; i++) // 遍历单元
  {

      element = i;

      bpu = 0.0; // 初始化
      bpv = 0.0;
      bpw = 0.0;

      app = 0.0; // 离散化方程的系数

      n = 0; // 记录每个单元的界面数

      /*  计算单元速度梯度
         gradup = Gradient (&xu, &xuf, LOGICAL_TRUE, element);
         gradvp = Gradient (&xv, &xvf, LOGICAL_TRUE, element);
         gradwp = Gradient (&xw, &xwf, LOGICAL_TRUE, element);
       */

      densp = V_GetCmp (&dens, element + 1); // 单元的密度

      for (j = 0; j < elements[element].nbfaces; j++) // 遍历单元界面
	{

	  face = elements[element].face[j];

	  pair = faces[face].pair;

	  if (pair != -1) // 非边界界面
	    {

	      neighbor = faces[pair].element; // 相邻单元的编号

              /*
	      dNf = GeoMagVector (GeoSubVectorVector (elements[neighbor].celement, faces[face].cface));
	      dPf = GeoMagVector (GeoSubVectorVector (elements[element].celement,  faces[face].cface));

	      lambda = dPf / (dPf + dNf);
              */
          lambda = 0.5; // 插值调节因子

          // 界面的粘度 插值获得
	      viscj = V_GetCmp (&visc, element + 1) * (1.0 - lambda) + V_GetCmp (&visc, neighbor + 1) * lambda;

	      // Convection 对流项
	      if (parameter.scheme[iu] == UDS)
		{

		  // UDS  上风差分格式
		  if (V_GetCmp (&uf, face + 1) > 0.0)
		    xsi = 0.0;
		  else
		    xsi = 1.0;

		}
	      else
		{
		  // CDS  中心差分格式
		  xsi = lambda;
		}

	      // Convection - 对流项  1/Vp*求和[ (1-xsi)*dens*Uf*Af ]
	      app += (1.0 - xsi) * densp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;

	      // Diffusion - 扩散项  1/Vp*求和[ viscf * Af / |df| ]
	      app += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

	      // Convection - 对流项  1/Vp*求和[ xsi*dens*Uf*Af ]
	      apn[n] = xsi * densp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;

	      // Diffusion - 扩散项  1/Vp*求和[ - viscf * Af / |df| ]
	      apn[n] += -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

	      ani[n] = neighbor; // 存储相邻单元的编号
	      n++;

	      /*  计算非正交修正项，归入源项, 扩散项中界面面心速度梯度grad(Uf) ——产生—— 非正交修正项
	       * 单元P miu*Af/|df| * dot[grad(Up), rpl - rp]
	       * 单元N -miu*Af/|df| * dot[grad(Un), rnl - rn]

	         if (parameter.orthof != 0.0)
	         {
	         gradun = Gradient (&xu, &xuf, LOGICAL_TRUE, neighbor);
	         gradvn = Gradient (&xv, &xvf, LOGICAL_TRUE, neighbor);
	         gradwn = Gradient (&xw, &xwf, LOGICAL_TRUE, neighbor);

	         // Non-orthogonal correction terms           
	         bpu += parameter.orthof * 
	             -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp *
	         (GeoDotVectorVector (gradun, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement)) -
	         GeoDotVectorVector (gradup, GeoSubVectorVector (faces[face].rpl, elements[element].celement)));

	         bpv += parameter.orthof * 
                 -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp *
	         (GeoDotVectorVector (gradvn, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement)) -
	         GeoDotVectorVector (gradvp, GeoSubVectorVector (faces[face].rpl, elements[element].celement)));

	         bpw += parameter.orthof * 
	             -viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp *
	         (GeoDotVectorVector (gradwn, GeoSubVectorVector (faces[face].rnl, elements[neighbor].celement)) -
		     GeoDotVectorVector (gradwp, GeoSubVectorVector (faces[face].rpl, elements[element].celement)));
	         }
	       */

	    }
	  else // 边界界面
	    {

	      if (faces[face].bc != EMPTY)
		{

		  viscj = V_GetCmp (&visc, element + 1);

		  // Diffusion
		  if ((faces[face].bc != PERMEABLE || V_GetCmp (&xs, element + 1) > 0.5) && faces[face].bc != SLIP)
		    {
		      // 1/Vp * miuf*Af/|df| , 边界处归入扩散项的部分
		      app += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp;

		      // 1/Vp * miuf*Af/|df|*Uf , 边界处归入源项的部分
		      bpu += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) * V_GetCmp (&xuf, face + 1) / elements[element].Vp;
		      bpv += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) * V_GetCmp (&xvf, face + 1) / elements[element].Vp;
		      bpw += viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) * V_GetCmp (&xwf, face + 1) / elements[element].Vp;
		    }

		  // Convection
		  if (parameter.scheme[iu] == UDS)
		  {

		      // UDS 迎风格式
		      if (V_GetCmp (&uf, face + 1) > 0.0) // 使用内部网格单元的变量
			{
		      // 1/Vp*densf*uf*Af , 可归入对流项
			  app += densp * V_GetCmp (&uf, face + 1) * faces[face].Aj / elements[element].Vp;

			}
		      else // 使用边界处的变量
			{

			  // Convection    1/Vp*densf*uf*Af*Ubf , 直接归入源项
			  bpu += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xuf, face + 1) / elements[element].Vp;
			  bpv += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xvf, face + 1) / elements[element].Vp;
			  bpw += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xwf, face + 1) / elements[element].Vp;

			}

		  }
		  else
		  {

		      // CDS 中心差分  1/Vp*densf*uf*Af*Ubf  直接归入源项
		      bpu += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xuf, face + 1) / elements[element].Vp;
		      bpv += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xvf, face + 1) / elements[element].Vp;
		      bpw += -densp * V_GetCmp (&uf, face + 1) * faces[face].Aj * V_GetCmp (&xwf, face + 1) / elements[element].Vp;

		  }

		  // Non-orthogonal correction terms
		  /*  扩散项中界面面心速度梯度grad(Uf) ——产生—— 非正交修正项  miu*Af/|df| * dot[grad(Up), rpl - rp]

		     bpu +=
		     viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * 
		     GeoDotVectorVector (gradup, GeoSubVectorVector(faces[face].rpl, elements[element].celement));

		     bpv +=
		     viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * 
		     GeoDotVectorVector (gradvp, GeoSubVectorVector(faces[face].rpl, elements[element].celement));

		     bpw +=
		     viscj * faces[face].Aj / (faces[face].dj + faces[face].kj) / elements[element].Vp * 
		     GeoDotVectorVector (gradwp, GeoSubVectorVector(faces[face].rpl, elements[element].celement));
		   */

		}

	    }

	}

      if (dt > 0) // 瞬态过程，方程中有非稳态项
	{

	  // Unsteady term - Euler 非稳态项

	  app += densp / dt * parameter.st; // dens / dt

	  //  dens * U0 / dt
	  bpu += densp / dt * V_GetCmp (&xu0, element + 1);
	  bpv += densp / dt * V_GetCmp (&xv0, element + 1);
	  bpw += densp / dt * V_GetCmp (&xw0, element + 1);

	}

      // Source - viscous term
      /*  计算源项中的 黏性应力项 dot[grad(Up),grad(miup)]

         gradvisc = Gradient (&visc, NULL, LOGICAL_FALSE, element);

         bpu += GeoDotVectorVector (gradup, gradvisc);
         bpv += GeoDotVectorVector (gradvp, gradvisc);
         bpw += GeoDotVectorVector (gradwp, gradvisc);
       */

      // Source - gravity  体积力源项 dens * g
      bpu += densp * g.x;
      bpv += densp * g.y;
      bpw += densp * g.z;

      // Initialize H with source contribution without pressure , 使用无压力左右的源项初始化 Hu
      // 此时 Hu = Su ,在 CalculateCorrectionFactors函数中 进一步完成Hu的计算
      V_SetCmp (&hu, element + 1, bpu);
      V_SetCmp (&hv, element + 1, bpv);
      V_SetCmp (&hw, element + 1, bpw);

      // Source - pressure 
      /* 计算源项中的压力梯度项 grad(P)，直接使用Gradient函数计算压力梯度

      gradp = Gradient (&xp, &xpf, LOGICAL_TRUE, element);
      bpu += -gradp.x;
      bpv += -gradp.y;
      bpw += -gradp.z;
      */

      if (app == 0.0 || app != app) // 若 app类型改变 或 取值为0
	{
	  printf ("\nError: Problem setting up momentum matrix\n");
	  exit (LOGICAL_ERROR);
	}

      V_SetCmp (&ap, element + 1, app); // 设置计算的参数 ap

      Q_SetLen (&Am, element + 1, n + 1); // 设置 Am 矩阵的规格

      Q_SetEntry (&Am, element + 1, 0, element + 1, app); // 将单元的app值放入矩阵Am的对角位置

      for (j = 0; j < n; j++)
	{
	  Q_SetEntry (&Am, element + 1, j + 1, ani[j] + 1, apn[j]); // 将apn[j] 的值放入Am 的其他位置
	}

      // 设置bu向量对应单元处的值
      if (parameter.calc[iu] == LOGICAL_TRUE)
	{
	  V_SetCmp (&bu, element + 1, bpu);
	}

      if (parameter.calc[iv] == LOGICAL_TRUE)
	{
	  V_SetCmp (&bv, element + 1, bpv);
	}

      if (parameter.calc[iw] == LOGICAL_TRUE)
	{
	  V_SetCmp (&bw, element + 1, bpw);
	}

  }

} // 组装离散化动量守恒方程的代数矩阵 AX=b

void
CorrectVelocity (char *var, int *fiter, double dt, double maxCp, // 修正所有网格单元的速度场
                 int verbose,int pchecks)
{

  if (parameter.calc[ip] == LOGICAL_TRUE) // 求解时计算压力场
    {

      // Correct face values 修正单元界面面心的速度值
      CorrectFaceUVW ();

      // Correct cell center 修正网格单元中心的速度值
      CorrectVelocityField ();

   }
   
}

void
CalculateVelocity (char *var, int *fiter, double dt, double maxCp,  // 求解流场的速度场分布
                   int verbose, int pchecks)
{

  double mres; // 矩阵求解的残差
  int miter; // 矩阵求解的迭代次数
  double mtime; // 矩阵求解的计算时间

  // 构造 QMatrix 类型的矩阵Am , Vector类型的向量 bu bv nw
  Q_Constr (&Am, "Momentum matrix", nbelements, False, Rowws, Normal, True);
  V_Constr (&bu, "Momentum source x-component", nbelements, Normal, True);
  V_Constr (&bv, "Momentum source y-component", nbelements, Normal, True);
  V_Constr (&bw, "Momentum source z-component", nbelements, Normal, True);

  // Store previous time step values  将上一轮迭代计算出的变量值存储在 U0 中
  Asgn_VV (&xu0, &xu);
  Asgn_VV (&xv0, &xv);
  Asgn_VV (&xw0, &xw);

  // Build three momentum matrices for u, v, w velocity components
  if (parameter.calc[ip] == LOGICAL_TRUE)
  {
      BuildMomentumMatrix (dt); // 创建u v w的三个动量方程的矩阵格式
  }

  if (pchecks == LOGICAL_TRUE)
  {
      if (!CheckIfDiagonalMatrix (&Am)) // 检查 Am 是否为 对角占优矩阵
	{
	  printf ("\nWarning: Momentum matrix is not diagonal dominant\n");

	  // 将矩阵 Am 向量 bu bv bw写入文件中
	  WriteMatrix (&Am, LOGICAL_FALSE);
	  WriteVector (&bu);
	  WriteVector (&bv);
	  WriteVector (&bw);
	  //exit (LOGICAL_ERROR);
	}
  }

  // 矩阵求解 信息输出
  if (parameter.calc[iu] == LOGICAL_TRUE) // u方向的动量矩阵
  {
      fiter[iu]++; // u动量方程 迭代求解的次数+1

      // Set matrix solution accuracy  设置求解的终止残差
      SetRTCAccuracy (parameter.mtol[iu]);

      // Solve matrix for u velocity component  求解含 u 的动量方程
      SolveMatrix (&Am, &xu, &bu, &miter, &mres, &mtime, parameter.msolver[iu],
                   parameter.mprecond[iu], parameter.miter[iu]);

      if (verbose == LOGICAL_TRUE) // 在命令行中打印求解信息(求解一次 u动量矩阵)
          printf("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
                 var[iu], miter, mres, mtime);

      // 矩阵求解结束的残差大于设定的终止残差 或 求解达到最大迭代次数
      if ((mres > parameter.mtol[iu] && miter == parameter.miter[iu])
            || LASResult () != LASOK) // 线性代数运算的结果标志
      {
          printf ("\nError: Problem solving matrix %c\n", var[iu]); // 打印 "矩阵求解出现问题"
          exit (LOGICAL_ERROR); // 程序退出，并返回 LOGICAL_ERROR
      }

  }

  if (parameter.calc[iv] == LOGICAL_TRUE) // v方向的动量矩阵
  {

      fiter[iv]++;

      // Set matrix solution accuracy
      SetRTCAccuracy (parameter.mtol[iv]);

      // Solve matrix for v velocity component
      SolveMatrix (&Am, &xv, &bv, &miter, &mres, &mtime, parameter.msolver[iv],
                   parameter.mprecond[iv], parameter.miter[iv]);

      if (verbose == LOGICAL_TRUE)
          printf("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
	   var[iv], miter, mres, mtime);

      if (mres > parameter.mtol[iv] && miter == parameter.miter[iv])
      {
          printf ("\nError: Problem solving matrix %c\n", var[iv]);
          exit (LOGICAL_ERROR);
	  }

  }

  if (parameter.calc[iw] == LOGICAL_TRUE) // w方向的动量矩阵
  {

      fiter[iw]++;

      // Set matrix solution accuracy
      SetRTCAccuracy (parameter.mtol[iw]);

      // Solve matrix for w velocity component
      SolveMatrix (&Am, &xw, &bw, &miter, &mres, &mtime, parameter.msolver[iw],
                   parameter.mprecond[iw], parameter.miter[iw]);

      if (verbose == LOGICAL_TRUE)
          printf("\nMatrix %c Number of iterations: %d Residual: %+E Time: %+E\n",
	   var[iw], miter, mres, mtime);

      if ((mres > parameter.mtol[iw] && miter == parameter.miter[iw])
            ||LASResult () != LASOK)
      {
          printf ("\nProblem solving matrix %c\n", var[iw]);
          exit (LOGICAL_ERROR);
      }

  }

  // Calculate correction factors  速度场求解完毕后，根据当前速度场变量计算 Hu
  if (parameter.calc[ip] == LOGICAL_TRUE)
    {
      // SIMPLE算法中，压力场更新后，进行下一步迭代时 Hu 需要重新计算
      CalculateCorrectionFactors (&Am, &xu, &xv, &xw, &hu, &hv, &hw);
    }

  Q_Destr (&Am); // 计算完毕，释放迭代计算所创建矩阵、向量
  V_Destr (&bu);
  V_Destr (&bv);
  V_Destr (&bw);

}
