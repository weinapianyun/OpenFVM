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
Gradient (Vector * phi, Vector * phif, int bound, int element)
{

  // Cell based 梯度计算的 CG单元基法

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

  for (j = 0; j < elements[element].nbfaces; j++) // 遍历所计算单元的所有界面
  {

      face = elements[element].face[j];

      pair = faces[face].pair;

      if (pair != -1) // 非边界单元
	{

	  neighbor = faces[pair].element;

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
	  phij = V_GetCmp (phi, neighbor + 1) * lambda +
	          V_GetCmp (phi, element + 1) * (1.0 - lambda);

	  factor =  phij / elements[element].Vp; // 中间变量 Uf / Vp

	  // Element center gradient , grad(U) = 求和（Uf * Af ）/ Vp
	  rv.x += factor * faces[face].A.x;
	  rv.y += factor * faces[face].A.y;
	  rv.z += factor * faces[face].A.z;

	}
      else
	{

	  // Element face variable

	  if (bound == LOGICAL_TRUE)
	    phij = V_GetCmp (phif, face + 1); // 边界界面直接取边界条件的值
	  else
	    phij = V_GetCmp (phi, element + 1); // 直接取网格单元中心的值作为面上的值

	  factor =  phij / elements[element].Vp;

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

      v1 = V_GetCmp (phi, element + 1) + rv.x * (faces[face].cface.x - elements[element].celement.x) + 
					 rv.y * (faces[face].cface.y - elements[element].celement.y) + 
					 rv.z * (faces[face].cface.z - elements[element].celement.z);

      v1min = LMIN(v1min, v1);
      v1max = LMAX(v1max, v1);

      pair = faces[face].pair;

      if (pair != -1)
	{

	  neighbor = faces[pair].element;

          v2 = (V_GetCmp (phi, element + 1) + V_GetCmp (phi, neighbor + 1)) * 0.5;
	  
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
GradientN (Vector * phin, Vector * phif, int bound, int element)
{

  // Node based 梯度计算的 CG 顶点基法

  int j, k;

  int node, face, pair; // 节点编号，界面编号，界面共面状态

  double phij; // 界面面心的变量值

  double factor; // 中间变量

  msh_vector rv; // 存储变量值

  rv.x = 0.0;
  rv.y = 0.0;
  rv.z = 0.0;

  for (j = 0; j < elements[element].nbfaces; j++)
    {

      face = elements[element].face[j];

      pair = faces[face].pair;

      if (pair != -1) // 非边界界面
	{

	  phij = 0.0;

	  for (k = 0; k < faces[face].nbnodes; k++) // 遍历界面包含的节点
	    {

	      node = faces[face].node[k];

	      phij += V_GetCmp (phin, node + 1); // phin 为节点的变量值

	    }

	  // Element face variable

	  if (faces[face].nbnodes > 0)
	      phij /= faces[face].nbnodes; // Uf = 求和(Ufn) / Nf

	  factor =  phij / elements[element].Vp;

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
	      phij = V_GetCmp (phif, face + 1);
	    }
	  else // 由节点变量值计算界面变量值
	    {
	      phij = 0.0;

	      for (k = 0; k < faces[face].nbnodes; k++)
		{

		  node = faces[face].node[k];

		  phij += V_GetCmp (phin, node + 1);

		}

	      // Element face variable

	      if (faces[face].nbnodes > 0)
		phij /= faces[face].nbnodes;

	    }

	  factor =  phij / elements[element].Vp;

	  // Element center gradient

	  rv.x += factor * faces[face].A.x;
	  rv.y += factor * faces[face].A.y;
	  rv.z += factor * faces[face].A.z;

	}

    }

  return rv;

}
