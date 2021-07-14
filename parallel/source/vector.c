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

// PETSc  基于PETSc的函数集类型创建新的函数
#include "petscvec.h"

#include "globals.h"
#include "vector.h"

#include "mesh.h"

void
V_Constr (Vec * v, int n, int sequential)
{

  if (sequential == 1)
    VecCreateSeq (PETSC_COMM_SELF, n, v); // 创建一个标准的、顺序数组样式的向量
  else
    {
      // 在每个处理器上创建一个带有 ghost 填充的并行向量
      VecCreateGhost (PETSC_COMM_WORLD, n, PETSC_DECIDE, nbghosts, ghosts, v);
    }

  VecSetFromOptions (*v); // 从选项数据库配置向量

}

void
V_Destr (Vec * v)
{

  VecDestroy (*v);


}

void
V_SetCmp (Vec * v, int ind, double value)
{

  VecSetValue (*v, ind, value, INSERT_VALUES); // 将向量v 的 ind行 直接设置为 value

}

void
V_SetAllCmp (Vec * v, double value)
{

  VecSet (*v, value); // 将一个向量v 全部设置为 value

  // 在完成对 VecSetValues() 的所有调用后，应调用下面两个函数
  VecAssemblyBegin (*v); // 开始组装向量
  VecAssemblyEnd (*v); // 完成组装向量

}

double
V_GetCmp (Vec * v, int ind)
{

  double value;

  // 从向量的某些位置获取值,目前只能在同一个处理器上获取值
  VecGetValues (*v, 1, &ind, &value);
  // gets y[i] = x[ix[i]], for i=0,...,ni-1.

  return value;

}
