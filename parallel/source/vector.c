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
    VecCreateSeq (PETSC_COMM_SELF, n, v);
  else
    {
      VecCreateGhost (PETSC_COMM_WORLD, n, PETSC_DECIDE, nbghosts, ghosts, v);
    }

  VecSetFromOptions (*v);

}

void
V_Destr (Vec * v)
{

  VecDestroy (*v);


}

void
V_SetCmp (Vec * v, int ind, double value)
{

  VecSetValue (*v, ind, value, INSERT_VALUES);

}

void
V_SetAllCmp (Vec * v, double value)
{

  VecSet (*v, value);

  VecAssemblyBegin (*v);
  VecAssemblyEnd (*v);

}

double
V_GetCmp (Vec * v, int ind)
{

  double value;
// VecGetValues(Vec x,PetscInt ni,const PetscInt ix[],PetscScalar y[])
// gets y[i] = x[ix[i]], for i=0,...,ni-1.
  VecGetValues (*v, 1, &ind, &value);//输入v,1,ind,输出value

  return value;

}
