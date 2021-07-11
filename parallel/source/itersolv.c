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

#include "petscksp.h" //线性系统求解器

#include "itersolv.h"

void
GMRESIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter, double maxtol,
	   double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp); //创建求解器对象

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN); //设置与线性系统相关的矩阵

  KSPSetType (ksp, KSPFGMRES); //设置要使用的 Krylov 子空间方法

  KSPGetPC (ksp, &pc); //设置预处理器矩阵

  PCSetType (pc, PrecondProc); //设置前处理方法

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE); //假设给定的解向量的初始值为零
  //残差范数相对于右侧范数的减少 rtol、残差范数的绝对大小 atol 和残差的相对增加量 dtol
  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp); //访问 KSP的大部分功能

  KSPSolve (ksp, *b, *x); //设置右侧向量和解向量并执行

  KSPGetIterationNumber (ksp, iter); //得到迭代过程停止的迭代次数

  KSPGetResidualNorm (ksp, res); //得到残差数值

  KSPDestroy (ksp); //销毁求解器对象

}

void
CRIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter, double maxtol,
	 double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetType (ksp, KSPCR);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);

}
void
CGIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter, double maxtol,
	double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetType (ksp, KSPCG);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);
    
}

void
CGSIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter, double maxtol,
	 double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetType (ksp, KSPCGS);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);

}

void
BiCGIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter,
	      double maxtol, double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetType (ksp, KSPBICG);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);

}

void
BiCGSTABIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter,
	      double maxtol, double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetType (ksp, KSPBCGS);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);

}
