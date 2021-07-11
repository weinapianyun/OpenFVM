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

enum
{
  UDS = 0,
  CDS
} par_scheme;

enum
{
  EXPLICITEULER = 0, //显式
  IMPLICITEULER, //隐式
  CRANKNICOLSON // C-N格式
} par_timemethod; //时间离散方法

enum
{
  sGMRES = 3,
  sCR,
  sCG,
  sCGS,
  sBICG,
  sBICGS
} par_solver; //代数矩阵迭代算法(系列共轭梯度算法)

enum
{
  pNONE = 0,
  pJACOBI,  //Jacobi迭代
  pSOR,  //松弛迭代
  pILU,  //LU分解
  pASM  
} par_precond; //预处理子

typedef struct
{

  int phi;

  int fixed;
  int along;

  float value;

} par_probe; //？？

typedef struct
{

  char ulength[3];
  char umass[3];
  char utime[3];
  char uenergy[3];
  char utemperature[3];

  int inertia; //惯性
 
  float dfactor; 

  float st;

  int restart;

  int wbinary;

  int steady;
  int scheme[6];
  int adjdt;

  float maxCp;

  float mtol[6];
  int miter[6];

  int northocor;
  float orthof;  //正交程度

  float ftol[6];

  int ncicsamsteps;
  int ncicsamcor;
  float kq;

  int nsav;

  int calc[6];

  int savflux;

  int fsav[6];
  int csav[6];

  int fvec;
  int cvec;

  int smooth;

  int vortex[3];
  int streamf;

  float t0, t1, dt;
  float g[3];

  int probe[6];

  int msolver[6];
  int mprecond[6];

  int timemethod[6];
  float ef[6];

  int fill;
  float pf;

} par_parameter; //求解计算参数设置

par_parameter parameter;

int ParImportPAR (char *file); //计算参数文件(.par)的读取
