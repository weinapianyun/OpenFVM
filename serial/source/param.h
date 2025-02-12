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
  EXPLICITEULER = 0,
  IMPLICITEULER,
  CRANKNICOLSON
} par_timemethod;

enum
{
  sJACOBI = 0,
  sSOR,
  sQMR,
  sGMRES,
  sCG,
  sCGN,
  sCGS,
  sBICG,
  sBICGS
} par_solver; // 迭代标志
// (0-Jacobi, 1-SOR, 2-QMR, 3-GMRES, 4-CG, 5-CGN, 6-CGS, 7-BiCG, 8-BiCGStab)

enum
{
  pNONE = 0,
  pJACOBI,
  pSOR,
  pILU,
  pASM
} par_precond;  // 预处理标志 (0-Null, 1-Jacobi, 2-SOR, 3-ILU)

typedef struct
{

  int phi;

  int fixed;
  int along;

  float value;

} par_probe;

typedef struct
{

  char ulength[3];
  char umass[3];
  char utime[3];
  char uenergy[3];
  char utemperature[3];

  int inertia; // 惯性

  float dfactor;

  float st; // 猜测：时间放缩因子

  int scheme[6]; // 迭代方法

  int restart; // 重启迭代次数(时间步数)

  int wbinary; // 写二进制文件

  int steady; // 流场流动的稳定判断
  int adjdt; // 时间步长调节

  float maxCp; // 最大courant number

  // 判断内外迭代收敛的条件
  float mtol[6]; // 内迭代的终止残差 eps
  int miter[6]; // 动量方程的总迭代次数

  int northocor; // 非正交修正次数，默认为 0
  float orthof; // 网格的正交程度，0 为正交

  float ftol[6]; // 各变量外迭代偏差允许值

  int ncicsamsteps; // CICSAM方法的计算步数
  int ncicsamcor; // CICSAM方法的修正次数
  float kq;
  int nsav; // 保存次数

  int calc[6]; // 判断各流场变量是否计算

  int savflux;

  int fsav[6];
  int csav[6];

  int fvec;
  int cvec;

  int smooth;

  int vortex[3];
  int streamf;

  float t0, t1, dt; // 开始、结束计算的时间，时间步
  float g[3];

  int probe[6];

  int msolver[6]; // 迭代方法
  int mprecond[6]; // 预处理方法

  int timemethod[6]; // 各流场变量的时间离散方法
  float ef[6];

  int fill; // 充满状态
  float pf; // 充填百分比

  int vofastemp; // 用相函数 s 的值初始化 温度 T

} par_parameter;

par_parameter parameter;

int ParImportPAR (char *file);
