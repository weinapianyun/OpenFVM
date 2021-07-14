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

//PETSc
#include "petscksp.h"

Vec cex, cey, cez; // 单元中心坐标
Vec cexl, ceyl, cezl; // 单元中心坐标中间变量

Vec Co; // 柯朗数
Vec Col;

Vec uf; //面上的速度标量值
Vec dens, visc, spheat, thcond; // 密度、粘度、比热、热导率
Vec densl, viscl, spheatl, thcondl;

// A x = b
Vec xu0, xv0, xw0, xp0, xT0, xs0; // 解向量的上一迭代步的值
Vec xu0l, xv0l, xw0l, xp0l, xT0l, xs0l;

Vec xu, xv, xw, xp, xT, xs; //解向量的当前迭代步待求值
Vec xul, xvl, xwl, xpl, xTl, xsl;

Vec xuf, xvf, xwf, xpf, xTf, xsf;//解向量的当前迭代步 单元界面上的待求值

Mat Am, Ac, Ae, As; // 动量方程、连续性方程、能量方程、VOF对流方程的系数矩阵

Vec bu, bv, bw, bp, bT, bs; // 代数方程组的右端向量

Vec hu, hv, hw; // SIMPLE法中的 H 项
Vec hul, hvl, hwl;

Vec ap; // 离散化方程中的系数项
Vec apl;

Vec xpp, xTp;

Vec betaf; // CICSAM格式中的 beta 因子

Vec temp1, temp2; // 中间变量
