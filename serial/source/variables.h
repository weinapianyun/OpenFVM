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

#include "../laspack/laspack.h"

Vector Co; // 柯朗数
Vector uf; // Uf，单元界面上的速度矢量

Vector dens, visc, thcond, spheat; // 密度、粘度、热导率、比热

// A X = b
Vector xu0, xv0, xw0, xp0, xT0, xs0; // 解向量各项的上一时间步的 单元中心值
Vector xu, xv, xw, xp, xT, xs; // 解向量各项的当前迭代步 单元中心值
Vector xuf, xvf, xwf, xpf, xTf, xsf; // 解向量各项的当前迭代步 单元界面值

QMatrix Am, Ac, Ae, As; // 动量、连续性、能量、VOF 方程的系数矩阵
Vector bu, bv, bw, bp, bT, bs; // 各项的右端向量

Vector hu, hv, hw; // SIMPLE法中的 H 向量
Vector ap; // 半离散动量方程中的 amp , 便于求解连续性方程时对amp参数的调用

Vector betaf; // CICSAM格式中的 beta 因子

Vector xpp, xTp;
//Vector xsm;


