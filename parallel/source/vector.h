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

void V_Constr (Vec * v, int n, int sequential); //创建一个长度为 n 的向量
void V_Destr (Vec * v); // 销毁一个向量
void V_SetCmp (Vec * v, int ind, double value); // 将向量v 的 ind行 直接设置为 value
void V_SetAllCmp (Vec * v, double value); // 将向量v 设置为value值
double V_GetCmp (Vec * v, int ind); // // 根据ind中的1个索引 从向量v中获取1个值，放入value中
