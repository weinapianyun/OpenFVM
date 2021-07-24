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

void PrintVector (Vector * x); // 在命令行中打印向量X
void WriteVector (Vector * x); // 将向量写入vector.m文件中
void WriteMatrix (QMatrix * A, int symmetric); // 将矩阵A 写入文件Matrix.m中
int CheckIfDiagonalMatrix (QMatrix * A); // 检查矩阵A是否 弱对角占优
void SolveMatrix (QMatrix * A, Vector * x, Vector * b, int *iter, double *res,
		  double *time, int msolver, int mprecond, int miter); // 求解线性代数矩阵
