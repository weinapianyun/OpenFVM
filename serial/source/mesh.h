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
  BEAM,
  TRIANGLE,
  QUADRANGLE,
  TETRAHEDRON,
  HEXAHEDRON,
  PRISM
} msh_type;

typedef struct
{

  double x, y, z;

} msh_vector;

typedef struct
{

  int index;

  int type; // 类型

  int nbnodes; // 节点数目
  int *node; // 指向节点

  int element; // 所属网格编号

  msh_vector cface; // 质心位置
  int pair; // 共面界面的编号

  msh_vector n; // 面的法向单位矢量

  msh_vector A; // 面的法向向量
  double Aj; // 面的面积

  msh_vector d; // 单元 P、N 辅助节点间的位移矢量 rp'-rN'
  double dj; // rp'-rN' 矢量的模
  double kj; // 非正交修正量，默认为 0.0
  msh_vector rpl; // 单元 P 辅助节点P'的位矢
  msh_vector rnl; // 单元 N 辅助节点N'的位矢

  int physreg; // 物理编号
  int elemreg; // 网格编号

  int bc; // 界面边界类型

} msh_face;

typedef struct
{

  int index;

  int type;

  msh_vector normal;
  msh_vector celement; // 网格中心

  int nbnodes;
  int *node;

  int nbfaces;
  int *face;

  double dp;
  double Lp;
  double Ap;
  double Vp;

  int physreg; // 物理编号
  int elemreg; // 网格编号

  int process; // 进程编号

  int bc; // 边界条件类型判别

  double *b;
  double *c;
  double *d;

} msh_element;

typedef struct
{

  msh_vector normal;
  double D;

} msh_plane;

// Mesh

int nbnodes;
int nbfaces;
int nbelements;
int nbpatches;

int nbtris;
int nbquads;
int nbtetras;
int nbhexas;
int nbprisms;

msh_vector *nodes;
msh_face *faces;
msh_element *elements;
msh_face *patches;

int nod_correlation_malloced;
int ele_correlation_malloced;

int *nod_correlation;
int *ele_correlation;

int MshImportMSH (char *file);
int MshExportMSH (char *file);

int MshExportDecomposedMSH (char *file, int region, int nbregions);

int MshFreeMemory ();
