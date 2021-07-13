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
  BEAM, //横梁
  TRIANGLE, //三角形
  QUADRANGLE, //四边形
  TETRAHEDRON, //四面体
  HEXAHEDRON, //六面体
  PRISM //棱柱(五面体)
} msh_type; //网格种类

typedef struct
{

  double x, y, z;

} msh_vector; //网格单元方向矢量

typedef struct
{

  int index; //索引

  int type; //种类

  int nbnodes; //节点数目
  int *node;

  int element; //单元

  msh_vector cface;
  int pair; //相邻单元数目

  msh_vector n;

  msh_vector A; //网格面矢量
  double Aj;

  msh_vector d;
  double dj;
  double kj;
  msh_vector rpl; //中心位移矢量
  msh_vector rnl;

  int physreg; //物理编号
  int elemreg; //单元编号
  int partition; //网格分割信息

  int ghost; // 精灵单元，进程间网格边界的通信

  int bc;

} msh_face; //网格面的信息

typedef struct
{

  int index;

  int type;

  msh_vector normal;
  msh_vector celement;

  int nbnodes;//节点数目
  int *node;//节点指针

  int nbfaces;//界面数目
  int *face;//网格界面的指针

  double dp;
  double Lp;
  double Ap;//界面面积
  double Vp;//网格体积

  int physreg;
  int elemreg;
  int partition;

  int process; //进程编号

  int bc;//边界条件

  double *b; //分配初始化后 网格单元内存空间
  double *c;
  double *d;

} msh_element; //网格单元的信息

typedef struct
{

  msh_vector normal;
  double D;

} msh_plane;

// Mesh

int nbnodes; //nb == number
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
msh_face *patches;
msh_element *elements;

// Ghosts
int nbghosts;
int *ghosts;

int nod_correlation_malloced;
int ele_correlation_malloced;

int *nod_correlation;
int *ele_correlation;

int MshImportMSH (char *file); //读取导入网格文件(.msh)的函数
int MshExportMSH (char *file); //导出网格文件(.msh)的函数
//导出分解后的网格文件的函数
int MshExportDecomposedMSH (char *file, int region, int nbregions); 

void MshFreeMemory (); //网格内存释放函数