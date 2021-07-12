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
  NONE = 0,
  EMPTY, //空的
  CYCLIC, //循环的
  PROCESSOR, //处理器
  OPEN, //开
  INLET, //入口
  PRESSUREINLET, //PRESSUREINLET
  OUTLET, //出口
  ADIABATICWALL, //绝热壁
  MOVINGWALL, //移动墙
  WALL, //墙
  SLIP, //滑移
  PERMEABLE, //透水 
  CONSTRAINTU,
  CONSTRAINTV,
  CONSTRAINTW,
  CONSTRAINT, //约束
  PRESSURE, //压力
  SURFACE, //面
  VOLUME //体积
} bcd_type;

typedef struct
{
  float x, y, z;

} bcd_vector; //边界向量

typedef struct
{

  int physreg;

  int bc;

  char *fu, *fv, *fw;
  char *fp;
  char *fT;
  char *fs;

} bcd_surface; //面物理条件

typedef struct
{

  int physreg;

  int bc;

  char *fu, *fv, *fw;
  char *fp;
  char *fT;
  char *fs;

} bcd_volume; //体物理条件

// Boundary conditions
int nbbcsurfaces; //面的数目
bcd_surface *bcsurfaces;

int nbbcvolumes; //体的数目
bcd_volume *bcvolumes;

int BcdImportBCD (char *file); //读取边界条件(.bcd文件)的函数
