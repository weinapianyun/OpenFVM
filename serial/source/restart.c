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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "variables.h"
#include "ioutils.h"
#include "globals.h"
#include "mesh.h"

int
GetInt (FILE * fp) // 从 fp中读取四个字符char 并转换为 int 类型
{

  int value; // int-4个字节，32位
  // fgetc函数：从fp指向的文件中读取一个字符char 并转换为int型返回，然后把位置标识符向下一个字符移动
  value = fgetc (fp) & 0xFF; // 得到 0~8位
  value |= (fgetc (fp) & 0xFF) << 0x08; // 得到 8~16位
  value |= (fgetc (fp) & 0xFF) << 0x10; // 得到 16~24位
  value |= (fgetc (fp) & 0xFF) << 0x18; // 得到 24~32位

  return (value);

}

float
GetFloat (FILE * fp) // // 从 fp中读取四个字符char 并转换为 float 类型
{

  union // 将不同类型数据存放与一处内存地址，所有的成员都从偏移地址零开始存储
  {
    int int_value; // int-4个字节，32位
    float float_value; // float-4个字节，32位
  } value;

  value.int_value = fgetc (fp) & 0xFF; // 转换成 int型变量存储
  value.int_value |= (fgetc (fp) & 0xFF) << 0x08;
  value.int_value |= (fgetc (fp) & 0xFF) << 0x10;
  value.int_value |= (fgetc (fp) & 0xFF) << 0x18;

  return (value.float_value); // 用 float型变量返回

}

void
PutInt (FILE * fp, int value_in) // 将 int型的 value_in写入 fp中
{

  int new_value;

  union
  {
    int int_value; // int-4个字节，32位
    char char_value[4]; // char-1个字节 x4 , 8位 x4
  } value;

  value.int_value = value_in; // 用 int_value存储 value_in的值

  new_value = value.char_value[0] & 0xFF; // 用char[]分离 int的高低字节进行操作
  new_value |= (value.char_value[1] & 0xFF) << 0x08;
  new_value |= (value.char_value[2] & 0xFF) << 0x10;
  new_value |= (value.char_value[3] & 0xFF) << 0x18;

  fwrite (&new_value, sizeof (int), 1, fp); // 将 new_value的值写入 fp

}

void
PutFloat (FILE * fp, float value_in) // 将 float型的 value_in写入 fp中
{

  int new_value;

  union
  {
    float float_value; // int-4个字节，32位
    char char_value[4]; // char-1个字节 x4 , 8位 x4
  } value;

  value.float_value = value_in; // 用 float_value存储 value_in的值

  new_value = value.char_value[0] & 0xFF; // 用char[]分离 float的高低字节进行操作
  new_value |= (value.char_value[1] & 0xFF) << 0x08;
  new_value |= (value.char_value[2] & 0xFF) << 0x10;
  new_value |= (value.char_value[3] & 0xFF) << 0x18;

  fwrite (&new_value, sizeof (int), 1, fp); // 将 new_value的值写入 fp

}

void
WriteRestart (FILE * fp, double curtime) // 将当前时刻的数据存储到重启文件中
{

  int i;

  int face; // 单元界面编号
  int element; // 网格单元编号

  PutFloat (fp, curtime); // 存储 当前的计算时刻

  PutInt (fp, nbelements); // 存储 网格单元数目ne
  PutInt (fp, nbfaces); // 存储 单元界面数目nf

  for (i = 0; i < nbelements; i++) // 遍历所有网格单元
  {

      element = i;

      PutInt (fp, elements[element].index); // 存储网格单元的编号

      PutFloat (fp, (float) V_GetCmp (&xu0, element + 1)); // 将单元的上一轮迭代变量值 存储到 fp中
      PutFloat (fp, (float) V_GetCmp (&xv0, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xw0, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xp0, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xT0, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xs0, element + 1));

      PutFloat (fp, (float) V_GetCmp (&xu, element + 1)); // 将单元的当前变量值 存储到 fp中
      PutFloat (fp, (float) V_GetCmp (&xv, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xw, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xp, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xT, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xs, element + 1));

  }

  for (i = 0; i < nbfaces; i++) // 遍历所有界面
  {

      face = i;

      PutFloat (fp, (float) V_GetCmp (&uf, face + 1)); // 将单元界面的速度值 存储到 fp中

      PutFloat (fp, (float) V_GetCmp (&xuf, face + 1)); // 将单元界面面心的变量值存储到 fp中
      PutFloat (fp, (float) V_GetCmp (&xvf, face + 1));
      PutFloat (fp, (float) V_GetCmp (&xwf, face + 1));
      PutFloat (fp, (float) V_GetCmp (&xpf, face + 1));
      PutFloat (fp, (float) V_GetCmp (&xTf, face + 1));
      PutFloat (fp, (float) V_GetCmp (&xsf, face + 1));

  }

}

void
ReadRestart (FILE * fp, double *curtime) // 读取重启文件的数据，并更新当前时刻
{

  int i;

  float t; // 文件中数据对应的时刻
  int nf, ne; // 网格单元数目ne 单元界面数目nf

  int face;
  int element;

  int eindex;

  t = GetFloat (fp); // 从 fp中读取 文件中数据对应的时刻

  ne = GetInt (fp); // 从 fp中读取 网格单元数目ne
  nf = GetInt (fp); // 从 fp中读取 单元界面数目nf

  if (nbelements != ne || nbfaces != nf) // 若文件中网格单元、单元界面的数目与正在计算的算例参数不符
  {
      printf ("\nWarning: Problem with restart file\n"); // 打印 重启文件有问题
      return;
  }

  for (i = 0; i < nbelements; i++) // 遍历正在执行的算例的网格单元
    {

      element = i;

      eindex = GetInt (fp); // 将文件中网格的索引读出，位置标识符向下一个字符移动

      V_SetCmp (&xu0, element + 1, GetFloat (fp)); // 从 fp中读取单元的上一轮迭代变量值
      V_SetCmp (&xv0, element + 1, GetFloat (fp));
      V_SetCmp (&xw0, element + 1, GetFloat (fp));
      V_SetCmp (&xp0, element + 1, GetFloat (fp));
      V_SetCmp (&xT0, element + 1, GetFloat (fp));
      V_SetCmp (&xs0, element + 1, GetFloat (fp));

      V_SetCmp (&xu, element + 1, GetFloat (fp)); // 从 fp中读取单元的当前变量值
      V_SetCmp (&xv, element + 1, GetFloat (fp));
      V_SetCmp (&xw, element + 1, GetFloat (fp));
      V_SetCmp (&xp, element + 1, GetFloat (fp));
      V_SetCmp (&xT, element + 1, GetFloat (fp));
      V_SetCmp (&xs, element + 1, GetFloat (fp));

    }

  for (i = 0; i < nbfaces; i++) // 遍历正在执行的算例的单元界面
    {

      face = i;

      V_SetCmp (&uf, face + 1, GetFloat (fp)); // 从 fp中读取单元界面的速度

      V_SetCmp (&xuf, face + 1, GetFloat (fp)); // 从 fp中读取单元界面面心的变量值
      V_SetCmp (&xvf, face + 1, GetFloat (fp));
      V_SetCmp (&xwf, face + 1, GetFloat (fp));
      V_SetCmp (&xpf, face + 1, GetFloat (fp));
      V_SetCmp (&xTf, face + 1, GetFloat (fp));
      V_SetCmp (&xsf, face + 1, GetFloat (fp));

    }

  *curtime = t; // 数据读取后，将 curtime更新为 文件中数据对应的时刻

}
