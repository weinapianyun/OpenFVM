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

void
GetLine (FILE * fp) // 从文件中读取一行字符串
{
    char one_char;
    one_char = fgetc (fp); // 从流 fp 获取下一个字符
    // 直到读取到最后一个字符，或者换行符
    do{
        one_char = fgetc (fp);
    } while (!feof (fp) && one_char != 10);
}


void
GetString (FILE * fp, char string[256]) // 从文件中读取一个字符串
{
    int i;
    char one_char;
    i = 0;
    do{
        one_char = fgetc (fp);
        if (one_char == 10)
            break;
        string[i] = (char) one_char;
        i++;
        if (i >= 512)
            break;
    } while (!feof (fp));
    string[i] = '\0';
}
