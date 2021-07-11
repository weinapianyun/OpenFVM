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
GetInt (FILE * fp)
{

  int value;
  value = fgetc (fp) & 0xFF;
  value |= (fgetc (fp) & 0xFF) << 0x08;
  value |= (fgetc (fp) & 0xFF) << 0x10;
  value |= (fgetc (fp) & 0xFF) << 0x18;

  return (value);

}

float
GetFloat (FILE * fp)
{

  union
  {
    int int_value;
    float float_value;
  } value;

  value.int_value = fgetc (fp) & 0xFF;
  value.int_value |= (fgetc (fp) & 0xFF) << 0x08;
  value.int_value |= (fgetc (fp) & 0xFF) << 0x10;
  value.int_value |= (fgetc (fp) & 0xFF) << 0x18;

  return (value.float_value);

}

void
PutInt (FILE * fp, int value_in)
{

  int new_value;

  union
  {
    int int_value;
    char char_value[4];
  } value;

  value.int_value = value_in;

  new_value = value.char_value[0] & 0xFF;
  new_value |= (value.char_value[1] & 0xFF) << 0x08;
  new_value |= (value.char_value[2] & 0xFF) << 0x10;
  new_value |= (value.char_value[3] & 0xFF) << 0x18;

  fwrite (&new_value, sizeof (int), 1, fp);

}

void
PutFloat (FILE * fp, float value_in)
{

  int new_value;

  union
  {
    float float_value;
    char char_value[4];
  } value;

  value.float_value = value_in;

  new_value = value.char_value[0] & 0xFF;
  new_value |= (value.char_value[1] & 0xFF) << 0x08;
  new_value |= (value.char_value[2] & 0xFF) << 0x10;
  new_value |= (value.char_value[3] & 0xFF) << 0x18;

  fwrite (&new_value, sizeof (int), 1, fp);

}

void
WriteRestart (FILE * fp, double curtime)
{

  int i;

  int face;
  int element;

  PutFloat (fp, curtime);

  PutInt (fp, nbelements);
  PutInt (fp, nbfaces);

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      PutInt (fp, elements[element].index);

      PutFloat (fp, (float) V_GetCmp (&xu0, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xv0, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xw0, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xp0, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xT0, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xs0, element + 1));

      PutFloat (fp, (float) V_GetCmp (&xu, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xv, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xw, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xp, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xT, element + 1));
      PutFloat (fp, (float) V_GetCmp (&xs, element + 1));

    }

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      PutFloat (fp, (float) V_GetCmp (&uf, face + 1));

      PutFloat (fp, (float) V_GetCmp (&xuf, face + 1));
      PutFloat (fp, (float) V_GetCmp (&xvf, face + 1));
      PutFloat (fp, (float) V_GetCmp (&xwf, face + 1));
      PutFloat (fp, (float) V_GetCmp (&xpf, face + 1));
      PutFloat (fp, (float) V_GetCmp (&xTf, face + 1));
      PutFloat (fp, (float) V_GetCmp (&xsf, face + 1));

    }

}

void
ReadRestart (FILE * fp, double *curtime)
{

  int i;

  float t;
  int nf, ne;

  int face;
  int element;

  int eindex;

  t = GetFloat (fp);

  ne = GetInt (fp);
  nf = GetInt (fp);

  if (nbelements != ne || nbfaces != nf)
    {
      printf ("\nWarning: Problem with restart file\n");
      return;
    }

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      eindex = GetInt (fp);

      V_SetCmp (&xu0, element + 1, GetFloat (fp));
      V_SetCmp (&xv0, element + 1, GetFloat (fp));
      V_SetCmp (&xw0, element + 1, GetFloat (fp));
      V_SetCmp (&xp0, element + 1, GetFloat (fp));
      V_SetCmp (&xT0, element + 1, GetFloat (fp));
      V_SetCmp (&xs0, element + 1, GetFloat (fp));

      V_SetCmp (&xu, element + 1, GetFloat (fp));
      V_SetCmp (&xv, element + 1, GetFloat (fp));
      V_SetCmp (&xw, element + 1, GetFloat (fp));
      V_SetCmp (&xp, element + 1, GetFloat (fp));
      V_SetCmp (&xT, element + 1, GetFloat (fp));
      V_SetCmp (&xs, element + 1, GetFloat (fp));

    }

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      V_SetCmp (&uf, face + 1, GetFloat (fp));

      V_SetCmp (&xuf, face + 1, GetFloat (fp));
      V_SetCmp (&xvf, face + 1, GetFloat (fp));
      V_SetCmp (&xwf, face + 1, GetFloat (fp));
      V_SetCmp (&xpf, face + 1, GetFloat (fp));
      V_SetCmp (&xTf, face + 1, GetFloat (fp));
      V_SetCmp (&xsf, face + 1, GetFloat (fp));

    }

  *curtime = t;

}
