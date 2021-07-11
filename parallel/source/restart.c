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
#include "vector.h"

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

  VecGhostGetLocalForm (xu0, &xu0l);
  VecGhostGetLocalForm (xv0, &xv0l);
  VecGhostGetLocalForm (xw0, &xw0l);
  VecGhostGetLocalForm (xp0, &xp0l);
  VecGhostGetLocalForm (xT0, &xT0l);
  VecGhostGetLocalForm (xs0, &xs0l);

  VecGhostGetLocalForm (xu, &xul);
  VecGhostGetLocalForm (xv, &xvl);
  VecGhostGetLocalForm (xw, &xwl);
  VecGhostGetLocalForm (xp, &xpl);
  VecGhostGetLocalForm (xT, &xTl);
  VecGhostGetLocalForm (xs, &xsl);

  PutFloat (fp, curtime);

  PutInt (fp, nbelements);
  PutInt (fp, nbfaces);

  for (i = 0; i < nbelements; i++)
    {

      element = i;

      PutInt (fp, elements[element].index);

      PutFloat (fp, (float) V_GetCmp (&xu0l, element));
      PutFloat (fp, (float) V_GetCmp (&xv0l, element));
      PutFloat (fp, (float) V_GetCmp (&xw0l, element));
      PutFloat (fp, (float) V_GetCmp (&xp0l, element));
      PutFloat (fp, (float) V_GetCmp (&xT0l, element));
      PutFloat (fp, (float) V_GetCmp (&xs0l, element));

      PutFloat (fp, (float) V_GetCmp (&xul, element));
      PutFloat (fp, (float) V_GetCmp (&xvl, element));
      PutFloat (fp, (float) V_GetCmp (&xwl, element));
      PutFloat (fp, (float) V_GetCmp (&xpl, element));
      PutFloat (fp, (float) V_GetCmp (&xTl, element));
      PutFloat (fp, (float) V_GetCmp (&xsl, element));

    }

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      PutFloat (fp, (float) V_GetCmp (&uf, face));

      PutFloat (fp, (float) V_GetCmp (&xuf, face));
      PutFloat (fp, (float) V_GetCmp (&xvf, face));
      PutFloat (fp, (float) V_GetCmp (&xwf, face));
      PutFloat (fp, (float) V_GetCmp (&xpf, face));
      PutFloat (fp, (float) V_GetCmp (&xTf, face));
      PutFloat (fp, (float) V_GetCmp (&xsf, face));

    }

  VecGhostRestoreLocalForm (xu0, &xu0l);
  VecGhostRestoreLocalForm (xv0, &xv0l);
  VecGhostRestoreLocalForm (xw0, &xw0l);
  VecGhostRestoreLocalForm (xp0, &xp0l);
  VecGhostRestoreLocalForm (xT0, &xT0l);
  VecGhostRestoreLocalForm (xs0, &xs0l);

  VecGhostRestoreLocalForm (xu, &xul);
  VecGhostRestoreLocalForm (xv, &xvl);
  VecGhostRestoreLocalForm (xw, &xwl);
  VecGhostRestoreLocalForm (xp, &xpl);
  VecGhostRestoreLocalForm (xT, &xTl);
  VecGhostRestoreLocalForm (xs, &xsl);

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

      V_SetCmp (&xu0, eindex, GetFloat (fp));
      V_SetCmp (&xv0, eindex, GetFloat (fp));
      V_SetCmp (&xw0, eindex, GetFloat (fp));
      V_SetCmp (&xp0, eindex, GetFloat (fp));
      V_SetCmp (&xT0, eindex, GetFloat (fp));
      V_SetCmp (&xs0, eindex, GetFloat (fp));

      V_SetCmp (&xu, eindex, GetFloat (fp));
      V_SetCmp (&xv, eindex, GetFloat (fp));
      V_SetCmp (&xw, eindex, GetFloat (fp));
      V_SetCmp (&xp, eindex, GetFloat (fp));
      V_SetCmp (&xT, eindex, GetFloat (fp));
      V_SetCmp (&xs, eindex, GetFloat (fp));

    }

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      V_SetCmp (&uf, face, GetFloat (fp));

      V_SetCmp (&xuf, face, GetFloat (fp));
      V_SetCmp (&xvf, face, GetFloat (fp));
      V_SetCmp (&xwf, face, GetFloat (fp));
      V_SetCmp (&xpf, face, GetFloat (fp));
      V_SetCmp (&xTf, face, GetFloat (fp));
      V_SetCmp (&xsf, face, GetFloat (fp));

    }

  VecAssemblyBegin (xu0);
  VecAssemblyEnd (xu0);
  VecAssemblyBegin (xv0);
  VecAssemblyEnd (xv0);
  VecAssemblyBegin (xw0);
  VecAssemblyEnd (xw0);
  VecAssemblyBegin (xp0);
  VecAssemblyEnd (xp0);
  VecAssemblyBegin (xT0);
  VecAssemblyEnd (xT0);
  VecAssemblyBegin (xs0);
  VecAssemblyEnd (xs0);

  VecAssemblyBegin (xu);
  VecAssemblyEnd (xu);
  VecAssemblyBegin (xv);
  VecAssemblyEnd (xv);
  VecAssemblyBegin (xw);
  VecAssemblyEnd (xw);
  VecAssemblyBegin (xp);
  VecAssemblyEnd (xp);
  VecAssemblyBegin (xT);
  VecAssemblyEnd (xT);
  VecAssemblyBegin (xs);
  VecAssemblyEnd (xs);

  *curtime = t;

}
