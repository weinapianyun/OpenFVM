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

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "variables.h"
#include "vector.h"
#include "globals.h"
#include "fill.h"

double
VolumeTotal ()
{

  int i;

  int element;

  double vol;

  vol = 0.0;

  for (i = 0; i < nbelements; i++)
    {
      element = i;

      vol += elements[element].Vp;
    }

  return vol;

}

double
VolumeFilled ()
{

  int i;

  int element;

  double vol;

  VecGhostGetLocalForm (xs, &xsl);

  vol = 0.0;

  for (i = 0; i < nbelements; i++)
    {
      element = i;

      vol += V_GetCmp (&xsl, element) * elements[element].Vp;
    }

  VecGhostRestoreLocalForm (xs, &xsl);

  return vol;

}

double
VolumeEntered (double dt)
{

  int i;

  int face;

  double vol;

  vol = 0.0;

  for (i = 0; i < nbfaces; i++)
    {

      face = i;

      if (faces[face].bc == INLET)
	{
	  vol += -V_GetCmp (&uf, face) * faces[face].Aj * dt;
	}

    }

  return vol;

}
