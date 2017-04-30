#pragma once

#include "../Core/DCoordinates3.h"
#include "../Core/Constants.h"

namespace cagd
{
    namespace TorusSurface
    {
        extern GLdouble a, c;

        extern GLdouble u_min, u_max,
                        v_min, v_max;

        DCoordinate3 d00(GLdouble u, GLdouble v);
        DCoordinate3 d10(GLdouble u, GLdouble v);
        DCoordinate3 d01(GLdouble u, GLdouble v);
    }
}
