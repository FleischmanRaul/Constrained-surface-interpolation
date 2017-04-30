#pragma once

#include "../Core/DCoordinates3.h"

namespace cagd
{
    namespace spiral_on_cone
    {
        extern GLdouble u_min, u_max;

        DCoordinate3 d0(GLdouble);
        DCoordinate3 d1(GLdouble);
        DCoordinate3 d2(GLdouble);
    }

    namespace SurfaceName
    {
        // possible shape parameters
        extern GLdouble parameter_1, parameter_2, parameter_3, parameter_4;

        // definition domain
        extern GLdouble u_min, u_max,
        v_min, v_max;

        DCoordinate3 d00(GLdouble u, GLdouble v); // zeroth order partial derivative, i.e. surface point
        DCoordinate3 d10(GLdouble u, GLdouble v); // first order partial derivative in direction u
        DCoordinate3 d01(GLdouble u, GLdouble v); // first order partial derivative in direction v
    }
}
