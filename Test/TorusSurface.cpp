#include <cmath>
#include "Core/Constants.h"
#include "TorusSurface.h"

using namespace cagd;

GLdouble TorusSurface::u_min = 0.0;
GLdouble TorusSurface::u_max = TWO_PI;
GLdouble TorusSurface::v_min = 0.0;
GLdouble TorusSurface::v_max = TWO_PI;

GLdouble TorusSurface::a = 1;
GLdouble TorusSurface::c = 0.65;

DCoordinate3 TorusSurface::d00(GLdouble u, GLdouble v)
{
    return DCoordinate3((a+c*sin(u))*cos(v), (a+c*sin(u))*sin(v), c*cos(u));
}

DCoordinate3 TorusSurface::d10(GLdouble u, GLdouble v)
{
    return DCoordinate3(-c*cos(u)*cos(v), -c*cos(u)*sin(v), c*sin(u));
}

DCoordinate3 TorusSurface::d01(GLdouble u, GLdouble v)
{
    return DCoordinate3(-(a+c*sin(u))*sin(v), (a+c*sin(u))*cos(v), 0.0);
}
