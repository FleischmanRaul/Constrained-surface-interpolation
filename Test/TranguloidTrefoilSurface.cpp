#include <cmath>
#include "Core/Constants.h"
#include "TranguloidTrefoilSurface.h"

using namespace cagd;

GLdouble TranguloidTrefoilSurface::u_min = -PI;
GLdouble TranguloidTrefoilSurface::u_max = PI;
GLdouble TranguloidTrefoilSurface::v_min = -PI;
GLdouble TranguloidTrefoilSurface::v_max = PI;


DCoordinate3 TranguloidTrefoilSurface::d00(GLdouble u, GLdouble v)
{
    return DCoordinate3(2.0*sin(3.0*u)/(2.0+cos(v)), 2.0*(sin(u)+2.0*sin(2.0*u))/(2.0+cos(v+2.0*PI/3.0)), (cos(u)-2.0*cos(2.0*u))*(2.0+cos(v))*(2.0+cos(v+2.0*PI/3.0))/4.0);
}

DCoordinate3 TranguloidTrefoilSurface::d10(GLdouble u, GLdouble v)
{
    return DCoordinate3(6.0*cos(3.0*u)/(cos(v)+2.0), 2.0*(cos(u)+4.0*cos(2.0*u))/(sin(v+PI/6.0)-2.0), 1/4*(sin(u)-4.0*sin(2.0*u))*(sin(v+PI/6.0)-2.0)*(cos(v)+2.0) );
}

DCoordinate3 TranguloidTrefoilSurface::d01(GLdouble u, GLdouble v)
{

    return DCoordinate3((2.0*sin(3.0*u)*sin(v))/((cos(v)+2.0)*(cos(v)+2.0)),  2.0*(sin(u)+2.0*sin(2.0*u))*cos(v+PI/6.0)/((sin(v+PI/6.0)-2.0)*(sin(v+PI/6.0)-2.0)), -1/4*(cos(u)-2*cos(2*u))*((cos(v)+2)*cos(v+PI/6)-sin(v)*(sin(v-PI/6)-2)));
}


