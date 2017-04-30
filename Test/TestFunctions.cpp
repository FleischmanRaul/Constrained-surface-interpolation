#include <cmath>
#include "TestFunctions.h"
#include "Core/Constants.h"

using namespace cagd;
using namespace std;

GLdouble spiral_on_cone::u_min = -2 * PI;
GLdouble spiral_on_cone::u_max = 2 * PI;

DCoordinate3 spiral_on_cone::d0(GLdouble u)
{
    return DCoordinate3(u * cos(u), u * sin(u), u);
}

DCoordinate3 spiral_on_cone::d1(GLdouble u)
{
    GLdouble c = cos(u), s = sin(u);
    return DCoordinate3(c - u * s, s + u * c, 1.0);
}

DCoordinate3 spiral_on_cone::d2(GLdouble u)
{
    GLdouble c = cos(u), s = sin(u);
    return DCoordinate3(-2.0 * s - u * c, 2.0 * c - u * s, 0);
}

//DCoordinate3 spiral_on_cone::d0(GLdouble u)
//{
//    return DCoordinate3(cos(u),sin(u), 0);
//}

//DCoordinate3 spiral_on_cone::d1(GLdouble u)
//{
//    GLdouble c = cos(u), s = sin(u);
//    return DCoordinate3(-s, c, 0);
//}

//DCoordinate3 spiral_on_cone::d2(GLdouble u)
//{
//    GLdouble c = cos(u), s = sin(u);
//    return DCoordinate3(-c, -s, 0);
//}

//DCoordinate3 spiral_on_cone::d0(GLdouble u)
//{
//    return DCoordinate3((2+cos(2*u/3))*cos(u),
//                        (2+cos(2*u/3))*sin(u),
//                        sin(2*u/3));
//}

//DCoordinate3 spiral_on_cone::d1(GLdouble u)
//{
//    GLdouble c = cos(u), s = sin(u), p = 2*u/3;
//    return DCoordinate3(-(2.0/3)*c*sin(p)-(2+cos(p))*s,
//                        -2.0/3*sin(p)*s + (2 + cos(p)) * c,
//                        2.0/3*cos(p));
//}

//DCoordinate3 spiral_on_cone::d2(GLdouble u)
//{
//    GLdouble c = cos(u), s = sin(u), p = 2*u/3;
//    return DCoordinate3(4.0/3*sin(p)*s-4.0/9*cos(p)*c-(cos(p)+2)*c,
//                        -1.0/18 * (-sin(u/3)-25*sin(5*u/3)),
//                        -4.0/9*sin(p));
//}
