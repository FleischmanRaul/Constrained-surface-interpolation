#include "TrigonometricBernsteinSurfaces.h"
#include <cmath>

using namespace std;

namespace cagd
{
    GLvoid TrigonometricBernsteinSurface3::_CalculateBinomialCoefficients(GLuint order, TriangularMatrix<GLdouble> &bc)
    {
        bc.ResizeRows(order + 1);
        bc(0, 0) = 1.0;

        for(GLuint r = 1; r<=order; ++r)
        {
            bc(r, 0) = 1.0;
            bc(r, r) = 1.0;

            for(GLuint i = 1; i<=r/2; ++i)
            {
                bc(r, i) = bc(r-1, i-1) + bc(r-1, i);
                bc(r, r-i) = bc(r, i);
            }
        }
    }

    GLboolean TrigonometricBernsteinSurface3::_CalculateNormalizingCoefficients(GLuint order, GLdouble alpha, RowMatrix<GLdouble> &c)
    {
        if(!order)
        {
            c.ResizeColumns(0);
            return GL_FALSE;
        }
        GLuint size = 2 * order + 1;
        c.ResizeColumns(size);

        GLdouble sa = pow(sin(alpha / 2.0), (GLint)(2*order));
        GLdouble ca = 2.0 * cos(alpha / 2.0);

    for(GLuint i=0; i <= order ; ++i)
    {
        c[i] = 0.0;
        for (GLuint r = 0; r <= i / 2 ; ++r)
        {
            c[i] += _bc( order , i - r ) * _bc( i - r , r ) * pow( ca , (GLint)( i - 2 * r ));
        }
    c[i] /= sa;
    c[size - 1 - i] = c[i];
    }
    return GL_TRUE;
    }
    //atirni az m-et orderre?
    TrigonometricBernsteinSurface3::TrigonometricBernsteinSurface3(GLdouble alpha, GLuint n, GLdouble beta, GLuint m):TensorProductSurface3(0.0, alpha, 0.0, beta, 2*n+1, 2*m+1),
    _alpha(alpha), _beta(beta), _n(n), _m(m)
    {
        GLuint max_order = (_n > -m ? _n : _m);
        _CalculateBinomialCoefficients(max_order, _bc);
        _CalculateNormalizingCoefficients(_n, _alpha, _u_c);
        _CalculateNormalizingCoefficients(_m, _beta, _v_c);
    }

    GLboolean TrigonometricBernsteinSurface3::UBlendingFunctionValues(GLdouble u, RowMatrix<GLdouble> &values) const
    {
        if(u < _u_min || u > _u_max)
        {
            values.ResizeColumns(0);
            return GL_FALSE;
        }

        GLuint size = 2 * _n + 1;
        values.ResizeColumns(size);

        for(GLuint i = 0; i<size; ++i)
        {
            values[i] = 0.0;
        }

        if(u == _u_max)
        {
            values[0] = 1.0;
        }
        else
        {
            if(u == _u_max)
            {
                values[size - 1] = 1.0;
            }
            else
            {
                GLdouble sau = sin((_alpha - u) / 2.0), su = sin(u/2.0);
                GLdouble factor = su / sau;
                GLdouble sau_order = pow(sau, (GLint)(size - 1));

                values[0] = sau_order;

                for(GLuint i = 1; i < size; ++i)
                {
                    values[i] = values[i-1] * factor;
                }
                for(GLuint i = 1; i < size; ++i)
                {
                    values[i] *= _u_c[i];
                }
            }
        }
        return GL_TRUE;
    }

    GLboolean TrigonometricBernsteinSurface3::VBlendingFunctionValues(GLdouble v, RowMatrix<GLdouble> &values) const
    {
        if(v < _v_min || v > _v_max)
        {
            values.ResizeColumns(0);
            return GL_FALSE;
        }

        GLuint size = 2 * _m + 1;
        values.ResizeColumns(size);

        for(GLuint i = 0; i<size; ++i)
        {
            values[i] = 0.0;
        }

        if(v == _v_max)
        {
            values[0] = 1.0;
        }
        else
        {
            if(v == _v_max)
            {
                values[size - 1] = 1.0;
            }
            else
            {
                GLdouble sbv = sin((_beta - v) / 2.0), sv = sin(v/2.0);
                GLdouble factor = sv / sbv;
                GLdouble sbv_order = pow(sbv, (GLint)(size - 1));

                values[0] = sbv_order;

                for(GLuint j = 1; j < size; ++j)
                {
                    values[j] = values[j-1] * factor;
                }
                for(GLuint j = 1; j < size; ++j)
                {
                    values[j] *= _v_c[j];
                }
            }
        }
        return GL_TRUE;
    }

    GLboolean TrigonometricBernsteinSurface3::CalculatePartialDerivatives(GLdouble u, GLdouble v, PartialDerivatives &pd) const
    {
        std::cout << "TrigonometricBernsteinSurface3::CalculatePartialDerivatives" << std::endl;
        pd.LoadNullVectors();
        if(u < _u_min || u > _u_max || v < _v_min || v > _v_max)
        {
            return GL_FALSE;
        }
        GLuint u_size = 2 * _n + 1, v_size = 2 * _m + 1;
        RowMatrix<GLdouble> Au, Av;

        if(!UBlendingFunctionValues(u, Au) || !VBlendingFunctionValues(v, Av))
        {
            return GL_FALSE;
        }
        GLdouble sua = 2.0 * sin(_alpha / 2.0), tua = tan(_alpha / 2.0);
        GLdouble sva = 2.0 * sin(_beta / 2.0), tva = tan(_beta / 2.0);

        for( GLuint i = 0; i < u_size; ++i)
        {
            DCoordinate3 diff_00, diff_01;
            for(GLuint j = 0; j < v_size; j++)
            {
                diff_00 += _data(i, j) * Av[j];
                diff_01 += _data(i, j) * (((j>0) ? (_v_c[j]/_v_c[j-1] * (GLdouble)j / sva * Av[j-1]): 0.0) - ((GLdouble)_m - (GLdouble)j) / tva * Av[j] -
                                          ((j < 2 * _m) ? (_v_c[j+1] * (2.0 * (GLdouble)_m - (GLdouble)j)/sva)*Av[j+1]:0.0));
            }
            pd(0, 0) += diff_00 * Au[i];
            pd(1, 0) += diff_01 * Au[i];

            pd(1, 1) += diff_00 * (((i>0) ? (_u_c[i] /_u_c[i-1] * (GLdouble)i / sua * Au[i-1]) : 0.0 ) -
                         ((GLdouble)_n - (GLdouble)i) / tua * Au[i]
                         -((i<2 * _n) ? (_u_c[i] /_u_c[i+1] * ( 2.0 * (GLdouble)_n - (GLdouble)i ) / sua ) * Au[i+1]:0.0));
        }
        return GL_TRUE;
    }
}



