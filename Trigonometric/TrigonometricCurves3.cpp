#include "Trigonometric/TrigonometricCurves3.h"
#include "Core/Constants.h"
#include <iostream>
#include <cmath>

using namespace std;

namespace cagd
{
    GLvoid TrigonometricCurve3::_CalculateBinomialCoefficients(GLuint order, RealSquareMatrix &bc)
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

    GLboolean TrigonometricCurve3::_CalculateNormalizingCoefficients(GLuint order, GLdouble alpha, RowMatrix<GLdouble> &c)
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


    TrigonometricCurve3::TrigonometricCurve3(GLuint n,GLdouble alpha, GLenum data_usage_flag):
        LinearCombination3(0.0, alpha, 2 * n + 1, data_usage_flag),
        _n(n),
        _alpha(alpha)
    {
        _CalculateBinomialCoefficients(_n, _bc);
        _CalculateNormalizingCoefficients(_n, _alpha, _c_n);
    }

    GLboolean TrigonometricCurve3::BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble>& values) const
    {
        if(u < 0.0 || u > _alpha)
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

        if(u == 0.0)
        {
            values[0] = 1.0;
        }
        else
        {
            if(u == _alpha)
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
                    values[i] *= _c_n[i];
                }
            }
        }
        return GL_TRUE;
    }
    GLboolean TrigonometricCurve3::CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives& d) const
    {
        if(u < 0.0 || u > _alpha)
        {
            return GL_FALSE;
        }
        d.ResizeRows(max_order_of_derivatives + 1);
        d.LoadNullVectors();

        GLuint u_size = 2 * _n + 1;
        RowMatrix<GLdouble> Au;

        if(!BlendingFunctionValues(u, Au))
        {
            return GL_FALSE;
        }
        GLdouble sua = 2.0 * sin(_alpha / 2.0), tua = tan(_alpha / 2.0);

        Matrix<GLdouble> dAu(max_order_of_derivatives+1, u_size);

        dAu.SetRow(0, Au);
        for( GLuint i = 0; i < u_size; ++i)
        {
           d[0] += dAu(0, i) * _data[i];
        }

        for(GLuint r = 1; r<=max_order_of_derivatives; r++)
        {
            for( GLuint i = 0; i < u_size; ++i)
            {
               dAu(r, i) = ((i>0) ? _c_n[i] /_c_n[i-1] * (GLdouble)i / sua * dAu(r-1,i-1) : 0.0 ) -
                       (_n - (GLdouble)i) / tua * dAu(r-1,i)
                       -((i<2 * _n) ? (_c_n[i] /_c_n[i+1] * ( 2.0 * _n - (GLdouble)i ) / sua ) * Au[i+1]:0.0);

               d[r] += dAu(r, i) * _data[i];
            }
        }

        return GL_TRUE;
    }
    TrigonometricCurve3* TrigonometricCurve3::Clone()
    {
        return new TrigonometricCurve3(*this);
    }
}
