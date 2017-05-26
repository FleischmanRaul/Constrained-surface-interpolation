#include <iostream>
#include <cmath>
#include "CyclicSurfaces3.h"
#include "Core/Constants.h"

using namespace std;

namespace cagd
{
    GLvoid CyclicSurface3::_CalculateBinomialCoefficients(GLuint order, TriangularMatrix<GLdouble> &bc)
    {
        bc.ResizeRows(order + 1);
        bc(0, 0) = 1.0;

        for(GLuint r = 1; r <= order; ++r)
        {
            bc(r, 0) = 1.0;
            bc(r, r) = 1.0;

            for(GLuint i = 1; i <= r/2; ++i)
            {
                bc(r, i) = bc(r - 1, i - 1) + bc(r - 1, i);
                bc(r, r - i) = bc(r, i);
            }
        }
    }
    GLboolean CyclicSurface3::_CalculateNormalizingCoefficients(GLuint order, RowMatrix<GLdouble> &c)
    {
        if(!order)
        {
            c.ResizeColumns(0);
            return GL_FALSE;
        }
        GLuint size = 2 * order + 1;
        c.ResizeColumns(size);

        c[0] = 1.0 / 3.0;
        for(GLuint i = 1; i < size ; ++i)
        {
            c[i] = c[i-1] * (GLdouble)(i+1) / (GLdouble)(2*i + 3);
        }
        return GL_TRUE;
    }


    CyclicSurface3::CyclicSurface3(GLuint n, GLuint m):
        TensorProductSurface3(0.0,TWO_PI,0.0,TWO_PI,2*n+1,2*m+1), _n(n), _m(m),_lambda_n(TWO_PI / (2 * n + 1)),_lambda_m(TWO_PI / (2 * m + 1))
    {
        GLuint max_order = (_n > _m ? _n : _m);
        _CalculateBinomialCoefficients(max_order, _bc);
        _CalculateNormalizingCoefficients(_n, _u_c);
        _CalculateNormalizingCoefficients(_m, _v_c);
    }

    GLboolean CyclicSurface3::UBlendingFunctionValues(GLdouble u, RowMatrix<GLdouble>& values) const
    {
        if(u < 0.0 || u > TWO_PI)
        {
            values.ResizeColumns(0);
            return GL_FALSE;
        }
        GLuint size = 2 * _n + 1;
        values.ResizeColumns(size);

        for(GLuint i = 0; i < 2* _n + 1; ++i)
        {
            values[i] = _u_c[i] * pow(1.0 + cos(u - i * _lambda_n), (GLint)_n);
        }
        return GL_TRUE;
    }

    GLboolean CyclicSurface3::VBlendingFunctionValues(GLdouble v, RowMatrix<GLdouble>& values) const
    {
        if(v < 0.0 || v > TWO_PI)
        {
            values.ResizeColumns(0);
            return GL_FALSE;
        }
        GLuint size = 2 * _m + 1;
        values.ResizeColumns(size);

        for(GLuint i = 0; i < 2* _m + 1; ++i)
        {
            values[i] = _v_c[i] * pow(1.0 + cos(v - i * _lambda_m), (GLint)_m);
        }
        return GL_TRUE;
    }

    GLboolean CyclicSurface3::CalculatePartialDerivatives(GLuint maximum_order_of_partial_derivatives, GLdouble u, GLdouble v, PartialDerivatives &pd) const
    {
        if(u < _u_min || u > _u_max || v < _v_min || v > _v_max)
        {
            pd.ResizeRows(0);
            return GL_FALSE;
        }

        pd.ResizeRows(maximum_order_of_partial_derivatives + 1);
        pd.LoadNullVectors();

        DCoordinate3 u_centroid;
        for(GLuint i = 0; i <= 2 * _n; ++i)
        {
            u_centroid += _data(0,i);
        }
        u_centroid /= (GLdouble)(2 * _n + 1);

        DCoordinate3 v_centroid;
        for(GLuint i = 0; i <= 2 * _m; ++i)             //proly nem jo
        {
            v_centroid += _data(0,i);
        }
        v_centroid /= (GLdouble)(2 * _m + 1);
        GLuint u_size = 2 * _n + 1, v_size = 2 * _m + 1;
        Matrix<GLdouble> dAu(maximum_order_of_partial_derivatives + 1, u_size);
        Matrix<GLdouble> dAv(maximum_order_of_partial_derivatives + 1, v_size);

//        for(GLuint r = 0; r <= max_order_of_derivatives; ++r)
//        {
//            for(GLuint i = 0; i <= 2 * _n; ++i)
//            {
//                GLdouble sum_k = 0.0;

//                for(GLuint k = 0; k <= _n - 1; ++k)
//                {
//                    sum_k += pow(_n - k, (GLint)r) *
//                            _bc(2 * _n, k) *
//                            cos((_n - k) * (u - i * _lambda_n) + r * PI / 2.0);
//                }

//                d[r] += sum_k * _data[i];
//            }

//            d[r] *= 2.0;
//            d[r] /= (GLdouble)(2 * _n + 1);
//            d[r] /= _bc(2 * _n, _n);

//            if (r == 0)
//            {
//                d[r] += centroid;
//            }
//        }


        return GL_TRUE;
    }

    CyclicSurface3* CyclicSurface3::Clone()
    {
        return new CyclicSurface3(*this);
    }
}
