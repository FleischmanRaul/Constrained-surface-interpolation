#pragma once

#include "Core/Matrices.h"
#include "Core/RealSquareMatrices.h"
#include "Core/TensorProductSurfaces3.h"

namespace cagd
{
    class CyclicSurface3: public TensorProductSurface3
    {
    protected:
        GLuint              _n, _m;
        GLdouble            _lambda_n, _lambda_m;          //phase change
        RowMatrix<GLdouble> _u_c, _v_c;

        TriangularMatrix<GLdouble> _bc;

        GLvoid    _CalculateBinomialCoefficients(GLuint order, TriangularMatrix<GLdouble> &bc);
        GLboolean _CalculateNormalizingCoefficients(GLuint order, RowMatrix<GLdouble> &c);

    public:
        CyclicSurface3(GLuint n, GLuint m);

        GLboolean UBlendingFunctionValues(GLdouble u, RowMatrix<GLdouble>& values) const;

        GLboolean VBlendingFunctionValues(GLdouble v, RowMatrix<GLdouble>& values) const;

        GLboolean CalculatePartialDerivatives(GLuint maximum_order_of_partial_derivatives, GLdouble u, GLdouble v, PartialDerivatives &pd) const;

        virtual CyclicSurface3* Clone();
    };
}
