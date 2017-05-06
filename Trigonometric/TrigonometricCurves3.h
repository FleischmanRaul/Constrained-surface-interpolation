#pragma once
#include "../Core/Matrices.h"
#include "../Core/LinearCombination3.h"
#include "../Core/RealSquareMatrices.h"

namespace cagd
{
    class TrigonometricCurve3: public LinearCombination3
    {
    protected:
        GLuint              _n;                 //order
        GLdouble            _alpha;             //shape parameter
        RowMatrix<GLdouble> _c_n;               //normalizing constants
        RealSquareMatrix    _bc;                //binomial coefficients

        GLvoid              _CalculateBinomialCoefficients(GLuint n, RealSquareMatrix &bc);
        GLboolean           _CalculateNormalizingCoefficients(GLuint order, GLdouble alpha, RowMatrix<GLdouble> &c);

    public:
        //special constructor
        TrigonometricCurve3(GLuint n, GLdouble alpha, GLenum data_usage_flag = GL_STATIC_DRAW);

        //redeclare and define inherited pure virutal methods from LinearCombination3
        GLboolean BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble>& values) const;
        GLboolean CalculateDerivatives(GLuint max_order_of_derivatives, GLdouble u, Derivatives& d) const;

        //redeclare and define inherited pure virutal methods from LocalCurveEnergy
        virtual TrigonometricCurve3* Clone();
    };
}
