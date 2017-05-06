#pragma once
#include "../Core/LinearCombination3.h"
#include "../Core/RealSquareMatrices.h"

namespace cagd
{
    class CyclicCurve3: public LinearCombination3
    {
    protected:
        GLuint              _n;                 //order
        GLdouble            _c_n;               //normalizing constants
        GLdouble            _lambda_n;          //phase change

        RealSquareMatrix    _bc;                //binomial coefficients

        GLdouble            _CalculateNormalizingCoefficient(GLuint n);

        GLvoid              _CalculateBinomialCoefficients(GLuint m, RealSquareMatrix &bc);
    public:
        //special constructor
        CyclicCurve3(GLuint n, GLenum data_usage_flag = GL_STATIC_DRAW);

        //redeclare and define inherited pure virutal methods from LinearCombination3
        GLboolean BlendingFunctionValues(GLdouble u, RowMatrix<GLdouble>& values) const;
        GLboolean CalculateDerivatives(GLuint max_order_of_derivates, GLdouble u, Derivatives& d) const;

        //redeclare and define inherited pure virutal methods from LocalCurveEnergy
        virtual CyclicCurve3* Clone();
    };
}
