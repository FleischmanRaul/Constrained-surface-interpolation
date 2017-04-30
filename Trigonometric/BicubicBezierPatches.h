#ifndef BICUBICBEZIERPATCHES_H
#define BICUBICBEZIERPATCHES_H

#include "Core/Matrices.h"
#include "Core/RealSquareMatrices.h"
#include "Core/TensorProductSurfaces3.h"

namespace cagd
{
    class BicubicBezierPatch: public TensorProductSurface3
    {
    public:
        BicubicBezierPatch();

        GLboolean UBlendingFunctionValues(GLdouble u_knot, RowMatrix<GLdouble>& blending_values) const;

        GLboolean VBlendingFunctionValues(GLdouble v_knot, RowMatrix<GLdouble>& blending_values) const;

        GLboolean CalculatePartialDerivatives(GLdouble u, GLdouble v, PartialDerivatives& pd) const;
    };
}


#endif // BICUBICBEZIERPATCHES_H
