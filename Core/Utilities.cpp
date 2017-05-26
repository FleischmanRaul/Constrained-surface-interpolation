#include "Utilities.h"

#include <omp.h>
#include <GL/glew.h>

namespace cagd
{

    Color4 coldToHotColormap(GLfloat value, GLfloat min_value, GLfloat max_value)
    {
        Color4 color(1.0, 1.0, 1.0);

        if (value < min_value)
        {
           value = min_value;
        }

        if (value > max_value)
        {
           value = max_value;
        }

        float dv = max_value - min_value;

        if (value < (min_value + 0.25f * dv))
        {
           color.r() = 0.0;
           color.g() = 4.0f * (value - min_value) / dv;
        }
        else
        {
            if (value < (min_value + 0.5f * dv))
            {
               color.r() = 0.0f;
               color.b() = 1.0f + 4.0f * (min_value + 0.25f * dv - value) / dv;
            }
            else
            {
                if (value < (min_value + 0.75f * dv))
                {
                   color.r() = 4.0f * (value - min_value - 0.5f * dv) / dv;
                   color.b() = 0.0f;
                }
                else
                {
                   color.g() = 1.0f + 4.0f * (min_value + 0.75f * dv - value) / dv;
                   color.b() = 0.0f;
                }
            }
        }

        return color;
    }
}
