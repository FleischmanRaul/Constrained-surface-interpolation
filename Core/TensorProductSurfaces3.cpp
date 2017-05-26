#include "TensorProductSurfaces3.h"
#include "RealSquareMatrices.h"
#include "Core/TriangulatedMeshes3.h"
#include "Utilities.h"
using namespace cagd;
using namespace std;

// done: special constructor
TensorProductSurface3::TensorProductSurface3(
        GLdouble u_min, GLdouble u_max,
        GLdouble v_min, GLdouble v_max,
        GLuint row_count, GLuint column_count,
        GLboolean is_closed):
        _is_closed(is_closed),
        _vbo_data(0),
        _u_min(u_min), _u_max(u_max),
        _v_min(v_min), _v_max(v_max),
        _data(row_count, column_count)
{
}

// done: copy constructor
TensorProductSurface3::TensorProductSurface3(const TensorProductSurface3& surface):
        _is_closed(surface._is_closed),
        _vbo_data(0),
        _u_min(surface._u_min), _u_max(surface._u_max),
        _v_min(surface._v_min), _v_max(surface._v_max),
        _data(surface._data)
{
    if (surface._vbo_data)
    {
        UpdateVertexBufferObjectsOfData();
    }
}


TriangulatedMesh3* TensorProductSurface3::GenerateImage(GLuint u_div_point_count, GLuint v_div_point_count, enum TensorProductSurface3::ImageColorScheme color_scheme, GLenum usage_flag) const
{
    if (u_div_point_count <= 1 || v_div_point_count <= 1)
        return GL_FALSE;

    // calculating number of vertices, unit normal vectors and texture coordinates
    GLuint vertex_count = u_div_point_count * v_div_point_count;

    // calculating number of triangular faces
    GLuint face_count = 2 * (u_div_point_count - 1) * (v_div_point_count - 1);

    TriangulatedMesh3 *result = nullptr;
    result = new TriangulatedMesh3(vertex_count, face_count, usage_flag);

    if (!result)
        return nullptr;

    // uniform subdivision grid in the definition domain
    GLdouble du = (_u_max - _u_min) / (u_div_point_count - 1);
    GLdouble dv = (_v_max - _v_min) / (v_div_point_count - 1);

    // uniform subdivision grid in the unit square
    GLfloat sdu = 1.0f / (u_div_point_count - 1);
    GLfloat tdv = 1.0f / (v_div_point_count - 1);

    GLuint maximum_order_of_partial_derivatives;

    switch (color_scheme)
    {
    case DEFAULT_NULL_FRAGMENT:
    case NORMAL_LENGTH_FRAGMENT:
        maximum_order_of_partial_derivatives = 1;
        break;

    case WILLMORE_ENERGY_FRAGMENT:
    case LOGARITHMIC_WILLMORE_ENERGY_FRAGMENT:
    case UMBILIC_DEVIATION_ENERGY_FRAGMENT:
    case LOGARITHMIC_UMBILIC_DEVIATION_ENERGY_FRAGMENT:
    case TOTAL_CURVATURE_ENERGY_FRAGMENT:
    case LOGARITHMIC_TOTAL_CURVATURE_ENERGY_FRAGMENT:
        maximum_order_of_partial_derivatives = 2;
        break;
    }

    std::vector<GLdouble> fragments(vertex_count);

    // for face indexing
    GLuint current_face = 0;

    // partial derivatives of order 0, 1, 2, and 3
    PartialDerivatives pd(2);

    for (GLuint i = 0; i < u_div_point_count; ++i)
    {
        GLdouble u = _u_min + i * du;
        GLfloat  s = i * sdu;
        for (GLuint j = 0; j < v_div_point_count; ++j)
        {
            GLdouble v = _v_min + j * dv;
            GLfloat  t = j * tdv;

            /*
                    3-2
                    |/|
                    0-1
                */
            GLuint index[4];

            index[0] = i * v_div_point_count + j;
            index[1] = index[0] + 1;
            index[2] = index[1] + v_div_point_count;
            index[3] = index[2] - 1;

            // calculating all needed surface data
            CalculatePartialDerivatives(maximum_order_of_partial_derivatives, u, v, pd);

            // surface point
            (*result)._vertex[index[0]] = pd(0, 0);

            // unit surface normal
            (*result)._normal[index[0]] = pd(1, 0);
            (*result)._normal[index[0]] ^= pd(1, 1);

            GLdouble normal_length = (*result)._normal[index[0]].length();

            if (normal_length && normal_length != 1.0)
            {
                (*result)._normal[index[0]] /= normal_length;
            }

            // coefficients of first and second fundamental forms
            GLdouble e_1 = 0.0, f_1 = 0.0, g_1 = 0.0;
            GLdouble e_2 = 0.0, f_2 = 0.0, g_2 = 0.0;

            // possible Gaussian curvature
            GLdouble K = 0.0;

            // possible mean curvature
            GLdouble H = 0.0;

            if (color_scheme >= WILLMORE_ENERGY_FRAGMENT)
            {
                e_1 = pd(1,0) * pd(1,0);
                f_1 = pd(1,0) * pd(1,1);
                g_1 = pd(1,1) * pd(1,1);

                e_2 = (*result)._normal[index[0]] * pd(2,0);
                f_2 = (*result)._normal[index[0]] * pd(2,1);
                g_2 = (*result)._normal[index[0]] * pd(2,2);
            }

            if (color_scheme >= WILLMORE_ENERGY_FRAGMENT)
            {
                H = (g_1 * e_2 - 2.0 * f_1 * f_2 + e_1 * g_2) / (e_1 * g_1 - f_1 * f_1);
            }

            if (color_scheme >= UMBILIC_DEVIATION_ENERGY_FRAGMENT)
            {
                K = (e_2 * g_2 - f_2 * f_2) / (e_1 * g_1 - f_1 * f_1);
            }

            switch (color_scheme)
            {
            case DEFAULT_NULL_FRAGMENT:
                fragments[index[0]] = 0.0;
                break;

            case NORMAL_LENGTH_FRAGMENT:
                fragments[index[0]] = normal_length;
                break;

            case WILLMORE_ENERGY_FRAGMENT:
                fragments[index[0]] = H * H * normal_length;
                break;

            case LOGARITHMIC_WILLMORE_ENERGY_FRAGMENT:
                fragments[index[0]] = log(1.0 + H * H) * normal_length;
                break;

            case UMBILIC_DEVIATION_ENERGY_FRAGMENT:
                fragments[index[0]] = 4.0 * (H * H - K) * normal_length;
                break;

            case LOGARITHMIC_UMBILIC_DEVIATION_ENERGY_FRAGMENT:
                fragments[index[0]] = log(1.0 + 4.0 * (H * H - K)) * normal_length;
                break;

            case TOTAL_CURVATURE_ENERGY_FRAGMENT:
                fragments[index[0]] = (1.5 * H * H - 0.5 * K) * normal_length;
                break;

            case LOGARITHMIC_TOTAL_CURVATURE_ENERGY_FRAGMENT:
                fragments[index[0]] = log(1.0 + (1.5 * H * H - 0.5 * K)) * normal_length;
                break;
            }

            // texture coordinates
            (*result)._tex[index[0]].s() = s;
            (*result)._tex[index[0]].t() = t;

            // faces
            if (i < u_div_point_count - 1 && j < v_div_point_count - 1)
            {
                (*result)._face[current_face][0] = index[0];
                (*result)._face[current_face][1] = index[1];
                (*result)._face[current_face][2] = index[2];
                ++current_face;

                (*result)._face[current_face][0] = index[0];
                (*result)._face[current_face][1] = index[2];
                (*result)._face[current_face][2] = index[3];
                ++current_face;
            }
        }
    }

    if (color_scheme > DEFAULT_NULL_FRAGMENT)
    {
        GLdouble min_fragment =  numeric_limits<GLdouble>::max();
        GLdouble max_fragment = -numeric_limits<GLdouble>::max();

        for (GLuint i = 0; i < fragments.size(); i++)
        {
            if (min_fragment > fragments[i])
            {
                min_fragment = fragments[i];
            }

            if (max_fragment < fragments[i])
            {
                max_fragment = fragments[i];
            }
        }

        for (GLuint i = 0; i < vertex_count; i++)
        {
            (*result)._color[i] = coldToHotColormap(fragments[i], min_fragment, max_fragment);
        }
    }


    return result;
}

// ensures interpolation, i.e. s(u_i, v_j) = d_{ij}
GLboolean TensorProductSurface3::UpdateDataForInterpolation(const RowMatrix<GLdouble>& u_knot_vector, const ColumnMatrix<GLdouble>& v_knot_vector, Matrix<DCoordinate3>& data_points_to_interpolate)
{
    GLuint row_count = _data.GetRowCount();
    if (!row_count)
        return GL_FALSE;
    GLuint column_count = _data.GetColumnCount();
    if (!column_count)
        return GL_FALSE;

    if (u_knot_vector.GetColumnCount() != row_count || v_knot_vector.GetRowCount() != column_count || data_points_to_interpolate.GetRowCount() != row_count || data_points_to_interpolate.GetColumnCount() != column_count)
        return GL_FALSE;

    // 1: calculate the u-collocation matrix and perfom LU-decomposition on it
    RowMatrix<GLdouble> u_blending_values;

    RealSquareMatrix u_collocation_matrix(row_count);

    for (GLuint i = 0; i < row_count; ++i)
    {
        if (!UBlendingFunctionValues(u_knot_vector(i), u_blending_values))
            return GL_FALSE;
        u_collocation_matrix.SetRow(i, u_blending_values);
    }
    if (!u_collocation_matrix.PerformLUDecomposition())
        return GL_FALSE;

    // 2: calculate the v-collocation matrix and perform LU-decomposition on it
    RowMatrix<GLdouble> v_blending_values;

    RealSquareMatrix v_collocation_matrix(column_count);

    for (GLuint j = 0; j < column_count; ++j)
    {
        if (!VBlendingFunctionValues(v_knot_vector(j), v_blending_values))
            return GL_FALSE;
        v_collocation_matrix.SetRow(j, v_blending_values);
    }

    if (!v_collocation_matrix.PerformLUDecomposition())
            return GL_FALSE;

    // 3:   for all fixed j in {0, 1,..., column_count} determine control points
    //
    //      a_k(v_j) = sum_{l=0}^{column_count} _data(l, j) G_l(v_j), k = 0, 1,..., row_count
    //
    //      such that
    //
    //      sum_{k=0}^{row_count} a_k(v_j) F_k(u_i) = data_points_to_interpolate(i, j),
    //
    //      for all i = 0, 1,..., row_count.
    Matrix<DCoordinate3> a(row_count, column_count);
    if (!u_collocation_matrix.SolveLinearSystem(data_points_to_interpolate, a))
        return GL_FALSE;

    // 4:   for all fixed i in {0, 1,..., row_count} determine control point
    //
    //      _data[i][j], j = 0, 1,..., column_count
    //
    //      such that
    //
    //      sum_{l=0}^{column_count} _data(i, l) G_l(v_j) = a_i(v_j)
    //
    //      for all j = 0, 1,..., column_count.
    if (!v_collocation_matrix.SolveLinearSystem(a, _data, GL_FALSE))
        return GL_FALSE;

    return GL_TRUE;
}

// done: assignment operator
TensorProductSurface3& TensorProductSurface3::operator =(const TensorProductSurface3& surface)
{
    if(this != &surface)
    {
        DeleteVertexBufferObjectsOfData();

        _is_closed = surface._is_closed;
        _vbo_data = surface._vbo_data;
        _u_min = surface._u_min;
        _u_max = surface._u_max;
        _v_min = surface._v_min;
        _v_max = surface._v_max;
        _data = surface._data;

        // !
        if (surface._vbo_data)
        {
            UpdateVertexBufferObjectsOfData();
        }
    }
    return *this;
}

// done: set/get the definition domain of the surface
GLvoid TensorProductSurface3::SetUInterval(GLdouble u_min, GLdouble u_max)
{
    _u_min = u_min;
    _u_max = u_max;
}
GLvoid TensorProductSurface3::SetVInterval(GLdouble v_min, GLdouble v_max)
{
    _v_min = v_min;
    _v_max = v_max;
}

GLvoid TensorProductSurface3::GetUInterval(GLdouble& u_min, GLdouble& u_max) const
{
    u_min = _u_min;
    u_max = _u_max;
}
GLvoid TensorProductSurface3::GetVInterval(GLdouble& v_min, GLdouble& v_max) const
{
    v_min = _v_min;
    v_max = _v_max;
}

GLvoid TensorProductSurface3::GetDefinitionDomain(GLdouble& u_min, GLdouble& u_max, GLdouble& v_min, GLdouble& v_max)const
{
    u_min = _u_min;
    u_max = _u_max;
    v_min = _v_min;
    v_max = _v_max;
}

// done: set coordinates of a selected data point
GLboolean TensorProductSurface3::SetData(GLuint row, GLuint column, GLdouble x, GLdouble y, GLdouble z)
{
    if(row < _data.GetRowCount() && column < _data.GetColumnCount())
    {
        DCoordinate3 &cp = _data(row,column);
        cp[0] = x;
        cp[1] = y;
        cp[2] = z;
        return GL_TRUE;
    }
    return GL_FALSE;
}

GLboolean TensorProductSurface3::SetData(GLuint row, GLuint column, const DCoordinate3& point)
{
    if(row < _data.GetRowCount() && column < _data.GetColumnCount())
    {
        _data(row,column) = point;
        return GL_TRUE;
    }
    return GL_FALSE;
}

// done: get coordinates of a selected data point
GLboolean TensorProductSurface3::GetData(GLuint row, GLuint column, GLdouble& x, GLdouble& y, GLdouble& z) const
{
    if(row < _data.GetRowCount() && column < _data.GetColumnCount())
    {
        const DCoordinate3& coord = _data(row, column);
        x = coord.x();
        y = coord.y();
        z = coord.z();
        return GL_TRUE;
    }
    return GL_FALSE;
}

GLboolean TensorProductSurface3::GetData(GLuint row, GLuint column, DCoordinate3& point) const
{
    if(row < _data.GetRowCount() && column < _data.GetColumnCount())
    {
        point = _data(row, column);
        return GL_TRUE;
    }
    return GL_FALSE;
}

// done: get data by value
DCoordinate3 TensorProductSurface3::operator ()(GLuint row, GLuint column) const
{
    return _data(row, column);
}

// done: get data by reference
DCoordinate3& TensorProductSurface3::operator ()(GLuint row, GLuint column)
{
    return _data(row, column);
}

// done: VBO handling methods (assume that Matrix<DCoordinate3> _data corresponds to a control net!)
GLvoid TensorProductSurface3::DeleteVertexBufferObjectsOfData()
{
    if (_vbo_data)
    {
        glDeleteBuffers(1, &_vbo_data);
        _vbo_data = 0;
    }
}
GLboolean TensorProductSurface3::RenderData(GLenum render_mode) const
{
    if (!_vbo_data)
        return GL_FALSE;
    if (render_mode != GL_LINE_STRIP && render_mode != GL_LINE_LOOP && render_mode != GL_POINTS)
    {
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableClientState(GL_VERTEX_ARRAY);
        return GL_FALSE;
    }

    // enable client states of the vertex array
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, _vbo_data);
    glVertexPointer(3, GL_FLOAT, 0, (const GLvoid*)0);

    GLuint first = 0;
    //we must render every row and every column
    for (GLuint i = 0; i < _data.GetRowCount(); ++i)
    {
        //glDrawArrays(GLenum mode, GLint first, GLsizei count);
        //render primitives from array data
        glDrawArrays(render_mode, first, _data.GetColumnCount());
        first += _data.GetColumnCount();
    }

    for (GLuint j = 0; j < _data.GetColumnCount(); ++j)
    {
        //render primitives from array data
        glDrawArrays(render_mode, first, _data.GetRowCount());
        first += _data.GetRowCount();
    }

    // unbind any buffer object previously bound and restore client memory usage
    // for these buffer object targets
    glDisableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    return GL_TRUE;
}


GLboolean TensorProductSurface3::UpdateVertexBufferObjectsOfData(GLenum usage_flag)
{
    if (usage_flag != GL_STREAM_DRAW  && usage_flag != GL_STREAM_READ  && usage_flag != GL_STREAM_COPY
     && usage_flag != GL_STATIC_DRAW  && usage_flag != GL_STATIC_READ  && usage_flag != GL_STATIC_COPY
     && usage_flag != GL_DYNAMIC_DRAW && usage_flag != GL_DYNAMIC_READ && usage_flag != GL_DYNAMIC_COPY)
        return GL_FALSE;

    // deleting old vertex buffer objects
    DeleteVertexBufferObjectsOfData();

    glGenBuffers(1, &_vbo_data);

    if (!_vbo_data)
        return GL_FALSE;

    // For efficiency reasons we convert all GLdouble coordinates
    // to GLfloat coordinates: we will use auxiliar pointers for
    // buffer data loading and functions glMapBuffer/glUnmapBuffer.

    // Notice that multiple buffers can be mapped simultaneously.

    GLuint data_count = _data.GetRowCount() * _data.GetColumnCount();
    GLuint data_byte_size = 6 * data_count * sizeof(GLfloat);

    glBindBuffer(GL_ARRAY_BUFFER, _vbo_data);
    glBufferData(GL_ARRAY_BUFFER, data_byte_size, 0, usage_flag);

    GLfloat *data_coordinate = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

    if (!data_coordinate)
    {
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        DeleteVertexBufferObjectsOfData();
        return GL_FALSE;
    }

    for (GLuint i = 0; i < _data.GetRowCount(); ++i)
        {
            for (GLuint j = 0; j < _data.GetColumnCount(); ++j)
            {
                DCoordinate3 &cp = _data(i, j);
                for (GLuint k = 0; k < 3; ++k)
                {
                    *data_coordinate = (GLfloat)cp[k];
                    ++data_coordinate;
                }
            }
        }

        for (GLuint j = 0; j < _data.GetColumnCount(); ++j)
        {
            for (GLuint i = 0; i < _data.GetRowCount(); ++i)
            {
                DCoordinate3 &cp = _data(i, j);
                for (GLuint k = 0; k < 3; ++k)
                {
                    *data_coordinate = (GLfloat)cp[k];
                    ++data_coordinate;
                }
            }
        }

        if (!glUnmapBuffer(GL_ARRAY_BUFFER))
        {
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            DeleteVertexBufferObjectsOfData();
            return GL_FALSE;
        }

        glBindBuffer(GL_ARRAY_BUFFER, 0);

        return GL_TRUE;
    }

// done: destructor
TensorProductSurface3::~TensorProductSurface3()
{
    DeleteVertexBufferObjectsOfData();
}


// special/default constructor
//PartialDerivatives(GLuint maximum_order_of_partial_derivatives = 0);

TensorProductSurface3::PartialDerivatives::PartialDerivatives(GLuint maximum_order_of_partial_derivatives):
    TriangularMatrix<DCoordinate3>(maximum_order_of_partial_derivatives + 1)
{
}

// when called, all inherited Descartes coordinates are set to the origin
GLvoid TensorProductSurface3::PartialDerivatives::LoadNullVectors()
{
    for (GLuint r = 0; r < _data.size(); r++)
    {
        for (GLuint c = 0; c <= r; c++)
        {
            DCoordinate3 &reference = _data[r][c];

            for (GLuint j = 0; j < 3; j++)
            {
                reference[j] = 0.0;
            }
        }
    }
}

GLdouble TensorProductSurface3::Fitness(GLuint u_div, GLuint v_div, SurfaceEnergyType etype)
{
    GLdouble result = 0.0;
    if(u_div % 2)           //u_div and v_div must be even numbers
    {
        ++u_div;
    }

    if(v_div % 2)
    {
        ++v_div;
    }

    GLdouble u_max, u_min, v_max, v_min;

    GetUInterval(u_min, u_max);
    GetVInterval(v_min, v_max);

    GLdouble du = (u_max - u_min) / u_div;
    GLdouble dv = (v_max - v_min) / v_div;

    RowMatrix<GLdouble> u(u_div + 1);

    for(GLuint i = 0; i < u_div; ++i)
    {
        u[i] = u_min + i * du;
    }
    u[u_div] = u_max;

    RowMatrix<GLdouble> v(v_div + 1);

    for(GLuint j = 0; j < u_div; ++j)
    {
        v[j] = v_min + j * dv;
    }
    v[v_div] = v_max;

    Matrix<GLdouble> phi(u_div+1, v_div+1);

    for(GLuint i = 0; i<= u_div; ++i)
    {
        for(GLuint j = 0; j<= v_div; ++j)
        {
            phi(i,j) = SurfaceFunctionals(u[i], v[j], etype);
        }
    }

    Matrix<GLdouble> weight(3,3);
    weight(0, 0) = weight(0, 2) = weight(2, 0) = weight(2, 2) = 1.0;
    weight(0, 1) = weight(1, 0) = weight(1, 2) = weight(2, 1) = 4.0;
    weight(1, 1) = 16.0;

    for(GLint i = 1; i <= (GLint)u_div; i+=2)
    {
        for(GLint j = 1; j <= (GLint)v_div; j+=2)
        {
            GLdouble S = 0.0;
            for(GLint k = -1.0; k <= 1; ++k)
            {
                for(GLint l = -1; l <= 1; ++l)
                {
                    S += weight(k+1, l+1) * phi(i+k,j+l);
                }
            }
            result += S;
        }
    }

    result *= du * dv / 9.0;
    if(result > 0)              //the energy must be negative
    {
        result *= -1;
    }
    return result;
}

GLdouble TensorProductSurface3::SurfaceFunctionals(GLdouble u, GLdouble v, const SurfaceEnergyType etype)
{
    GLdouble result;
    PartialDerivatives pd;
    GLuint order_of_derivatives_needed = 2;

    if(etype == MEHLUM_TARROU)
    {
        cout << "Not finished yet!" << endl;
        order_of_derivatives_needed = 3;
    }

    if(!CalculatePartialDerivatives(order_of_derivatives_needed,u,v,pd))
    {
        cout << "Could not calculate partial derivetives!" << endl;
        return 0.0;
    }

    GLdouble e_1 = pd(1,0) * pd(1,0);
    GLdouble f_1 = pd(1,0) * pd(1,1);
    GLdouble g_1 = pd(1,1) * pd(1,1);
    DCoordinate3 n0 = (pd(1,0) ^ pd(1,1));
    n0.normalize();
    GLdouble e_2 = n0 * pd(2,0);
    GLdouble f_2 = n0 * pd(2,1);
    GLdouble g_2 = n0 * pd(2,2);
    GLdouble H = (g_1 * e_2 - 2.0 * f_1 * f_2 + e_1 * g_2) / (e_1 * g_1 - f_1 * f_1);
    GLdouble K = (e_2 * g_2 - f_2 * f_2) / (e_1 * g_1 - f_1 * f_1);

    switch (etype)
    {
    case WILLMORE:
        result = H * H;
        break;
    case UMBILIC_DEVIATION:
        result = 4 * (H * H - K);
        break;
    case TOTAL_CURVATURE:
        result = 1.5 * H * H - 0.5 * K;
        break;
    case MEHLUM_TARROU:
        result = 0.0;
        break;
    default:
        cout << "Unknown functional type!" << endl;
        break;
    }
    return result;
}
