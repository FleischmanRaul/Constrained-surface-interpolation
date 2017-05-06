#pragma once

#include <GL/glew.h>
#include <QGLWidget>
#include <QGLFormat>
#include <QTimer>
#include "Parametric/ParametricCurves3.h"
#include "Parametric/ParametricSurfaces3.h"
#include "Cyclic/CyclicCurves3.h"
#include "Core/DCoordinates3.h"
#include "Core/Exceptions.h"
#include "Core/Matrices.h"
#include "Test/TestFunctions.h"
#include "Core/GenericCurves3.h"
#include "Core/Constants.h"
#include "Core/TriangulatedMeshes3.h"
#include "Core/TriangularFaces.h"
#include "Core/Lights.h"
#include "Core/HCoordinates3.h"
#include "Core/Colors4.h"
#include "Trigonometric/TrigonometricBernsteinSurfaces.h"
#include "Trigonometric/BicubicBezierPatches.h"
#include "Trigonometric/TrigonometricCurves3.h"
#include "GA/InterpolatingCurveBiology.h"


#include <stdlib.h>     //for using the function sleep
#include <unistd.h>
namespace cagd
{
class GLWidget: public QGLWidget
{
    Q_OBJECT

private:

    // variables defining the projection matrix
    float       _aspect;            // aspect ratio of the rendering window
    float       _fovy;              // field of view in direction y
    float       _z_near, _z_far;    // distance of near and far clipping planes

    // variables defining the model-view matrix
    float       _eye[3], _center[3], _up[3];

    // variables needed by transformations
    int         _angle_x, _angle_y, _angle_z;
    double      _zoom;
    double      _trans_x, _trans_y, _trans_z;

    bool        _evolve;

    QTimer      *_curveTimer;
    GLuint      _curvePopulationGeneration;

    // your other declarations
    //testing parametric curves
    cagd::ParametricCurve3* _pc;
    cagd::GenericCurve3* _image_of_pc;

    //testing generic curves
    GLuint                     _n, _mod, _div;
    RowMatrix<GLdouble>        _u;
    cagd::CyclicCurve3*        _cc;
    cagd::TrigonometricCurve3* _tc;
    cagd::GenericCurve3*       _icc;
    cagd::GenericCurve3*       _icc2;
    cagd::GenericCurve3*       _icc3;

    //testing CurveGA
    GLuint                  curve_iterations = 0;
    GLuint                  curve_generations = 0;
    GLuint                  curve_subgenerations = 0;
    GLdouble                mutation_radius;
    GLdouble                new_fitness;
    GLdouble                previous_fitness;
    cagd::CurveIndividual*  ci;
    cagd::CurvePopulation*  cp;
    ColumnMatrix<DCoordinate3> datatoint;
    ColumnMatrix<DCoordinate3> controllPoints;
    //testing triangulated meshes
    QTimer *_timer;
    GLfloat _angle;
    cagd::TriangulatedMesh3 _model;
    //creating light source
    DirectionalLight *dl;

    TriangularMatrix<ParametricSurface3::PartialDerivative> _pd;
    GLdouble _u_min, u_max, _v_min, _v_max;

    cagd::ParametricSurface3 * _ps;

    GLuint _u_div_point_count, _v_div_point_count;
    cagd::TriangulatedMesh3 * _ips;

    cagd::BicubicBezierPatch _patch;
    cagd::TriangulatedMesh3 *_before_interpolation, *_after_interpolation;

public:
    // special and default constructor
    // the format specifies the properties of the rendering window
    GLWidget(QWidget* parent = 0, const QGLFormat& format = QGL::Rgba | QGL::DepthBuffer | QGL::DoubleBuffer);

    // redeclared virtual functions
    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);

public slots:
    // public event handling methods/slots
    void set_angle_x(int value);
    void set_angle_y(int value);
    void set_angle_z(int value);

    void set_zoom_factor(double value);

    void set_trans_x(double value);
    void set_trans_y(double value);
    void set_trans_z(double value);
    void set_evolve_state(bool value);

private slots:
    void _evolveCurves();
};
}
