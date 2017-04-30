#include "GLWidget.h"

#include <iostream>
using namespace std;

#include <GL/glu.h>
#include <Core/Exceptions.h>

namespace cagd
{
//--------------------------------
// special and default constructor
//--------------------------------
GLWidget::GLWidget(QWidget *parent, const QGLFormat &format): QGLWidget(format, parent)
{
    _curveTimer = new QTimer(this);
    _curveTimer->setInterval(1);
    connect(_curveTimer, SIGNAL(timeout()), this, SLOT(_evolveCurves()));
}

//--------------------------------------------------------------------------------------
// this virtual function is called once before the first call to paintGL() or resizeGL()
//--------------------------------------------------------------------------------------
DirectionalLight * initializeDirectionalLight()
{
    DirectionalLight * dl = 0;

    HCoordinate3 direction(0.0,0.0,1.0,0.0);
    Color4      ambient   (0.2,0.4,0.4,1.0);
    Color4      diffuse   (0.8,0.8,0.8,1.0);
    Color4      specular   (1.0,1.0,1.0,1.0);

    dl = new DirectionalLight(GL_LIGHT0, direction, ambient, diffuse, specular);
    return dl;
}

//--------------------------------------------------------------------------------------
// this virtual function is called once before the first call to paintGL() or resizeGL()
//--------------------------------------------------------------------------------------
void GLWidget::initializeGL()
{
    // creating a perspective projection matrix
    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();

    _aspect = (float)width() / (float)height();
    _z_near = 1.0;
    _z_far = 1000.0;
    _fovy = 45.0;

    gluPerspective(_fovy, _aspect, _z_near, _z_far);

    // setting the model view matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    _eye[0] = _eye[1] = 0.0, _eye[2] = 6.0;
    _center[0] = _center[1] = _center[2] = 0.0;
    _up[0] = _up[2] = 0.0, _up[1] = 1.0;

    gluLookAt(_eye[0], _eye[1], _eye[2], _center[0], _center[1], _center[2], _up[0], _up[1], _up[2]);

    // enabling depth test
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POLYGON_SMOOTH) ;
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST ) ;
    glHint(GL_PERSPECTIVE_CORRECTION_HINT , GL_NICEST);


    // setting the color of background
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    // initial values of transformation parameters
    _angle_x = _angle_y = _angle_z = 0.0;
    _trans_x = _trans_y = _trans_z = 0.0;
    _zoom = 1.0;

    try
    {
        // initializing the OpenGL Extension Wrangler library
        GLenum error = glewInit();

        if (error != GLEW_OK)
        {
            throw Exception("Could not initialize the OpenGL Extension Wrangler Library!");
        }

        if (!glewIsSupported("GL_VERSION_2_0"))
        {
            throw Exception("Your graphics card is not compatible with OpenGL 2.0+! "
                            "Try to update your driver or buy a new graphics adapter!");
        }

        // create and store your geometry in display lists or vertex buffer objects
        _pc = 0;
        _image_of_pc = 0;
        dl = 0;
        _ips = 0;
        _cc = 0;
        _icc = 0;
        _before_interpolation = 0;
        _after_interpolation = 0;
        // create and store your geometry in display lists or vertex buffer objects


        _model.LoadFromOFF("Models/sphere.off", 1);
        _model.UpdateVertexBufferObjects();

        _n = 5;

        ColumnMatrix<GLdouble> chromosome(2*_n+1);
        GLdouble step = TWO_PI / (2 * _n + 1);

        _u.ResizeColumns(2*_n+1);
        for(GLuint i=0; i<2*_n+1; i++)
        {
            _u[i] = min(i*step, TWO_PI);
            chromosome[i] = _u[i];
        }
        chromosome[6]+=0.4;
        datatoint.ResizeRows(2*_n+1);
        for(GLuint i=0;i<2*_n+1;i++)
        {
            DCoordinate3 &p = datatoint[i];
            p[0] = cos(_u[i]);
            p[1] = sin(_u[i]);
            //p[2] = -1.0 + 2.0 * (GLdouble)rand() / (GLdouble)RAND_MAX;
        }

        // ////////////////////////////////////////////////////////////////////////////////////

//        _cc = new CyclicCurve3(_n);

//        if (!_cc)
//        {
//            throw Exception("Could not generate the cyclic curve!");
//        }
//        for(GLuint i=0;i<2*_n+1;i++)
//        {
//            DCoordinate3 &cp = (*_cc)[i];
//            cp[0] = cos(_u[i]);
//            cp[1] = sin(_u[i]);
//        }
//        if (!_cc->UpdateVertexBufferObjectsOfData())
//        {
//            throw Exception("Could not update the VBO of the cyclic curve's control polygon!");
//        }
//        _mod = 2;
//        _div = 200;
//        _icc3 = _cc->GenerateImage(_mod, _div);
//        cout << _cc->Curvature(100) << endl;
//        cout << _cc->Length(20) << endl;

//        if (!_icc3)
//        {
//            throw Exception("Could not generate the image of the cyclic curve!");
//        }

//        if (!_icc3->UpdateVertexBufferObjects())
//        {
//            throw Exception("Could not update the VBO of the cyclic curve's image!");
//        }
        // ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //        ci = new CurveIndividual(CYCLIC,2*_n+1);
        //        ci->SetChromosome(chromosome);
        //        ci->CalculateFitness(datatoint);
        //        ci->GenerateImage();
        //        cout << "Fitness: " << ci->Fitness() << endl;

        //cp = new CurvePopulation(100,CYCLIC,11,100,100);
        cp = new CurvePopulation(100,CYCLIC,LENGTH,datatoint,100,100);
        //cp->SetDataToInterpolate(datatoint);
        //cp->CalculateFitnesses();
        _icc2 = cp->ImageOfBestIndividual();
        _icc2->UpdateVertexBufferObjects();

        _icc = nullptr;
        _icc = cp->ImageOfBestIndividual();
        _icc->UpdateVertexBufferObjects();

//        _curveTimer->start();

    }
    catch (Exception &e)
    {
        cout << e << endl;
    }
}

//-----------------------
// the rendering function
//-----------------------
void GLWidget::paintGL()
{
    // clears the color and depth buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // stores/duplicates the original model view matrix
    glPushMatrix();

    // applying transformations
    glRotatef(_angle_x, 1.0, 0.0, 0.0);
    glRotatef(_angle_y, 0.0, 1.0, 0.0);
    glRotatef(_angle_z, 0.0, 0.0, 1.0);
    glTranslated(_trans_x, _trans_y, _trans_z);
    glScaled(_zoom, _zoom, _zoom);

    // render your geometry (this is oldest OpenGL rendering technique, later we will use some advanced methods)

    glColor3f(5.0, 5.0, 1.0);
    //    ci->_cc->UpdateVertexBufferObjectsOfData();
    //    if(!ci->_cc->RenderData(GL_POINTS))
    //    {
    //        cout << "False" << endl;
    //    }


    //if (_icc)
    //{
    glColor3f(1.0,1.0,1.0);                         //a feher a kezdeti
    //ci->_icc->RenderDerivatives(0, GL_LINE_STRIP);
    //_icc2->RenderDerivatives(0, GL_LINE_STRIP);
    //sleep(5000);
//    _cc->RenderData(GL_POINTS);
//    _icc3->RenderDerivatives(0,GL_LINE_LOOP);
//    _icc3->RenderDerivatives(1,GL_LINES);
//    _icc3->RenderDerivatives(2,GL_LINES);

    if (_icc)
    {
        //_cc->RenderData(GL_POINTS);
        glColor3f(0.4,0.5,1.0);
        _icc->RenderDerivatives(0,GL_LINE_LOOP);
        glColor3f(0.4,0.7,0.1);
        //_icc->RenderDerivatives(1,GL_LINES);
    }

    //glColor3f(0.7,0.1,0.3);
    //ci->_icc->RenderDerivatives(1,GL_LINES);
    //ci->_icc->RenderDerivatives(1,GL_POINTS);
    //glColor3f(0.2,0.5,0.9);
    //ci->_icc->RenderDerivatives(2,GL_LINES);
    //ci->_icc->RenderDerivatives(2,GL_POINTS);
    //}
    glPointSize(10.0f);
    for(GLuint i=0;i<=2*_n;i++)
    {
        DCoordinate3 &p = datatoint[i];
        glBegin(GL_POINTS); //starts drawing of points
        glVertex3d(p[0],p[1],p[2]);
        glEnd();//end drawing of points
    }

    // pops the current matrix stack, replacing the current matrix with the one below it on the stack,
    // i.e., the original model view matrix is restored
    glPopMatrix();
}

//----------------------------------------------------------------------------
// when the main window is resized one needs to redefine the projection matrix
//----------------------------------------------------------------------------
void GLWidget::resizeGL(int w, int h)
{
    // setting the new size of the rendering context
    glViewport(0, 0, w, h);

    // redefining the projection matrix
    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();

    _aspect = (float)w / (float)h;

    gluPerspective(_fovy, _aspect, _z_near, _z_far);

    // switching back to the model view matrix
    glMatrixMode(GL_MODELVIEW);

    updateGL();
}

//-----------------------------------
// implementation of the public slots
//-----------------------------------

void GLWidget::set_angle_x(int value)
{
    if (_angle_x != value)
    {
        _angle_x = value;
        updateGL();
    }
}

void GLWidget::set_angle_y(int value)
{
    if (_angle_y != value)
    {
        _angle_y = value;
        updateGL();
    }
}

void GLWidget::set_angle_z(int value)
{
    if (_angle_z != value)
    {
        _angle_z = value;
        updateGL();
    }
}

void GLWidget::set_zoom_factor(double value)
{
    if (_zoom != value)
    {
        _zoom = value;
        updateGL();
    }
}

void GLWidget::set_trans_x(double value)
{
    if (_trans_x != value)
    {
        _trans_x = value;
        updateGL();
    }
}

void GLWidget::set_trans_y(double value)
{
    if (_trans_y != value)
    {
        _trans_y = value;
        updateGL();
    }
}

void GLWidget::set_trans_z(double value)
{
    if (_trans_z != value)
    {
        _trans_z = value;
        updateGL();
    }
}

void GLWidget::_evolveCurves()
{
    cout << iterations << endl;
    iterations++;
    cp->Evolve(0.6,0.35,0.4,10);
    if (_icc)
    {
        delete _icc, _icc = nullptr;
    }

    _icc = cp->ImageOfBestIndividual();
    _icc->UpdateVertexBufferObjects();
    updateGL();
}

void GLWidget::set_evolve_state(bool value)
{
    if (_evolve != value)
    {
        _evolve = value;

        if (_evolve)
        {
            _curveTimer->start();
        }
        else
        {
            _curveTimer->stop();
        }

        updateGL();
    }
}

}
