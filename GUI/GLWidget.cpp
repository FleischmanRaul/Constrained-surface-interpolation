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
    _surfaceTimer = new QTimer(this);
    _surfaceTimer->setInterval(0);
    connect(_curveTimer,   SIGNAL(timeout()), this, SLOT(_evolveCurves()));
    connect(_surfaceTimer, SIGNAL(timeout()), this, SLOT(_evolveSurfaces()));
}

//--------------------------------------------------------------------------------------
// this virtual function is called once before the first call to paintGL() or resizeGL()
//--------------------------------------------------------------------------------------
DirectionalLight * initializeDirectionalLight()
{
    DirectionalLight * _dl = 0;

    HCoordinate3 direction(0.0,0.0,1.0,0.0);
    Color4      ambient   (0.4,0.4,0.4,1.0);
    Color4      diffuse   (0.8,0.8,0.8,1.0);
    Color4      specular   (1.0,1.0,1.0,1.0);

    _dl = new DirectionalLight(GL_LIGHT0, direction, ambient, diffuse, specular);
    return _dl;
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
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH) ;
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
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
        _dl = 0;
        _ips = 0;
        _cc = 0;
        _icc = 0;
        _icc = nullptr;
        _before_interpolation = nullptr;
        _after_interpolation = 0;
        // create and store your geometry in display lists or vertex buffer objects


//        _model.LoadFromOFF("Models/sphere.off", 1);
//        _model.UpdateVertexBufferObjects();

        _n = 4;

        ColumnMatrix<GLdouble> chromosome  (2*_n+1);
        RowMatrix   <GLdouble> u_chromosome(2*_n+1);
        ColumnMatrix<GLdouble> v_chromosome(2*_n+1);
        GLdouble tart = PI;
        GLdouble step = tart / (2*_n);

        _u.ResizeColumns(2*_n+1);
        for(GLuint i=0; i<2*_n+1; i++)
        {
            _u[i] = min(i*step, tart);
            chromosome[i] = _u[i];
        }

        for(GLuint i=0; i<2*_n+1; i++)
        {
            _u[i] = min(i*step, tart);
            u_chromosome[i] = _u[i];
            v_chromosome[i] = _u[i];
        }

//        u_chromosome[0] = 0;
//        u_chromosome[1] = 0.255996;
//        u_chromosome[2] = 0.327781; 0.419147 0.681159 1.21572 1.29602 1.31108 2.02642 2.38215 3.14159

        //chromosome[0]+=0.2;
        datatoint.ResizeRows(2*_n+1);
        for(GLuint i=0;i<2*_n+1;i++)
        {
            DCoordinate3 &p = datatoint[i];
            p[0] = cos(_u[i]);
            p[1] = sin(_u[i]);
            //p[2] = -1.0 + 2.0 * (GLdouble)rand() / (GLdouble)RAND_MAX;
        }

        // ////////////////////////////////////////////////////////////////////////////////////

        controll_nett.ResizeColumns(2*_n+1);
        controll_nett.ResizeRows(2*_n+1);
        _Bsurface = new TrigonometricBernsteinSurface3(tart,_n,tart,_n);
        for(GLuint i = 0; i < 2*_n+1; i++)
        {
            for(GLuint j = 0; j<2*_n+1; j++)
            {

                //DCoordinate3 &p = (*_Bsurface)(i,j);//controll_nett(i,j);
                DCoordinate3 &p = controll_nett(i,j);
                p[0] = -_n * 0 +  (GLdouble)i/2;
                p[1] = -_n * 0 +  (GLdouble)j/2;//-1.0 + 2.0 * (GLdouble)rand() / (GLdouble)RAND_MAX;
                //p[2] = (GLdouble)rand() / (GLdouble)RAND_MAX;
            }
        }

        DCoordinate3 &p = controll_nett(_n,_n);
        p[2] = 1.34;

//        if (!_Bsurface->UpdateDataForInterpolation(u_chromosome,v_chromosome,controll_nett))
//        {
//            throw Exception("Could not update the VBO!");
//        }
//        _before_interpolation = _Bsurface->GenerateImage(100,100, TensorProductSurface3::LOGARITHMIC_UMBILIC_DEVIATION_ENERGY_FRAGMENT);
//        _Bsurface->UpdateVertexBufferObjectsOfData();
//        if(!_before_interpolation->UpdateVertexBufferObjects())
//        {
//            cout << "Error updating the vertex buffer objects" << endl;
//        }
//        cout << "Umbilic deviation: " << _Bsurface->Fitness(30,30,UMBILIC_DEVIATION) << endl;
//        cout << "Total curvature: "     << _Bsurface->Fitness(30,30,TOTAL_CURVATURE) << endl;
//        cout << "Willmore: " << _Bsurface->Fitness(30,30,WILLMORE) << endl;

//        si = new SurfaceIndividual(BERNSTEIN,TOTAL_CURVATURE,9,9);
//        if(!si->GenerateChromosomes())
//        {
//            throw Exception("error in generating the chromosomes!");
//        }
//        if(!si->CalculateFitness(controll_nett,30,30))
//        {
//            throw Exception("error in calculating fitness!");
//        }
//        cout << si->Fitness() << endl;
//        _before_interpolation = si->GenerateImage(TensorProductSurface3::LOGARITHMIC_TOTAL_CURVATURE_ENERGY_FRAGMENT);
//        cout << "hello" << endl;
        sp = new SurfacePopulation(100, BERNSTEIN, WILLMORE,controll_nett,20,TensorProductSurface3::LOGARITHMIC_WILLMORE_ENERGY_FRAGMENT);
        _before_interpolation = sp->ImageOfBestIndividual();
        cout << sp->FitnessOfBestIndividual() << endl;
        if(!_before_interpolation->UpdateVertexBufferObjects())
        {
            cout << "Error updating the vertex buffer objects" << endl;
        }

        if (!_shader.InstallShaders("Shaders/two_sided_lighting_color.vert", "Shaders/two_sided_lighting_color.frag"))
        {
            // throw
        }


    }
    catch (Exception &e)
    {
        cout << e << endl;
    }

    DirectionalLight * _dl2 = 0;

    HCoordinate3 direction(0.0,1.0,0.0,0.0);
    Color4      ambient   (0.4,0.4,0.4,1.0);
    Color4      diffuse   (0.8,0.8,0.8,1.0);
    Color4      specular   (1.0,1.0,1.0,1.0);

//    _dl2 = new DirectionalLight(GL_LIGHT2, direction, ambient, diffuse, specular);
    _dl = initializeDirectionalLight();
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHT0);
//    glEnable(GL_LIGHT2);
    _dl->Enable();
//    _dl2->Enable();
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

    if (_icc)
    {
        //_cc->RenderData(GL_POINTS);
        glColor3f(0.0,0.5,1.0);
        _icc->RenderDerivatives(0,GL_LINE_STRIP);
        //glColor3f(0.4,0.7,0.1);
        //_icc->RenderDerivatives(1,GL_LINES);
    }
    glPointSize(5.0f);
    glLineWidth(1);
    //_Bsurface->RenderData();

//    glDisable(GL_LIGHTING);
//    glColor3f(1.0, 1.0, 1.0);
//    _before_interpolation->RenderNormals();

//    glEnable(GL_LIGHTING);


    MatFBRuby.Apply();
    _shader.Enable();
    if(!_before_interpolation->Render(GL_TRIANGLES))
    {
        cout << "Error in rendering" << endl;
    }
    _shader.Disable();
    MatFBGold.Apply();
    for(GLuint i = 0; i<2*_n+1;++i)
    {
        for(GLuint j = 0; j<2*_n+1;++j)
        {
            glColor3f(0.2,0.5,0.9);
            DCoordinate3 &p = controll_nett(i,j);
            glBegin(GL_POINTS); //starts drawing of points
            glVertex3d(p[0],p[1],p[2]);
            glEnd();//end drawing of points
        }
    }
//    if(!_Bsurface->RenderData(GL_POINTS))
//    {
//        cout << "Error in rendering derivatives" << endl;
//    }
//    glColor3f(0.7,0.1,0.3);
//              --Curve interpolation
//    for(GLuint i=0;i<=2*_n;i++)
//    {
//        DCoordinate3 &p = datatoint[i];
//        glBegin(GL_POINTS); //starts drawing of points
//        glVertex3d(p[0],p[1],p[2]);
//        glEnd();//end drawing of points
//    }

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
    //cout << curve_iterations++ << endl;
    cp->Evolve(0.6,mutation_radius,0.4,10);
    if (_icc)
    {
        delete _icc, _icc = nullptr;
    }

    _icc = cp->ImageOfBestIndividual();
    _icc->UpdateVertexBufferObjects();
    new_fitness = cp->FitnessOfBestIndividual();
    //cout << -new_fitness << endl;
    if(new_fitness != previous_fitness)
    {
        curve_generations++;
        curve_subgenerations = 0;
    }
    else
    {
        if(curve_subgenerations > _maxMaturityLevel)
        {
            mutation_radius/=2;
            curve_subgenerations = 0;
            emit valueChanged(mutation_radius);
        }
        else
        {
            curve_subgenerations++;
        }
    }
    previous_fitness = new_fitness;
    cout << previous_fitness << endl;
    if(mutation_radius < _threshold)
    {
        //curveEvolveCheckBox -> setChecked(false);
        set_evolve_state(false);
    }
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
            cout << "threshold: " << _threshold << endl;
            _curveTimer->stop();
        }

        updateGL();
    }
}

void GLWidget::_evolveSurfaces()
{
    sp->Evolve(1.0,0.5,0.5,1);
    cout << "it: " <<curve_iterations++ << endl;
    _before_interpolation = sp->ImageOfBestIndividual();
    cout << sp->FitnessOfBestIndividual() << endl;
    if(!_before_interpolation->UpdateVertexBufferObjects())
    {
        cout << "Error updating the vertex buffer objects" << endl;
    }
}

void GLWidget::set_evolve_state2(bool value)
{
    if (_evolve != value)
    {
        _evolve = value;
        if (_evolve)
        {
            curve_iterations = 0;
            _surfaceTimer->start();
        }
        else
        {
            _surfaceTimer->stop();
        }
        updateGL();
    }
}

void GLWidget::set_curve_energy_type(int value)
{
    _curveEnergyType = (CurveEnergyType)value;
}
void GLWidget::set_curve_type(int value)
{
    _curveType = (CurveType)value;
}
void GLWidget::set_population_count(int value)
{
    _populationCount = (GLuint)value;
}
void GLWidget::set_max_maturity_level(int value)
{
    _maxMaturityLevel = (GLuint)value;
}
void GLWidget::set_mating_pool_size(int value)
{
    _matingPoolCount = (GLuint)value;
}
void GLWidget::set_threshold(double value)
{
    _threshold = value/1000;
}
void GLWidget::create_curve_population()
{
    curve_iterations = 0;
    mutation_radius  = 1.0;
    previous_fitness = 0.0;
    RowMatrix<GLdouble> energyprop(2);
    energyprop[0] = 3;
    energyprop[1] = 1;

    cp = new CurvePopulation(_populationCount,_curveType,_curveEnergyType,datatoint,_matingPoolCount, energyprop);
    _icc = cp->ImageOfBestIndividual();
    _icc->UpdateVertexBufferObjects();
}
GLWidget::~GLWidget()
{
    if(_cc)
        delete _cc, _cc=0;
    if(_icc)
        delete _icc, _icc=0;
    if(_pc)
        delete _pc, _pc=0;
    if(_image_of_pc)
        delete _image_of_pc, _image_of_pc=0;
    if(_dl)
        delete _dl, _dl = 0;
    if(_ips)
        delete _ips, _ips = 0;
    if(_before_interpolation)
        delete _before_interpolation, _before_interpolation = 0;
    if(_after_interpolation)
        delete _after_interpolation, _after_interpolation = 0;
}
}
