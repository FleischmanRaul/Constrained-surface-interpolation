#pragma once

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <random>
#include <climits>
#include <algorithm>    // std::swap
#include <GL/glew.h>
#include "Core/Exceptions.h"
#include "../Core/Constants.h"
#include "../Core/Matrices.h"
#include "../Core/TensorProductSurfaces3.h"
#include "../Trigonometric/TrigonometricBernsteinSurfaces.h"
#include "../Cyclic/CyclicSurfaces3.h"
namespace cagd
{
    typedef GLdouble Gene;
    typedef RowMatrix   <Gene>UChromosome;
    typedef ColumnMatrix<Gene>VChromosome;
    enum    SurfaceType{BERNSTEIN, CYCLICSURF /*, etc.*/};


    class SurfaceIndividual
    {
    public:
    protected:
        SurfaceType          _type;
        GLdouble             _fitness;
        SurfaceEnergyType    _etype;
        TensorProductSurface3*   _cc;
    public:
        UChromosome          _uchromosome;
        VChromosome          _vchromosome;
        // default constructor
        SurfaceIndividual(GLuint u_geneNumber = 0, GLuint v_geneNumber = 0);
        // copy constructor
        SurfaceIndividual(const SurfaceIndividual& ci);
        // special constructor
        SurfaceIndividual(SurfaceType type, SurfaceEnergyType etype, GLuint u_geneNumber = 0, GLuint v_geneNumber = 0);
        //assignment operator
        SurfaceIndividual& operator =(const SurfaceIndividual& rhs);

        GLboolean SetUChromosome(RowMatrix   <GLdouble> uchromosome);
        GLboolean SetVChromosome(ColumnMatrix<GLdouble> vchromosome);
        GLboolean SetChromosomes(RowMatrix   <GLdouble> uchromosome, ColumnMatrix<GLdouble> vchromosome);

        Gene  operator [](GLuint i) const;      //for the U chromosome
        Gene& operator [](GLuint i);

        Gene  operator ()(GLuint i) const;      //for the V chromosome
        Gene& operator ()(GLuint i);

        GLvoid              GetDefinitionDomain(GLdouble& u_min, GLdouble& u_max, GLdouble& v_min, GLdouble& v_max) const;
        GLboolean           GenerateChromosomes();
        TriangulatedMesh3*  GenerateImage(TensorProductSurface3::ImageColorScheme colorscheme);
        GLboolean           CalculateFitness(Matrix<DCoordinate3> &dataToInterpolate, GLuint n, GLuint m);
        GLdouble            Fitness() const;

        virtual ~SurfaceIndividual();
    };

    class SurfacePopulation
    {
    protected:
        GLuint                      _u_geneNumber;
        GLuint                      _v_geneNumber;
        GLint                       _indexOfBestIndividual;
        Matrix<DCoordinate3>        _dataToInterpolate;
        GLuint                      _divPointCount;
        RowMatrix<SurfaceIndividual>_individual;
        TensorProductSurface3::ImageColorScheme _colorscheme;

    public:
        SurfacePopulation(GLuint individualCount, SurfaceType type, SurfaceEnergyType etype, Matrix<DCoordinate3> &dataToInterpolate,
                        GLuint divPointCount, TensorProductSurface3::ImageColorScheme colorscheme);


        GLvoid Mutation(GLdouble probability, GLdouble mutationRadiusPercentage);
        GLvoid Recombination(GLdouble probability);
        GLvoid TournamentSelection(GLuint poolSize);

        GLvoid Evolve(GLdouble mutationProbability, GLdouble mutationRadiusPercentage,
                       GLdouble recombinationProbability,
                       GLuint poolSize);

        GLboolean      SetDataToInterpolate(const Matrix<DCoordinate3>& dataToInterpolate);
        GLuint         FindBestIndividual(const RowMatrix<GLuint> &pool) const;
        GLdouble       FitnessOfBestIndividual();

        TriangulatedMesh3* ImageOfIndividual(GLuint i) const;
        TriangulatedMesh3* ImageOfBestIndividual();
    };

}
