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
#include "../Core/LinearCombination3.h"
#include "../Core/GenericCurves3.h"
#include "Cyclic/CyclicCurves3.h"
#include "Trigonometric/TrigonometricCurves3.h"
namespace cagd
{
    typedef GLdouble Gene;
    typedef ColumnMatrix<Gene> Chromosome;
    enum    CurveType{BSPLINE, CYCLIC, TRIGONOMETRIC /*, etc.*/};
    enum    CurveEnergyType{LENGTH, CURVATURE, MIXED};

    class CurveIndividual
    {
    public:
    protected:
        CurveType           _type;
        GLdouble            _fitness;
        CurveEnergyType     _etype;
        Chromosome          _chromosome;
        RowMatrix<GLdouble> _eProportion;

    public:
        LinearCombination3*   _cc;
        // default constructor
        CurveIndividual(GLuint geneNumber = 0);
        // copy constructor
        //CurveIndividual(const CurveIndividual& ci);
        // special constructor
        CurveIndividual(CurveType type, CurveEnergyType etype,GLuint geneNumber = 0);

        GLboolean SetChromosome(ColumnMatrix<GLdouble> chromosome);

        Gene  operator [](GLuint i) const;
        Gene& operator [](GLuint i);

        GLvoid          GetDefinitionDomain(GLdouble& u_min, GLdouble& u_max) const;
        GLvoid          SetEnergyProportions(const RowMatrix<GLdouble>& eProportion);
        GLboolean       GenerateChromosome();
        GenericCurve3*  GenerateImage();
        GLboolean       CalculateFitness(const ColumnMatrix<DCoordinate3> &dataToInterpolate, GLuint n);
        GLdouble        Fitness() const;

        virtual ~CurveIndividual();
    };

    class CurvePopulation
    {
    protected:
        GLuint                     _geneNumber;
        GLuint                     _indexOfBestIndividual;
        GLuint                     _maxMaturityLevel;
        ColumnMatrix<DCoordinate3> _dataToInterpolate;
        GLuint                     _divPointCount;
        GLdouble                   _threshold;
        RowMatrix<CurveIndividual> _individual;

    public:
        CurvePopulation(GLuint individualCount, CurveType type, CurveEnergyType etype, const ColumnMatrix<DCoordinate3> &dataToInterpolate,
                        GLuint maxMaturityLevel, GLuint _divPointCount,RowMatrix<GLdouble>& eProportion, GLdouble threshold = EPS);


        GLvoid Mutation(GLdouble probability, GLdouble mutationRadiusPercentage);
        GLvoid Recombination(GLdouble probability);
        GLvoid TournamentSelection(GLuint poolSize);

        GLvoid Evolve(GLdouble mutationProbability, GLdouble mutationRadiusPercentage,
                       GLdouble recombinationProbability,
                       GLuint poolSize);


        GLboolean      CalculateFitnesses();
        GLboolean      SetDataToInterpolate(const ColumnMatrix<DCoordinate3>& dataToInterpolate);
        GLuint         FindBestIndividual(const RowMatrix<GLuint> &pool) const;
        GLuint         GetMaxMaturityLevel();
        GLdouble       GetThreshold();
        GLdouble       FitnessOfBestIndividual();

        GenericCurve3* ImageOfIndividual(GLuint i) const;
        GenericCurve3* ImageOfBestIndividual();
    };
}
