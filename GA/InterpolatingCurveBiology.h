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
    enum    CurveType{CYCLIC, TRIGONOMETRIC /*, etc.*/};
    enum    CurveEnergyType{LENGTH, CURVATURE, MIXED /*, etc.*/};

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
        CurveIndividual(const CurveIndividual& ci);
        // special constructor
        CurveIndividual(CurveType type, CurveEnergyType etype,GLuint geneNumber = 0);

        CurveIndividual& operator =(const CurveIndividual& rhs);

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
        ColumnMatrix<DCoordinate3> _dataToInterpolate;
        GLuint                     _divPointCount;
        RowMatrix<CurveIndividual> _individual;

    public:
        CurvePopulation(GLuint individualCount, CurveType type, CurveEnergyType etype, const ColumnMatrix<DCoordinate3> &dataToInterpolate,
                        GLuint divPointCount,RowMatrix<GLdouble>& eProportion);


        GLvoid Mutation(GLdouble probability, GLdouble mutationRadiusPercentage);
        GLvoid Recombination(GLdouble probability);
        GLvoid TournamentSelection(GLuint poolSize);

        GLvoid Evolve(GLdouble mutationProbability, GLdouble mutationRadiusPercentage,
                       GLdouble recombinationProbability,
                       GLuint poolSize);


        GLboolean      CalculateFitnesses();
        GLboolean      SetDataToInterpolate(const ColumnMatrix<DCoordinate3>& dataToInterpolate);
        GLuint         FindBestIndividual(const RowMatrix<GLuint> &pool) const;
        GLdouble       FitnessOfBestIndividual();

        GenericCurve3* ImageOfIndividual(GLuint i) const;
        GenericCurve3* ImageOfBestIndividual();
    };
}
