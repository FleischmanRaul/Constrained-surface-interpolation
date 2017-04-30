#pragma once

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <random>
#include <algorithm>    // std::swap
#include <GL/glew.h>
#include "Core/Exceptions.h"
#include "../Core/Constants.h"
#include "../Core/Matrices.h"
#include "../Core/LinearCombination3.h"
#include "../Core/GenericCurves3.h"
#include "Cyclic/CyclicCurves3.h"
namespace cagd
{
    typedef GLdouble Gene;
    typedef ColumnMatrix<Gene> Chromosome;
    enum    CurveType{BSPLINE, CYCLIC /*, etc.*/};
    enum    CurveEnergyType{LENGTH, CURVATURE};

    class CurveIndividual
    {
    public:
    protected:
        CurveType       _type;
        GLdouble        _fitness;
        CurveEnergyType _etype;

    public:
        Chromosome            _chromosome;
        LinearCombination3*   _cc;
        // default constructor
        CurveIndividual(GLuint geneNumber = 0);
        // special constructor
        CurveIndividual(CurveType type, CurveEnergyType etype,GLuint geneNumber = 0);

        GLboolean SetChromosome(ColumnMatrix<GLdouble> chromosome);

        Gene  operator [](GLuint i) const;
        Gene& operator [](GLuint i);

        GLvoid          GetDefinitionDomain(GLdouble& u_min, GLdouble& u_max) const;
        GLboolean       GenerateChromosome();
        GenericCurve3*  GenerateImage();
        GLboolean       CalculateFitness(const ColumnMatrix<DCoordinate3> &dataToInterpolate, GLuint n);
        GLdouble        Fitness() const;
    };

    class CurvePopulation
    {
    protected:
        GLuint                     _geneNumber;
        GLuint                     _indexOfBestIndividual;
        GLuint                     _maxMaturityLevel;
        GLuint                     _divPointCount;
        GLdouble                   _threshold;
        ColumnMatrix<DCoordinate3> _dataToInterpolate;

    public:
        RowMatrix<CurveIndividual> _individual;
        CurvePopulation(GLuint individualCount, CurveType type, CurveEnergyType etype, const ColumnMatrix<DCoordinate3> &dataToInterpolate, GLuint maxMaturityLevel, GLuint _divPointCount, GLdouble threshold = EPS);


        GLvoid Mutation(GLdouble probability, GLdouble mutationRadiusPercentage);
        GLvoid Recombination(GLdouble probability);
        GLvoid TournamentSelection(GLuint poolSize);

        GLvoid Evolve(GLdouble mutationProbability, GLdouble mutationRadiusPercentage,
                      GLdouble recombinationProbability,
                      GLuint poolSize);

        GLboolean      CalculateFitnesses();
        GLboolean      SetDataToInterpolate(const ColumnMatrix<DCoordinate3>& dataToInterpolate);
        GLuint         FindBestIndividual(const RowMatrix<GLuint> &pool) const;
        GenericCurve3* ImageOfIndividual(GLuint i) const;
        GenericCurve3* ImageOfBestIndividual();
    };
}
