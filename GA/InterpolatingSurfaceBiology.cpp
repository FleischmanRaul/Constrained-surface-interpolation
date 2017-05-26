#include "GA/InterpolatingSurfaceBiology.h"
#include <vector>       // std::vector
#include <limits>

using namespace cagd;
using namespace std;

std::random_device rd2;
std::mt19937 generator2(rd2());
std::uniform_real_distribution<GLdouble> prob_distribution2(0.0,1.0);

//-------------------------- Implementation of SurfaceIndividual --------------------------
SurfaceIndividual::SurfaceIndividual(GLuint u_geneNumber, GLuint v_geneNumber)
{
    _uchromosome.ResizeColumns(u_geneNumber);
    _vchromosome.ResizeRows(v_geneNumber);
    _cc = nullptr;
}

SurfaceIndividual::SurfaceIndividual(SurfaceType type, SurfaceEnergyType etype, GLuint u_geneNumber, GLuint v_geneNumber):
    _type(type),_etype(etype)
{
    _uchromosome.ResizeColumns(u_geneNumber);
    _vchromosome.ResizeRows(v_geneNumber);
    switch(_type)
    {
    case BERNSTEIN:
        _cc = new TrigonometricBernsteinSurface3(PI, u_geneNumber/2, PI, v_geneNumber/2);
        break;
    case CYCLICSURF:
        _cc = new CyclicSurface3(u_geneNumber/2, v_geneNumber/2);
        break;
    default: std::cout << "Error in constructing " << _type << " surface!" << endl;
                 break;
    }
    _fitness = std::numeric_limits<GLdouble>::min();
}

SurfaceIndividual::SurfaceIndividual(const SurfaceIndividual& si):
_type(si._type),_fitness(si._fitness),_etype(si._etype),_uchromosome(si._uchromosome),_vchromosome(si._vchromosome)
{
    if(si._cc != nullptr)
    {
        _cc = si._cc->Clone();
    }
    else
    {
        _cc = nullptr;
    }
}

SurfaceIndividual& SurfaceIndividual::operator=(const SurfaceIndividual& rhs)
{
    if (this == &rhs) return *this; // handle self assignment

    _type        = rhs._type;
    _fitness     = rhs._fitness;
    _etype       = rhs._etype;
    _uchromosome  = rhs._uchromosome;
    _vchromosome  = rhs._vchromosome;

    if(_cc)
    {
        delete _cc;
    }
    _cc = nullptr;
    if(rhs._cc)
    {
        _cc = rhs._cc->Clone();
    }
    return *this;
}

GLboolean SurfaceIndividual::SetUChromosome(RowMatrix<GLdouble> uchromosome)
{
    if(_uchromosome.GetRowCount() != uchromosome.GetRowCount())
    {
        return GL_FALSE;
    }
    _uchromosome = uchromosome;
    return GL_TRUE;
}

GLboolean SurfaceIndividual::SetVChromosome(ColumnMatrix<GLdouble> vchromosome)
{
    if(_vchromosome.GetRowCount() != vchromosome.GetRowCount())
    {
        return GL_FALSE;
    }
    _vchromosome = vchromosome;
    return GL_TRUE;
}

GLboolean SurfaceIndividual::SetChromosomes(RowMatrix<GLdouble> uchromosome, ColumnMatrix<GLdouble> vchromosome)
{
    return SetUChromosome(uchromosome) && SetVChromosome(vchromosome);
}

Gene SurfaceIndividual::operator [](GLuint i) const
{
    return _uchromosome[i];
}


Gene& SurfaceIndividual::operator [](GLuint i)
{
    return _uchromosome[i];
}


Gene SurfaceIndividual::operator ()(GLuint i) const
{
    return _vchromosome[i];
}


Gene& SurfaceIndividual::operator ()(GLuint i)
{
    return _vchromosome[i];
}

GLvoid SurfaceIndividual::GetDefinitionDomain(GLdouble& u_min, GLdouble& u_max, GLdouble& v_min, GLdouble& v_max) const
{
    _cc->GetDefinitionDomain(u_min,u_max,v_min,v_max);
}

GLdouble SurfaceIndividual::Fitness() const
{
    return _fitness;
}

GLboolean SurfaceIndividual::GenerateChromosomes()
{
    GLuint u_size = _uchromosome.GetColumnCount();
    GLuint v_size = _vchromosome.GetRowCount();
    GLdouble u_min, u_max, v_min, v_max;
    GetDefinitionDomain(u_min,u_max, v_min, v_max);
    std::uniform_real_distribution<GLdouble> distribution(0.0, 0.6);

    GLdouble u_step = (u_max - u_min) / ((GLdouble)u_size - 1);
    GLdouble v_step = (v_max - v_min) / ((GLdouble)v_size - 1);

    for(GLuint i = 0; i < u_size; i++)
    {
        _uchromosome[i] = (GLdouble)i * u_step;
    }

    for(GLuint i = 0; i < u_size; i++)
    {
        _vchromosome[i] = (GLdouble)i * v_step;
    }

    for(GLuint i = 0; i < u_size; i++)          //perturbating the u knot vector
    {
        if(i == 0)
        {
            if(_type != BERNSTEIN)              //if Bernstein the first knot must stay the same
            {
                _uchromosome[i] += (_uchromosome[i+1] - _uchromosome[i]) * distribution(generator2);
            }
        }
        else if(i == (u_size - 1))
        {
            if(_type != BERNSTEIN)              //if Bernstein the last knot must stay the same
            {
                _uchromosome[i] -= (_uchromosome[i] - _uchromosome[i-1]) * distribution(generator2);
            }
        }
        else
        {
            if(distribution(generator2) < 0.3)  //move to the left
            {
                _uchromosome[i] -= (_uchromosome[i] - _uchromosome[i-1]) * distribution(generator2);
            }
            else                                //move to the right
            {
                _uchromosome[i] += (_uchromosome[i+1] - _uchromosome[i]) * distribution(generator2);
            }
        }
    }

    for(GLuint i = 0; i < v_size; i++)           //perturbating the v knot vector
    {
        if(i == 0)
        {
            if(_type != BERNSTEIN)               //if Bernstein the first knot must stay the same
            {
                _vchromosome[i] += (_vchromosome[i+1] - _vchromosome[i]) * distribution(generator2);
            }
        }
        else if(i == (u_size - 1))
        {
            if(_type != BERNSTEIN)               //if Bernstein the last knot must stay the same
            {
                _vchromosome[i] -= (_vchromosome[i] - _vchromosome[i-1]) * distribution(generator2);
            }
        }
        else
        {
            if(distribution(generator2) < 0.3)  //move to the left
            {
                _vchromosome[i] -= (_vchromosome[i] - _vchromosome[i-1]) * distribution(generator2);
            }
            else                                //move to the right
            {
                _vchromosome[i] += (_vchromosome[i+1] - _vchromosome[i]) * distribution(generator2);
            }
        }
    }
    return GL_TRUE;
}

TriangulatedMesh3* SurfaceIndividual::GenerateImage(TensorProductSurface3::ImageColorScheme colorscheme)
{
    return _cc->GenerateImage(50,50,colorscheme);
}

GLboolean SurfaceIndividual::CalculateFitness(Matrix<DCoordinate3> &dataToInterpolate, GLuint n, GLuint m)
{

    //cout << _uchromosome.GetColumnCount() << "  " << dataToInterpolate.GetColumnCount() << endl;
    if(n == 0 || m == 0)
    {
        return GL_FALSE;
    }
    if(_uchromosome.GetColumnCount() != dataToInterpolate.GetColumnCount() || _uchromosome.GetColumnCount() == 0 ||
       _vchromosome.GetRowCount() != dataToInterpolate.GetRowCount() || _vchromosome.GetRowCount() == 0)
    {
        return GL_FALSE;
    }

    if(!_cc->UpdateDataForInterpolation(_uchromosome, _vchromosome, dataToInterpolate))
    {
        return GL_FALSE;
    }
    _fitness = _cc->Fitness(n, m, _etype);
    return GL_TRUE;
}

SurfaceIndividual::~SurfaceIndividual()
{

}

//------------------------------------Implementation of SurfacePopulation--------------------------------
//GLuint                      _u_geneNumber;
//GLuint                      _v_geneNumber;
//GLuint                      _indexOfBestIndividual;
//Matrix<DCoordinate3>        _dataToInterpolate;
//GLuint                      _divPointCount;
//RowMatrix<SurfaceIndividual>_individual;

SurfacePopulation::SurfacePopulation(GLuint individualCount, SurfaceType type, SurfaceEnergyType etype, Matrix<DCoordinate3> &dataToInterpolate,
                                     GLuint divPointCount, TensorProductSurface3::ImageColorScheme colorscheme):
    _u_geneNumber(dataToInterpolate.GetRowCount()),
    _v_geneNumber(dataToInterpolate.GetColumnCount()),
    _indexOfBestIndividual(-1),
    _dataToInterpolate(dataToInterpolate),
    _divPointCount(divPointCount),
    _individual(individualCount),
    _colorscheme(colorscheme)
{
    for(GLuint i=0; i<individualCount; ++i)
    {
        _individual[i] = SurfaceIndividual(type,etype,_u_geneNumber,_v_geneNumber);
        _individual[i].GenerateChromosomes();
        if(!_individual[i].CalculateFitness(dataToInterpolate,_divPointCount,_divPointCount))
        {
            cout << "Error calculating fitness" << endl;
        }
    }
}

GLboolean SurfacePopulation::SetDataToInterpolate(const Matrix<DCoordinate3> &dataToInterpolate)
{
    if(dataToInterpolate.GetColumnCount() != _u_geneNumber)
    {
        return GL_FALSE;
    }
    if(dataToInterpolate.GetRowCount() != _v_geneNumber)
    {
        return GL_FALSE;
    }
    _dataToInterpolate = dataToInterpolate;

//    for(GLuint i=0; i<individualCount; ++i)
//    {
//        _individual[i].CalculateFitness(_dataToInterpolate,_divPointCount);
//    }

    return GL_TRUE;
}

GLdouble SurfacePopulation::FitnessOfBestIndividual()
{/*
    cout << _individual[_indexOfBestIndividual]._uchromosome << endl;

    cout << _individual[_indexOfBestIndividual]._vchromosome << endl;*/

    return _individual[_indexOfBestIndividual].Fitness();
}

TriangulatedMesh3* SurfacePopulation::ImageOfIndividual(GLuint i) const
{
    return _individual[i].GenerateImage(_colorscheme);
}

TriangulatedMesh3* SurfacePopulation::ImageOfBestIndividual()
{
    if(_individual.GetColumnCount() < 1)
    {
        cout << "Empty Curve Population!" << endl;
    }
    _indexOfBestIndividual = 0;
    for(GLuint i = 1; i<_individual.GetColumnCount(); ++i)
    {
        if(_individual[i].Fitness() > _individual[_indexOfBestIndividual].Fitness())
        {
            _indexOfBestIndividual = i;
        }
    }
    _individual[_indexOfBestIndividual].CalculateFitness(_dataToInterpolate,_divPointCount,_divPointCount);
    return  _individual[_indexOfBestIndividual].GenerateImage(_colorscheme);
}

GLvoid SurfacePopulation::Recombination(GLdouble probability)
{
    std::uniform_int_distribution <GLuint>   int_distribution(0, _individual.GetColumnCount() - 1);
    GLuint r,w;
    GLdouble u_min, u_max, v_min, v_max, value;
    _individual[0].GetDefinitionDomain(u_min, u_max, v_min, v_max);        //all individuals should have the same definition domain
    std::uniform_real_distribution<GLdouble> u_distribution(u_min,u_max);
    std::uniform_real_distribution<GLdouble> v_distribution(v_min,v_max);
    if(prob_distribution2(generator2) < probability)
    {
        //cout << "In recombination" << endl;
        r = int_distribution(generator2);
        w = int_distribution(generator2);
        if(r==w)                        //they cant have the same value so we change one of them
        {                               //staying between 0 and _individual.GetColumnCount()-1
            if(r==0)
            {
                r++;
            }
            else
            {
                r--;
            }
        }
        SurfaceIndividual rOffspring = _individual[r];
        SurfaceIndividual wOffspring = _individual[w];
        value = u_distribution(generator2);
        GLuint i = 0;

        while(rOffspring[i] < value && wOffspring[i] < value && i < _u_geneNumber-2)
        {
            swap(rOffspring[i],  wOffspring[i]);
            i++;
        }
        value = v_distribution(generator2);
        while(rOffspring(i) < value && wOffspring(i) < value && i < _v_geneNumber-2)
        {
            swap(rOffspring[i],  wOffspring[i]);
            i++;
        }

        rOffspring.CalculateFitness(_dataToInterpolate, _divPointCount, _divPointCount);
        wOffspring.CalculateFitness(_dataToInterpolate, _divPointCount, _divPointCount);

        if(rOffspring.Fitness() > _individual[r].Fitness())
        {
            swap(rOffspring,_individual[r]);
        }
        if(wOffspring.Fitness() > _individual[w].Fitness())
        {
            swap(wOffspring,_individual[w]);
        }
        if(rOffspring.Fitness() > _individual[w].Fitness())
        {
            swap(rOffspring,_individual[w]);
        }
        if(wOffspring.Fitness() > _individual[r].Fitness())
        {
            swap(wOffspring,_individual[r]);
        }
        //_individual[_indexOfBestIndividual].CalculateFitness(_dataToInterpolate, _divPointCount);
    }
}

GLvoid SurfacePopulation::Mutation(GLdouble probability, GLdouble mutationRadiusPercentage)
{
    std::uniform_int_distribution <GLuint>   surface_distribution(0, _individual.GetColumnCount() - 1);
    std::uniform_int_distribution <GLuint>   u_distribution(1, _u_geneNumber - 2);
    std::uniform_int_distribution <GLuint>   v_distribution(1, _v_geneNumber - 2);

    if(prob_distribution2(generator2) < probability)
    {
        GLuint parent = surface_distribution(generator2);
        SurfaceIndividual offspring(_individual[parent]);
        GLdouble range = 0;

        if(prob_distribution2(generator2) < 0.5)  //which chromosome to be mutated
        {
            GLuint geneToBeMutated = u_distribution(generator2);
            if(prob_distribution2(generator2) < 0.5)  //shift the gene to the left
            {
                range = offspring[geneToBeMutated - 1] - offspring[geneToBeMutated];
            }
            else                                    //shift to the right
            {
                range = offspring[geneToBeMutated + 1] - offspring[geneToBeMutated];
            }
            offspring[geneToBeMutated] += range * mutationRadiusPercentage * prob_distribution2(generator2);
        }
        else
        {
            GLuint geneToBeMutated = v_distribution(generator2);
            if(prob_distribution2(generator2) < 0.5)  //shift the gene to the left
            {
                range = offspring(geneToBeMutated - 1) - offspring(geneToBeMutated);
            }
            else                                    //shift to the right
            {
                range = offspring(geneToBeMutated + 1) - offspring(geneToBeMutated);
            }
            offspring(geneToBeMutated) += range * mutationRadiusPercentage * prob_distribution2(generator2);
        }

        offspring.CalculateFitness(_dataToInterpolate, _divPointCount, _divPointCount);
        if(offspring.Fitness() > _individual[parent].Fitness())
        {
            _individual[parent] =  offspring;
        }
    }
}

GLvoid SurfacePopulation::TournamentSelection(GLuint poolSize)
{
    GLuint individualNumber =   _individual.GetColumnCount();
    RowMatrix<GLuint>           preSelection(poolSize);
    RowMatrix<SurfaceIndividual>  winners(individualNumber);

    std::uniform_int_distribution<> distribution(0,individualNumber-1);

    for(GLuint i = 0; i < individualNumber; i++)
    {
        for(GLuint j = 0; j < poolSize; j++)
        {
            preSelection[j] = distribution(generator2);
        }
        winners[i] = _individual[FindBestIndividual(preSelection)];
    }
    swap(_individual, winners);
}

GLuint SurfacePopulation::FindBestIndividual(const RowMatrix<GLuint> &pool) const
{
    if(pool.GetColumnCount() < 1)
    {
        cout << "Empty pool!" << endl;
    }
    GLuint indexOfBestIndividual = pool[0];
    for(GLuint i = 0;i < pool.GetColumnCount();i++)
    {
        if(_individual[pool[i]].Fitness() > _individual[indexOfBestIndividual].Fitness())
        {
            indexOfBestIndividual = pool[i];
        }
    }
    return indexOfBestIndividual;
}

GLvoid SurfacePopulation::Evolve(GLdouble mutationProbability, GLdouble mutationRadiusPercentage,
              GLdouble recombinationProbability,
              GLuint poolSize)
{
    //TournamentSelection(poolSize);
    Mutation(mutationProbability,mutationRadiusPercentage);
    //Recombination(recombinationProbability);
}
