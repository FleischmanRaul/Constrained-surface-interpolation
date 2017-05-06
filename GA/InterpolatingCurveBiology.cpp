#include "GA/InterpolatingCurveBiology.h"
#include <vector>       // std::vector
#include <limits>

using namespace cagd;
using namespace std;

std::random_device rd;
std::mt19937 generator(rd());
std::uniform_real_distribution<GLdouble> prob_distribution(0.0,1.0);

//-------------------------- Implementation of CurveIndividual --------------------------
CurveIndividual::CurveIndividual(GLuint geneNumber)
{
    _chromosome.ResizeRows(geneNumber);
    _cc = nullptr;
}

//CurveIndividual::CurveIndividual(CurveIndividual& ci)
//_type(ci._type)
//{}

CurveIndividual::CurveIndividual(CurveType type,CurveEnergyType etype, GLuint geneNumber):
    _type(type),_etype(etype)
{
    _chromosome.ResizeRows(geneNumber);
    switch(_type)
    {
    case CYCLIC:
        _cc = new CyclicCurve3(geneNumber/2);
        break;
    case TRIGONOMETRIC:
        _cc = new TrigonometricCurve3(geneNumber/2, PI);
        break;
    default: std::cout << "Error in constructing " << _type << " curve!" << endl;
                 break;
    }
    _fitness = std::numeric_limits<GLdouble>::max();
}

//CurveIndividual::CurveIndividual(const CurveIndividual& ci):
//_type(ci._type),_fitness(ci._fitness),_etype(ci._etype),_chromosome(ci._chromosome)
//{
//    if(ci._cc != nullptr)
//    {
//        _cc = ci._cc->Clone();
//    }
//    else
//    {
//        _cc = nullptr;
//    }
//}

GLvoid CurveIndividual::SetEnergyProportions(const RowMatrix<GLdouble>& eProportion)
{
    _eProportion = eProportion;
}

Gene CurveIndividual::operator [](GLuint i) const
{
    return _chromosome[i];
}


Gene& CurveIndividual::operator [](GLuint i)
{
    return _chromosome[i];
}

GLdouble CurveIndividual::Fitness() const
{
    return _fitness;
}

GLboolean CurveIndividual::SetChromosome(ColumnMatrix<GLdouble> chromosome)
{
    if(_chromosome.GetRowCount() != chromosome.GetRowCount())
    {
        return GL_FALSE;
    }
    _chromosome = chromosome;
    return GL_TRUE;
}

GLboolean CurveIndividual::GenerateChromosome()
{
    GLuint size = _chromosome.GetRowCount();
    GLdouble u_min, u_max;
    GetDefinitionDomain(u_min,u_max);
    std::uniform_real_distribution<GLdouble> domain_distribution(u_min,u_max-0.0001);
    vector<GLdouble> myvector;
    for(GLuint i = 0; i<size;i++)
    {
        myvector.push_back(domain_distribution(generator));
    }
    std::sort(myvector.begin(),myvector.begin()+size);

    GLuint i = 0;
    for (std::vector<GLdouble>::iterator it=myvector.begin(); it!=myvector.end(); ++it)
    {
        _chromosome[i] = *it;
        i++;
    }
    if(_type == TRIGONOMETRIC)
    {
        _chromosome[0] = u_min;
        _chromosome[i-1] = u_max;
    }
    return GL_TRUE;
}

GLvoid CurveIndividual::GetDefinitionDomain(GLdouble& u_min, GLdouble& u_max) const
{
    _cc->GetDefinitionDomain(u_min,u_max);
}

GenericCurve3* CurveIndividual::GenerateImage()
{
    return _cc->GenerateImage(2,100);                   //maybe has to be canged to 1
}

GLboolean CurveIndividual::CalculateFitness(const ColumnMatrix<DCoordinate3> &dataToInterpolate, GLuint n)
{
    if(_chromosome.GetRowCount() != dataToInterpolate.GetRowCount() || _chromosome.GetRowCount() == 0)
    {
        return GL_FALSE;
    }


    if(!_cc->UpdateDataForInterpolation(_chromosome, dataToInterpolate))
    {
        return GL_FALSE;
    }
    switch(_etype)
    {
    case LENGTH:
        _eProportion.ResizeColumns(2);
        _eProportion[0] = 1;
        _eProportion[1] = 0;
        break;
    case CURVATURE:
        _eProportion.ResizeColumns(2);
        _eProportion[0] = 1;
        _eProportion[1] = 1;
        break;
    }

    _fitness = -_cc->Fitness(n, _eProportion);
    return GL_TRUE;
}

CurveIndividual::~CurveIndividual()
{

}



//-------------------------- Implementation of CurvePopulation --------------------------
CurvePopulation::CurvePopulation(GLuint individualCount, CurveType type, CurveEnergyType etype, const ColumnMatrix<DCoordinate3> &dataToInterpolate,
                                 GLuint maxMaturityLevel, GLuint divPointCount,RowMatrix<GLdouble>& eProportion, GLdouble threshold):
    _individual(individualCount),_dataToInterpolate(dataToInterpolate),_maxMaturityLevel(maxMaturityLevel),_divPointCount(divPointCount),_threshold(threshold)
{

    _geneNumber = dataToInterpolate.GetRowCount();
    for(GLuint i=0; i<individualCount; ++i)
    {
        _individual[i] = CurveIndividual(type,etype,_geneNumber);
        _individual[i].GenerateChromosome();
        _individual[i].SetEnergyProportions(eProportion);
        _individual[i].CalculateFitness(_dataToInterpolate,_divPointCount);
    }
}

GLuint CurvePopulation::GetMaxMaturityLevel()
{
    return _maxMaturityLevel;
}

GLdouble CurvePopulation::GetThreshold()
{
    return _threshold;
}

GLboolean CurvePopulation::SetDataToInterpolate(const ColumnMatrix<DCoordinate3> &dataToInterpolate)
{
    if(dataToInterpolate.GetRowCount() != _geneNumber)
    {
        return GL_FALSE;
    }
    _dataToInterpolate.ResizeRows(_geneNumber);
    _dataToInterpolate = dataToInterpolate;

//    for(GLuint i=0; i<individualCount; ++i)
//    {
//        _individual[i].CalculateFitness(_dataToInterpolate,_divPointCount);
//    }

    return GL_TRUE;
}

GenericCurve3* CurvePopulation::ImageOfIndividual(GLuint i) const
{
    return _individual[i].GenerateImage();
}

GenericCurve3* CurvePopulation::ImageOfBestIndividual()
{
    if(_individual.GetColumnCount() < 1)
    {
        throw Exception("Empty Curve Population!");
    }
    _indexOfBestIndividual = 0;
    for(GLuint i = 1; i<_individual.GetColumnCount(); ++i)
    {
        if(_individual[i].Fitness() > _individual[_indexOfBestIndividual].Fitness())
        {
            _indexOfBestIndividual = i;
        }
    }
    return  _individual[_indexOfBestIndividual].GenerateImage();
}

GLvoid CurvePopulation::TournamentSelection(GLuint poolSize)
{
    GLuint individualNumber =   _individual.GetColumnCount();
    RowMatrix<GLuint>           preSelection(poolSize);
    RowMatrix<CurveIndividual>  winners(individualNumber);

    std::uniform_int_distribution<> distribution(0,individualNumber-1);

    for(GLuint i = 0; i < individualNumber; i++)
    {
        for(GLuint j = 0; j < poolSize; j++)
        {
            preSelection[j] = distribution(generator);
        }
        winners[i] = _individual[FindBestIndividual(preSelection)];
    }
    swap(_individual, winners);
}

GLvoid CurvePopulation::Recombination(GLdouble probability)
{
    std::uniform_int_distribution <GLuint>   int_distribution(0, _individual.GetColumnCount() - 1);
    GLuint r,w;
    GLdouble u_min, u_max, value;
    _individual[0].GetDefinitionDomain(u_min,u_max);        //all individuals should have the same definition domain
    std::uniform_real_distribution<GLdouble> domain_distribution(u_min,u_max);
    if(prob_distribution(generator) < probability)
    {
        //cout << "In recombination" << endl;
        r = int_distribution(generator);
        w = int_distribution(generator);
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
        CurveIndividual rOffspring = _individual[r];
        CurveIndividual wOffspring = _individual[w];
        value = domain_distribution(generator);
        GLuint i = 0;

        while(rOffspring[i] < value && wOffspring[i] < value && i < _geneNumber-2)
        {
            swap(rOffspring[i],  wOffspring[i]);
            i++;
        }
        rOffspring.CalculateFitness(_dataToInterpolate, _divPointCount);
        wOffspring.CalculateFitness(_dataToInterpolate, _divPointCount);

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
        _individual[_indexOfBestIndividual].CalculateFitness(_dataToInterpolate, _divPointCount);
    }
}

GLvoid CurvePopulation::Mutation(GLdouble probability, GLdouble mutationRadiusPercentage)
{
    std::uniform_int_distribution <GLuint>   curve_distribution(0, _individual.GetColumnCount() - 1);
    std::uniform_int_distribution <GLuint>   gene_distribution(1, _geneNumber - 2);
    if(prob_distribution(generator) < probability)
    {
        //cout << "In mutation" << endl;
        GLuint parent = curve_distribution(generator);
        CurveIndividual offspring(_individual[parent]);
        GLuint geneToBeMutated = gene_distribution(generator);
        GLdouble range = 0;


        if(prob_distribution(generator) < 0.5)  //shift the gene to the left
        {
            range = offspring[geneToBeMutated - 1] - offspring[geneToBeMutated];
        }
        else                                    //shift to the right
        {
            range = offspring[geneToBeMutated + 1] - offspring[geneToBeMutated];
        }
        offspring[geneToBeMutated] += range * mutationRadiusPercentage * prob_distribution(generator);
        offspring.CalculateFitness(_dataToInterpolate, _divPointCount);
        //cout << offspring.Fitness() << "  " << _individual[parent].Fitness() << endl;
        if(offspring.Fitness() > _individual[parent].Fitness())
        {
            _individual[parent] =  offspring;
        }
        _individual[parent].CalculateFitness(_dataToInterpolate, _divPointCount);
    }
}

GLuint CurvePopulation::FindBestIndividual(const RowMatrix<GLuint> &pool) const
{
    if(pool.GetColumnCount() < 1)
    {
        throw Exception("Empty pool");
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
GLdouble CurvePopulation::FitnessOfBestIndividual()
{
    return _individual[_indexOfBestIndividual].Fitness();
}

GLboolean CurvePopulation::CalculateFitnesses()
{
    for(GLuint i = 0; i<_individual.GetColumnCount(); ++i)
    {
        if(!_individual[i].CalculateFitness(_dataToInterpolate, _divPointCount))
        {
            cout << "Could not calculate fitness!" << endl;
            return GL_FALSE;
        }
    }
    return GL_TRUE;
}

GLvoid CurvePopulation::Evolve(GLdouble mutationProbability, GLdouble mutationRadiusPercentage,
              GLdouble recombinationProbability,
              GLuint poolSize)
{
    TournamentSelection(poolSize);
    Mutation(mutationProbability,mutationRadiusPercentage);
    Recombination(recombinationProbability);

}
