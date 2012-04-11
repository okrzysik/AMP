namespace AMP {
namespace LinearAlgebra {


inline
void  MultiVariable::setUnits ( const std::string &units )
{
    Variable::setUnits ( units );
    iterator  curVar = beginVariable();
    while ( curVar != endVariable() )
    {
        (*curVar)->setUnits ( units );
        curVar++;
    }
}


inline
MultiVariable::iterator  MultiVariable::beginVariable()
{
    return d_vVariables.begin();
}


inline
MultiVariable::iterator  MultiVariable::endVariable()
{
    return d_vVariables.end();
}


inline
MultiVariable::const_iterator  MultiVariable::beginVariable() const
{
    return d_vVariables.begin();
}


inline
MultiVariable::const_iterator  MultiVariable::endVariable() const
{
    return d_vVariables.end();
}


}
}

