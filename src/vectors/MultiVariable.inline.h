
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

  inline
  bool  variableLessThan ( const Variable::shared_ptr  &left , 
                           const Variable::shared_ptr  &right )
  {
    return left->getName() < right->getName();
  }

  inline
  bool  variableEquals ( const Variable::shared_ptr  &left , 
                           const Variable::shared_ptr  &right )
  {
    return ((*left) == (*right));
  }

  inline
  void MultiVariable::removeDuplicateVariables ()
  {
    std::sort ( beginVariable() , endVariable() , variableLessThan );
    std::vector<Variable::shared_ptr>::iterator new_end = std::unique ( beginVariable() , endVariable(), variableEquals );
    d_vVariables.resize ( new_end - beginVariable() );
  }

}
}

