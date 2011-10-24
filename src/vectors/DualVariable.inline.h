
namespace AMP {
namespace LinearAlgebra {

  inline
  DualVariable::DualVariable ( Variable::shared_ptr one , 
                               Variable::shared_ptr two , 
                               std::string name ) 
                             : Variable ( name ) , 
                               d_Var1 ( one ) , 
                               d_Var2 ( two ) 
  {
  }

  inline
  DualVariable::~DualVariable () 
  {
  }

  inline
  Variable::shared_ptr  DualVariable::first () 
  { 
    return d_Var1; 
  }

  inline
  Variable::shared_ptr  DualVariable::second () 
  { 
    return d_Var2; 
  }

  inline
  Variable::shared_ptr  DualVariable::cloneVariable ( const std::string &name ) const
  {
    return Variable::shared_ptr ( new DualVariable ( d_Var1 , d_Var2 , name ) );
  }

}
}
