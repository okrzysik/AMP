#ifndef included_AMP_test_Environments
#define included_AMP_test_Environments

namespace  AMP {
namespace  unit_test {

// Not technically an environment class, this test is handy to have
// around to make sure unit test implementations work.
class  AllPassTest
{
   public:
      static const char * get_test_name () { return "test instantiation"; }

      template <typename UTILS>
      static  void run_test ( UTILS &utils )
      {
         utils.pass_test ();
      }
}
;



class  AllEnvironments
{
   public:
      static const char * get_env_name () { return "all environments"; }

      template <typename UTILS>
      static  bool verify_environment ( UTILS & ) { return true; }
}
;


class  SerialTest
{
   public:
      static const char * get_env_name () { return "serial only"; }

      template <typename UTILS>
      static  bool verify_environment ( UTILS &utils )
      {
         return utils.comm_size() == 1;
      }
}
;

class  ParallelTest
{
  public:
    static const char * get_env_name () { return "multiprocessor only"; }

    template <typename UTILS>
    static  bool verify_environment ( UTILS &utils )
    {
      return utils.comm_size() > 1;
    }
}
;

class  MeshHasBCs
{
  public:
    static const char * get_env_name () { return "mesh has boundary conditions"; }

    template <typename UTILS>
    static  bool verify_environment ( UTILS &utils )
    {
      if ( utils.mesh )
        return utils.mesh->getBoundaryIds().size() > 0;
      return true;  // Don't know yet....
    }
};


template <typename FIRST_CONDITION , typename SECOND_CONDITION>
class AndCombinedEnvironment
{
  public:
    static std::string get_env_name ()
    {
      std::stringstream  ans;
      ans << FIRST_CONDITION::get_env_name () << " AND " << SECOND_CONDITION::get_env_name();
      return ans.str();
    }

    template <typename UTILS>
    static  bool verify_environment ( UTILS &utils )
    {
      return FIRST_CONDITION::verify_environment ( utils ) &&
             SECOND_CONDITION::verify_environment ( utils );
    }
};

template <typename FIRST_CONDITION , typename SECOND_CONDITION>
class OrCombinedEnvironment
{
  public:
    static std::string get_env_name ()
    {
      std::stringstream  ans;
      ans << FIRST_CONDITION::get_env_name () << " OR " << SECOND_CONDITION::get_env_name();
      return ans.str();
    }

    template <typename UTILS>
    static  bool verify_environment ( UTILS &utils )
    {
      return FIRST_CONDITION::verify_environment ( utils ) ||
             SECOND_CONDITION::verify_environment ( utils );
    }
};

template <typename CONDITION>
class NotEnvironment
{
  public:
    static std::string get_env_name ()
    {
      std::stringstream  ans;
      ans << "NOT " << CONDITION::get_env_name();
      return ans.str();
    }

    template <typename UTILS>
    static  bool verify_environment ( UTILS &utils )
    {
      return !CONDITION::verify_environment ( utils );
    }
};



}
}


#endif
