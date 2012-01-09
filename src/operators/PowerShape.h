//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   operators/PowerShape.h
 * \author Kevin Clarno and Gokhan Yesilyurt
 * \brief  Power shape   
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: PowerShape.h,v 1.3 2010/06/15 12:00:00 h51 Exp $
//---------------------------------------------------------------------------//

#ifndef included_AMP_PowerShape_h 
#define included_AMP_PowerShape_h 

#include <vector>
#include "utils/Utilities.h"

/* AMP files */
#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "PowerShapeParameters.h"
#include "vectors/DualVariable.h"
#include "vectors/Variable.h"
#include "utils/InputDatabase.h"

/*Boost files */
#include "boost/shared_ptr.hpp"

/*LibMesh Includes */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature_gauss.h"

#include <vector>

namespace AMP {
namespace Operator {

//===========================================================================//
/*!
 * \class PowerShape
 * \brief Provides the specific power distribution by node for a given mesh.
 */
/*! 
 * \examples operators/test/testPowerShape-1.cc
 *           operators/test/testPowerShape-2.cc
 *           operators/test/testPowerShape-3.cc
 *           operators/test/testPowerShape-5.cc
 */
//===========================================================================//

class PowerShape : public  Operator {

    public:
      //typedef AMP::LinearAlgebra::VectorVariable<AMP::Mesh::IntegrationPointVariable, 8>      HexGaussPointVariable;
      //typedef boost::shared_ptr<HexGaussPointVariable>      SP_HexGaussPointVariable;
      typedef boost::shared_ptr<PowerShapeParameters>                  SP_Parameters;
      typedef boost::shared_ptr<OperatorParameters>            SP_OperatorParameters;
      typedef boost::shared_ptr<AMP::LinearAlgebra::Vector>                                    SP_Vector; 
      typedef std::vector<double>                                            Vec_Dbl;
      typedef boost::shared_ptr<Vec_Dbl>                                  SP_Vec_Dbl; 
      typedef boost::shared_ptr<AMP::Database>                           SP_Database;
    
  private:
    // Defines fission data types.
    enum PowerShape_Types
    {
        LINEAR,
        QUADRATIC,
        CUBIC,
        NUM_POWER_SHAPES
    };
    
    // use a spatially constant power distribution
    bool d_useFixedPower;

    // Coordinate system.
    std::string d_coordinateSystem;

    // Cylindrical coordinate system.
    std::string d_type;

    // Number of moments in the x direction.
    unsigned int d_numXmoments;
    
    // Moments in the x direction:  size= d_numXmoments.
    std::vector<double> d_Xmoments;

    // Number of moments in the y direction.
    unsigned int d_numYmoments;
    
    // Moments in the y direction:  size= d_numYmoments.
    std::vector<double> d_Ymoments;

    // Number of moments in the z direction.
    unsigned int d_numZmoments;
    
    // Moments in the z direction:  size= d_numZmoments.
    std::vector<double> d_Zmoments;

    // Number of moments in the Zernike basis function.
    unsigned int d_numMoments;
    unsigned int d_numMNmoments;
    
    // Moments in the Zernike basis function.
    std::vector<double> d_Moments;

    // Frapcon volume integral type.
    std::string d_frapconVolumeIntegral;

    // Frapcon constant.
    double d_frapconConstant;

    // Angular constant.
    double d_angularConstant;

    // Gaussian constants
    double d_muX;
    double d_muY;
    double d_sigmaX;
    double d_sigmaY;

  public:

      /*
       * A class for representing the neutronics source operator.
      */
      PowerShape(SP_Parameters parameters);

      /**
       * Empty destructor for PowerShape
       */
      virtual ~PowerShape();

      /**
       * Print out all members of integrator instance to given output stream.
       */
      void printClassData(std::ostream& os) const;

      /**
       * Write out state of object to given database.
       *
       * When assertion checking is active, the database pointer must be non-null.
       */
      void putToDatabase(SP_Database db);

      /**
        The function that computes the residual.
       * @param f: rhs vector for A(u)=f, this may be a null pointer if f=0. 
       * @param u: multivector of the state.
       * @param r: specific power in Watts per gram 
       The result of apply is
       * r = b*f+a*A(u)
       */
      void apply(const  SP_Vector & f, 
                 const  SP_Vector & u, 
                        SP_Vector & r,
                 const  double      a = 1.0,
                 const  double      b = 0.0);

      /**
        A function to reinitialize this object.
        */
      void reset(const SP_OperatorParameters & parameters);

      /*SP_HexGaussPointVariable createOutputVariable (const std::string & name) 
      {
        SP_HexGaussPointVariable var( new HexGaussPointVariable (name) );
        return var;
      }*/

      double evalFactorial(const int n);
      double choose(int, int);
      double evalZernike(const int m, const int n, const double rho, const double phi);
      double evalLegendre(const int n, const double x);
      double getVolumeIntegralAnalytical(double rmax);
      double getVolumeIntegralSum(double rmax, double cx, double cy);
      double getFrapconFr(double radius, double rmax);
      double getZernikeRadial(double rho);
      double getZernike(double rho, double phi);
      double getGaussianF(double x, double y);

    protected:

      /*
       * Read input data from specified database and initialize class members.
       * If run is from restart, a subset of the restart values may be replaced
       * with those read from input.
       *
       * When assertion checking is active, the database pointer must be non-null.
       */
      void getFromDatabase(SP_Database db);

      SP_Database d_db;

      //SP_HexGaussPointVariable d_Variable;

      boost::shared_ptr < ::FEType > d_feType;
      boost::shared_ptr < ::FEBase > d_fe;
      boost::shared_ptr < ::QBase >  d_qrule;

      AMP::Mesh::Mesh::shared_ptr d_Mesh;
    private:

  };

}
}

#endif
    
//---------------------------------------------------------------------------//
//              end of operator/PowerShape.h
//---------------------------------------------------------------------------//
