#include "AMP/operators/mechanics/MechanicsManufacturedSolutions.h"
#include "AMP/utils/PIO.h"


namespace AMP {
namespace MechanicsManufacturedSolution {


// Default constructor
MMS::MMS( double YM, double PR )
    : E( YM ),
      nu( PR ),
      Lx( 1.0 ),
      Ly( 1.0 ),
      Lz( 1.0 ),
      Mx( 1.0 ),
      My( 1.0 ),
      Mz( 1.0 ),
      axx( 0.0 ),
      bxx( 1.0 ),
      axy( 0.0 ),
      bxy( 1.0 ),
      axz( 0.0 ),
      bxz( 1.0 ),
      ayx( 0.0 ),
      byx( 1.0 ),
      ayy( 0.0 ),
      byy( 1.0 ),
      ayz( 0.0 ),
      byz( 1.0 ),
      azx( 0.0 ),
      bzx( 1.0 ),
      azy( 0.0 ),
      bzy( 1.0 ),
      azz( 0.0 ),
      bzz( 1.0 ),
      name( "unnamed" )
{
}


// Copy constructor
MMS::MMS( const MMS &mms ) : name( mms.name )
{
    E   = mms.E;
    nu  = mms.nu;
    Lx  = mms.Lx;
    Ly  = mms.Ly;
    Lz  = mms.Lz;
    Mx  = mms.Mx;
    My  = mms.My;
    Mz  = mms.Mz;
    axx = mms.axx;
    bxx = mms.bxx;
    axy = mms.axy;
    bxy = mms.bxy;
    axz = mms.axz;
    bxz = mms.bxz;
    ayx = mms.ayx;
    byx = mms.byx;
    ayy = mms.ayy;
    byy = mms.byy;
    ayz = mms.ayz;
    byz = mms.byz;
    azx = mms.azx;
    bzx = mms.bzx;
    azy = mms.azy;
    bzy = mms.bzy;
    azz = mms.azz;
    bzz = mms.bzz;
}


// Assignment operator
MMS &MMS::operator=( const MMS &mms )
{
    if ( this == &mms )
        return *this;
    E    = mms.E;
    nu   = mms.nu;
    Lx   = mms.Lx;
    Ly   = mms.Ly;
    Lz   = mms.Lz;
    Mx   = mms.Mx;
    My   = mms.My;
    Mz   = mms.Mz;
    axx  = mms.axx;
    bxx  = mms.bxx;
    axy  = mms.axy;
    bxy  = mms.bxy;
    axz  = mms.axz;
    bxz  = mms.bxz;
    ayx  = mms.ayx;
    byx  = mms.byx;
    ayy  = mms.ayy;
    byy  = mms.byy;
    ayz  = mms.ayz;
    byz  = mms.byz;
    azx  = mms.azx;
    bzx  = mms.bzx;
    azy  = mms.azy;
    bzy  = mms.bzy;
    azz  = mms.azz;
    bzz  = mms.bzz;
    name = mms.name;
    return *this;
}


// MMSLinear
double MMSLinear::getExactSolutionX( double x, double y, double z ) const
{
    return Mx * ( bxx + ( axx * x ) / Lx ) * ( bxy + ( axy * y ) / Ly ) *
           ( bxz + ( axz * z ) / Lz );
}
double MMSLinear::getExactSolutionY( double x, double y, double z ) const
{
    return My * ( byx + ( ayx * x ) / Lx ) * ( byy + ( ayy * y ) / Ly ) *
           ( byz + ( ayz * z ) / Lz );
}
double MMSLinear::getExactSolutionZ( double x, double y, double z ) const
{
    return Mz * ( bzx + ( azx * x ) / Lx ) * ( bzy + ( azy * y ) / Ly ) *
           ( bzz + ( azz * z ) / Lz );
}
double MMSLinear::getForcingTermX( double /* x */, double y, double z ) const
{
    return ( ( Mz * azx * azz * ( bzy + ( azy * y ) / Ly ) ) / ( Lx * Lz ) +
             ( My * ayx * ayy * ( byz + ( ayz * z ) / Lz ) ) / ( Lx * Ly ) ) *
               ( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) -
           ( E * Mz * azx * azz * ( bzy + ( azy * y ) / Ly ) ) / ( Lx * Lz * ( nu * 2.0 + 2.0 ) ) -
           ( E * My * ayx * ayy * ( byz + ( ayz * z ) / Lz ) ) / ( Lx * Ly * ( nu * 2.0 + 2.0 ) );
}
double MMSLinear::getForcingTermY( double x, double /* y */, double z ) const
{
    return ( ( Mz * azy * azz * ( bzx + ( azx * x ) / Lx ) ) / ( Ly * Lz ) +
             ( Mx * axx * axy * ( bxz + ( axz * z ) / Lz ) ) / ( Lx * Ly ) ) *
               ( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) -
           ( E * Mz * azy * azz * ( bzx + ( azx * x ) / Lx ) ) / ( Ly * Lz * ( nu * 2.0 + 2.0 ) ) -
           ( E * Mx * axx * axy * ( bxz + ( axz * z ) / Lz ) ) / ( Lx * Ly * ( nu * 2.0 + 2.0 ) );
}
double MMSLinear::getForcingTermZ( double x, double y, double /* z */ ) const
{
    return ( ( My * ayy * ayz * ( byx + ( ayx * x ) / Lx ) ) / ( Ly * Lz ) +
             ( Mx * axx * axz * ( bxy + ( axy * y ) / Ly ) ) / ( Lx * Lz ) ) *
               ( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) -
           ( E * My * ayy * ayz * ( byx + ( ayx * x ) / Lx ) ) / ( Ly * Lz * ( nu * 2.0 + 2.0 ) ) -
           ( E * Mx * axx * axz * ( bxy + ( axy * y ) / Ly ) ) / ( Lx * Lz * ( nu * 2.0 + 2.0 ) );
}
std::vector<double> MMSLinear::getStressTensor( double x, double y, double z ) const
{
    std::vector<double> stress;
    stress.push_back(
        -( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) *
            ( ( Mz * azz * ( bzx + ( azx * x ) / Lx ) * ( bzy + ( azy * y ) / Ly ) ) / Lz +
              ( My * ayy * ( byx + ( ayx * x ) / Lx ) * ( byz + ( ayz * z ) / Lz ) ) / Ly +
              ( Mx * axx * ( bxy + ( axy * y ) / Ly ) * ( bxz + ( axz * z ) / Lz ) ) / Lx ) +
        ( E * Mx * axx * ( bxy + ( axy * y ) / Ly ) * ( bxz + ( axz * z ) / Lz ) * 2.0 ) /
            ( Lx * ( nu * 2.0 + 2.0 ) ) );
    stress.push_back(
        -( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) *
            ( ( Mz * azz * ( bzx + ( azx * x ) / Lx ) * ( bzy + ( azy * y ) / Ly ) ) / Lz +
              ( My * ayy * ( byx + ( ayx * x ) / Lx ) * ( byz + ( ayz * z ) / Lz ) ) / Ly +
              ( Mx * axx * ( bxy + ( axy * y ) / Ly ) * ( bxz + ( axz * z ) / Lz ) ) / Lx ) +
        ( E * My * ayy * ( byx + ( ayx * x ) / Lx ) * ( byz + ( ayz * z ) / Lz ) * 2.0 ) /
            ( Ly * ( nu * 2.0 + 2.0 ) ) );
    stress.push_back(
        -( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) *
            ( ( Mz * azz * ( bzx + ( azx * x ) / Lx ) * ( bzy + ( azy * y ) / Ly ) ) / Lz +
              ( My * ayy * ( byx + ( ayx * x ) / Lx ) * ( byz + ( ayz * z ) / Lz ) ) / Ly +
              ( Mx * axx * ( bxy + ( axy * y ) / Ly ) * ( bxz + ( axz * z ) / Lz ) ) / Lx ) +
        ( E * Mz * azz * ( bzx + ( azx * x ) / Lx ) * ( bzy + ( azy * y ) / Ly ) * 2.0 ) /
            ( Lz * ( nu * 2.0 + 2.0 ) ) );
    stress.push_back(
        ( E *
          ( ( My * ayz * ( byx + ( ayx * x ) / Lx ) * ( byy + ( ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) /
                Lz +
            ( Mz * azy * ( bzx + ( azx * x ) / Lx ) * ( bzz + ( azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) /
                Ly ) *
          2.0 ) /
        ( nu * 2.0 + 2.0 ) );
    stress.push_back(
        ( E *
          ( ( Mx * axz * ( bxx + ( axx * x ) / Lx ) * ( bxy + ( axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) /
                Lz +
            ( Mz * azx * ( bzy + ( azy * y ) / Ly ) * ( bzz + ( azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) /
                Lx ) *
          2.0 ) /
        ( nu * 2.0 + 2.0 ) );
    stress.push_back(
        ( E *
          ( ( Mx * axy * ( bxx + ( axx * x ) / Lx ) * ( bxz + ( axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) /
                Ly +
            ( My * ayx * ( byy + ( ayy * y ) / Ly ) * ( byz + ( ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) /
                Lx ) *
          2.0 ) /
        ( nu * 2.0 + 2.0 ) );
    return stress;
}


// MMSTrigonometric
double MMSTrigonometric::getExactSolutionX( double x, double y, double z ) const
{
    return Mx * sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
           sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
           sin( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) );
}
double MMSTrigonometric::getExactSolutionY( double /* x */, double y, double z ) const
{
    return My * sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
           sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
           sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) );
}
double MMSTrigonometric::getExactSolutionZ( double /* x */, double y, double z ) const
{
    return Mz * sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
           sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
           sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) );
}
double MMSTrigonometric::getForcingTermX( double x, double y, double z ) const
{
    return 1 / ( Lx * Lx ) * ( M_PI * M_PI ) * Mx * ( axx * axx ) *
               sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
               sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
               sin( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) *
               ( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) *
               ( -1.0 / 4.0 ) +
           ( E * 1 / ( Lx * Lx ) * ( M_PI * M_PI ) * Mx * ( axx * axx ) *
             sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
             sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
             sin( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
               ( nu * 2.0 + 2.0 ) +
           ( E * 1 / ( Ly * Ly ) * ( M_PI * M_PI * M_PI * M_PI ) * Mx * ( axy * axy ) *
             sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
             sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
             sin( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
               ( nu * 2.0 + 2.0 ) +
           ( E * 1 / ( Lz * Lz ) * ( M_PI * M_PI * M_PI * M_PI ) * Mx * ( axz * axz ) *
             sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
             sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
             sin( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
               ( nu * 2.0 + 2.0 );
}
double MMSTrigonometric::getForcingTermY( double x, double y, double z ) const
{
    return ( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) *
               ( 1 / ( Ly * Ly ) * ( M_PI * M_PI ) * My * ( ayx * ayx ) *
                     sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                     sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                     sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( -1.0 / 4.0 ) -
                 1 / ( Ly * Ly ) * ( M_PI * M_PI * M_PI * M_PI ) * My * ( ayy * ayy ) *
                     sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                     sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                     sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) +
                 1 / ( Ly * Ly ) * ( M_PI * M_PI * M_PI ) * My * ayx * ayy *
                     cos( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                     sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                     cos( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) +
                 ( ( M_PI * M_PI * M_PI ) * Mx * axx * axy *
                   cos( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   cos( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
                     ( Lx * Ly ) +
                 ( ( M_PI * M_PI * M_PI ) * Mz * azx * azy *
                   cos( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   cos( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
                     ( Ly * Lz ) +
                 ( ( M_PI * M_PI * M_PI * M_PI ) * Mz * azy * azz *
                   cos( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   cos( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
                     ( Ly * Lz ) ) +
           ( E *
             ( 1 / ( Ly * Ly ) * ( M_PI * M_PI ) * My * ( ayx * ayx ) *
                   sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) +
               1 / ( Ly * Ly ) * ( M_PI * M_PI * M_PI * M_PI ) * My * ( ayy * ayy ) *
                   sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) -
               1 / ( Ly * Ly ) * ( M_PI * M_PI * M_PI ) * My * ayx * ayy *
                   cos( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   cos( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) *
             2.0 ) /
               ( nu * 2.0 + 2.0 ) -
           ( E *
             ( 1 / ( Lz * Lz ) * ( M_PI * M_PI * M_PI * M_PI ) * My * ( ayz * ayz ) *
                   sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( -1.0 / 8.0 ) +
               ( ( M_PI * M_PI * M_PI ) * Mz * azx * azy *
                 cos( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                 sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                 cos( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 8.0 ) ) /
                   ( Ly * Lz ) +
               ( ( M_PI * M_PI * M_PI * M_PI ) * Mz * azy * azz *
                 cos( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                 cos( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                 sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 8.0 ) ) /
                   ( Ly * Lz ) ) *
             2.0 ) /
               ( nu * 2.0 + 2.0 ) -
           ( E * ( M_PI * M_PI * M_PI ) * Mx * axx * axy *
             cos( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
             sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
             cos( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
               ( Lx * Ly * ( nu * 2.0 + 2.0 ) );
}
double MMSTrigonometric::getForcingTermZ( double x, double y, double z ) const
{
    return ( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) *
               ( 1 / ( Lz * Lz ) * ( M_PI * M_PI ) * Mz * ( azx * azx ) *
                     sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                     sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                     sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( -1.0 / 4.0 ) -
                 1 / ( Lz * Lz ) * ( M_PI * M_PI * M_PI * M_PI ) * Mz * ( azz * azz ) *
                     sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                     sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                     sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) +
                 1 / ( Lz * Lz ) * ( M_PI * M_PI * M_PI ) * Mz * azx * azz *
                     cos( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                     sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                     cos( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) +
                 ( ( M_PI * M_PI * M_PI ) * Mx * axx * axz *
                   cos( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   cos( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
                     ( Lx * Lz ) +
                 ( ( M_PI * M_PI * M_PI ) * My * ayx * ayz *
                   cos( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   cos( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
                     ( Ly * Lz ) +
                 ( ( M_PI * M_PI * M_PI * M_PI ) * My * ayy * ayz *
                   cos( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   cos( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
                     ( Ly * Lz ) ) +
           ( E *
             ( 1 / ( Lz * Lz ) * ( M_PI * M_PI ) * Mz * ( azx * azx ) *
                   sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) +
               1 / ( Lz * Lz ) * ( M_PI * M_PI * M_PI * M_PI ) * Mz * ( azz * azz ) *
                   sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) -
               1 / ( Lz * Lz ) * ( M_PI * M_PI * M_PI ) * Mz * azx * azz *
                   cos( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   cos( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) *
             2.0 ) /
               ( nu * 2.0 + 2.0 ) -
           ( E *
             ( 1 / ( Ly * Ly ) * ( M_PI * M_PI * M_PI * M_PI ) * Mz * ( azy * azy ) *
                   sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                   sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( -1.0 / 8.0 ) +
               ( ( M_PI * M_PI * M_PI ) * My * ayx * ayz *
                 cos( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                 sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                 cos( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 8.0 ) ) /
                   ( Ly * Lz ) +
               ( ( M_PI * M_PI * M_PI * M_PI ) * My * ayy * ayz *
                 cos( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                 cos( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                 sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 8.0 ) ) /
                   ( Ly * Lz ) ) *
             2.0 ) /
               ( nu * 2.0 + 2.0 ) -
           ( E * ( M_PI * M_PI * M_PI ) * Mx * axx * axz *
             cos( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
             sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
             cos( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
               ( Lx * Lz * ( nu * 2.0 + 2.0 ) );
}
std::vector<double> MMSTrigonometric::getStressTensor( double x, double y, double z ) const
{
    std::vector<double> stress;
    stress.push_back(
        -( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) *
            ( ( ( M_PI * M_PI ) * My * ayy *
                cos( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Ly +
              ( ( M_PI * M_PI ) * Mz * azz *
                cos( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Lz +
              ( M_PI * Mx * axx * sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                cos( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Lx +
              ( M_PI * My * ayx * sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                cos( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Ly +
              ( M_PI * Mz * azx * sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                cos( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Lz ) +
        ( E * M_PI * Mx * axx * sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
          sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
          cos( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) ) /
            ( Lx * ( nu * 2.0 + 2.0 ) ) );
    stress.push_back(
        -( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) *
            ( ( ( M_PI * M_PI ) * My * ayy *
                cos( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Ly +
              ( ( M_PI * M_PI ) * Mz * azz *
                cos( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Lz +
              ( M_PI * Mx * axx * sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                cos( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Lx +
              ( M_PI * My * ayx * sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                cos( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Ly +
              ( M_PI * Mz * azx * sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                cos( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Lz ) +
        ( E *
          ( ( ( M_PI * M_PI ) * My * ayy *
              cos( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
              sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
              sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                Ly +
            ( M_PI * My * ayx * sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
              sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
              cos( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                Ly ) *
          2.0 ) /
            ( nu * 2.0 + 2.0 ) );
    stress.push_back(
        -( ( E * ( 2.0 / 3.0 ) ) / ( nu * 2.0 + 2.0 ) + E / ( nu * 6.0 - 3.0 ) ) *
            ( ( ( M_PI * M_PI ) * My * ayy *
                cos( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Ly +
              ( ( M_PI * M_PI ) * Mz * azz *
                cos( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Lz +
              ( M_PI * Mx * axx * sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                cos( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Lx +
              ( M_PI * My * ayx * sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                cos( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Ly +
              ( M_PI * Mz * azx * sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                cos( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                  Lz ) +
        ( E *
          ( ( ( M_PI * M_PI ) * Mz * azz *
              cos( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
              sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
              sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                Lz +
            ( M_PI * Mz * azx * sin( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
              sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
              cos( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                Lz ) *
          2.0 ) /
            ( nu * 2.0 + 2.0 ) );
    stress.push_back(
        ( E *
          ( ( ( M_PI * M_PI ) * My * ayz *
              cos( M_PI * ( byz + ( M_PI * ayz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
              sin( M_PI * ( byy + ( M_PI * ayy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
              sin( M_PI * ( byx + ( ayx * y ) / Ly ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
                Lz +
            ( ( M_PI * M_PI ) * Mz * azy *
              cos( M_PI * ( bzy + ( M_PI * azy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
              sin( M_PI * ( bzz + ( M_PI * azz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
              sin( M_PI * ( bzx + ( azx * z ) / Lz ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 4.0 ) ) /
                Ly ) *
          2.0 ) /
        ( nu * 2.0 + 2.0 ) );
    stress.push_back( ( E * ( M_PI * M_PI ) * Mx * axz *
                        cos( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                        sin( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                        sin( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                      ( Lz * ( nu * 2.0 + 2.0 ) ) );
    stress.push_back( ( E * ( M_PI * M_PI ) * Mx * axy *
                        cos( M_PI * ( bxy + ( M_PI * axy * y ) / Ly ) * ( 1.0 / 2.0 ) ) *
                        sin( M_PI * ( bxz + ( M_PI * axz * z ) / Lz ) * ( 1.0 / 2.0 ) ) *
                        sin( M_PI * ( bxx + ( axx * x ) / Lx ) * ( 1.0 / 2.0 ) ) * ( 1.0 / 2.0 ) ) /
                      ( Ly * ( nu * 2.0 + 2.0 ) ) );
    return stress;
}


// MMSBuilder
std::shared_ptr<MMS> MMSBuilder::createMMS( std::shared_ptr<AMP::Database> mmsDatabase )
{
    std::shared_ptr<MMS> mms;
    std::string name = mmsDatabase->getWithDefault<std::string>( "name", "One" );
    if ( name == "Trigonometric" ) {
        mms = std::shared_ptr<MMS>( new MMSTrigonometric );
    } else if ( name == "Linear" ) {
        mms = std::shared_ptr<MMS>( new MMSLinear );
    } else {
        std::cerr << "Error: AMP::MechanicsManufacturedSolution::MMS" << name << " is not defined"
                  << std::endl;
        AMP_ASSERT( false );
    }
    if ( mmsDatabase->keyExists( "scale_x" ) ) {
        mms->scaleX( mmsDatabase->getScalar<double>( "scale_x" ) );
    }
    if ( mmsDatabase->keyExists( "scale_y" ) ) {
        mms->scaleY( mmsDatabase->getScalar<double>( "scale_y" ) );
    }
    if ( mmsDatabase->keyExists( "scale_z" ) ) {
        mms->scaleZ( mmsDatabase->getScalar<double>( "scale_z" ) );
    }
    if ( mmsDatabase->keyExists( "scale_xyz" ) ) {
        mms->scaleXYZ( mmsDatabase->getScalar<double>( "scale_xyz" ) );
    }
    if ( mmsDatabase->keyExists( "mult_x" ) ) {
        mms->multX( mmsDatabase->getScalar<double>( "mult_x" ) );
    }
    if ( mmsDatabase->keyExists( "mult_y" ) ) {
        mms->multY( mmsDatabase->getScalar<double>( "mult_y" ) );
    }
    if ( mmsDatabase->keyExists( "mult_z" ) ) {
        mms->multZ( mmsDatabase->getScalar<double>( "mult_z" ) );
    }
    if ( mmsDatabase->keyExists( "mult_xyz" ) ) {
        mms->multXYZ( mmsDatabase->getScalar<double>( "mult_xyz" ) );
    }
    return mms;
}


} // end namespace MechanicsManufacturedSolution
} // end namespace AMP
