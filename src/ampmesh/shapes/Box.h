#ifndef included_AMP_Geometry_Box
#define included_AMP_Geometry_Box

#include "AMP/ampmesh/Geometry.h"

#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
template<std::size_t NDIM>
class Box : public Geometry
{
public:
    /**
     * \brief Construct a box
     * \param range     The range of the box [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     * \param coord     Optional list of coordinates to adjust the logical - physical mapping
     */
    explicit Box( const std::vector<double> &range );

public: // Functions inherited from Geometry
    virtual uint8_t getDim() const override final { return NDIM; }
    virtual double distance( const Point<double> &pos,
                             const Point<double> &dir ) const override final;
    virtual bool inside( const Point<double> &pos ) const override final;
    virtual int surface( const Point<double> &x ) const override final;
    virtual Point<double> surfaceNorm( const Point<double> &x ) const override final;
    virtual Point<double> logical( const Point<double> &x ) const override final;
    virtual Point<double> physical( const Point<double> &x ) const override final;
    virtual void displaceMesh( const double *x ) override final;

protected:
    // Internal data
    double d_range[6];

private:
    // Private constuctor
    Box();
};


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
template<std::size_t NDIM>
class Grid : public Geometry
{
public:
    /**
     * \brief Construct a grid geometry
     * \param range     The range of the box [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     * \param coord     Optional list of coordinates to adjust the logical - physical mapping
     */
    explicit Grid( const std::vector<std::vector<double>> &coord );

public: // Functions inherited from Geometry
    virtual uint8_t getDim() const override final { return NDIM; }
    virtual double distance( const Point<double> &pos,
                             const Point<double> &dir ) const override final;
    virtual bool inside( const Point<double> &pos ) const override final;
    virtual int surface( const Point<double> &x ) const override final;
    virtual Point<double> surfaceNorm( const Point<double> &x ) const override final;
    virtual Point<double> logical( const Point<double> &x ) const override final;
    virtual Point<double> physical( const Point<double> &x ) const override final;
    virtual void displaceMesh( const double *x ) override final;

protected:
    // Internal data
    double d_range[6];
    std::vector<double> d_coord[NDIM];

private:
    // Private constuctor
    Grid();
};


} // namespace Geometry
} // namespace AMP

#endif
