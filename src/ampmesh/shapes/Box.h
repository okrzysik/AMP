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
    virtual std::string getName() const override final;
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int surface( const Point &x ) const override final;
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point logical( const Point &x ) const override final;
    virtual Point physical( const Point &x ) const override final;
    virtual Point centroid() const override final;
    virtual std::pair<Point, Point> box() const override final;
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
    virtual std::string getName() const override final;
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int surface( const Point &x ) const override final;
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point logical( const Point &x ) const override final;
    virtual Point physical( const Point &x ) const override final;
    virtual Point centroid() const override final;
    virtual std::pair<Point, Point> box() const override final;
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
