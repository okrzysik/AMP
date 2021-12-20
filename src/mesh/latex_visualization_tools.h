#ifndef LATEX_VISUALIZATION_TOOLS
#define LATEX_VISUALIZATION_TOOLS

#include "AMP/mesh/hex8_element_t.h"
#include "AMP/mesh/triangle_t.h"

#include <iostream>

void write_point( double const *p, std::ostream &os = std::cout );
void write_cycle( unsigned int n, double const **c, std::ostream &os = std::cout );
void write_face( double const **f, std::ostream &os = std::cout );
void write_triangle( double const **t, std::ostream &os = std::cout );

void draw_point( double const *p,
                 const std::string &option = "",
                 std::ostream &os          = std::cout,
                 const std::string &text   = "." );
void draw_line( double const *b,
                double const *e,
                const std::string &option = "",
                std::ostream &os          = std::cout );
void draw_line( unsigned int n,
                double const *l,
                const std::string &option = "",
                std::ostream &os          = std::cout );
void draw_triangle( triangle_t *t_ptr,
                    const std::string &option = "",
                    std::ostream &os          = std::cout );
void draw_face( hex8_element_t *e_ptr,
                unsigned int f,
                const std::string &option = "",
                std::ostream &os          = std::cout );


void draw_bounding_box( hex8_element_t *e_ptr,
                        double const *point_of_view,
                        std::ostream &os = std::cout );
void draw_bounding_polyhedron( hex8_element_t *e_ptr,
                               double const *point_of_view,
                               std::ostream &os = std::cout );
void draw_hex8_element( hex8_element_t *e_ptr,
                        double const *point_of_view,
                        std::ostream &os = std::cout );

#endif // LATEX_VISUALIZATION_TOOLS
