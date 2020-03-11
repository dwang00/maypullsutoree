#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"


/*======== void add_circle() ==========
  Inputs:   struct matrix * edges
            double cx
            double cy
            double r
            double step

  Adds the circle at (cx, cy) with radius r to edges
  ====================*/
  void add_circle( struct matrix *edges,
                   double cx, double cy, double cz,
                   double r, double step ) {
    double a = -1;
    double b = -1;
    for (double i = 0; i < 1 + step; i += step){
      double c =  r * cos(2 * M_PI * i) + cx;
      double d =  r * sin(2 * M_PI * i) + cy;
      if (a != -1){
        add_edge(edges, a, b, cz, c, d, cz);
      }
      a = c;
      b = d;
    }
  }

  /*======== void add_curve() ==========
  Inputs:   struct matrix *edges
           double x0
           double y0
           double x1
           double y1
           double x2
           double y2
           double x3
           double y3
           double step
           int type

  Adds the curve bounded by the 4 points passsed as parameters
  of type specified in type (see matrix.h for curve type constants)
  to the matrix edges
  ====================*/
  void add_curve( struct matrix *edges,
                  double x0, double y0,
                  double x1, double y1,
                  double x2, double y2,
                  double x3, double y3,
                  double step, int type ) {
    struct matrix *xval = generate_curve_coefs(x0, x1, x2, x3, type);
    double aX = xval->m[0][0];
    double bX = xval->m[1][0];
    double cX = xval->m[2][0];
    double dX = xval->m[3][0];
    struct matrix *yval = generate_curve_coefs(y0, y1, y2, y3, type);
    double aY = yval->m[0][0];
    double bY = yval->m[1][0];
    double cY = yval->m[2][0];
    double dY = yval->m[3][0];

    double a = -1;
    double b = -1;
    for (double i = 0; i <= 1; i += step){
      double c =  aX * pow(i, 3.0) + bX * pow(i, 2.0) + cX * i + dX;
      double d =  aY * pow(i, 3.0) + bY * pow(i, 2.0) + cY * i + dY;
      if (a != -1){
        add_edge(edges, a, b, 0, c, d, 0);
      }
      a = c;
      b = d;
    }
  }


/*======== void add_point() ==========
Inputs:   struct matrix * points
         int x
         int y
         int z
Returns:
adds point (x, y, z) to points and increment points.lastcol
if points is full, should call grow on points
====================*/
void add_point( struct matrix * points, double x, double y, double z) {

  if ( points->lastcol == points->cols )
    grow_matrix( points, points->lastcol + 100 );

  points->m[0][ points->lastcol ] = x;
  points->m[1][ points->lastcol ] = y;
  points->m[2][ points->lastcol ] = z;
  points->m[3][ points->lastcol ] = 1;
  points->lastcol++;
} //end add_point

/*======== void add_edge() ==========
Inputs:   struct matrix * points
          int x0, int y0, int z0, int x1, int y1, int z1
Returns:
add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
should use add_point
====================*/
void add_edge( struct matrix * points,
	       double x0, double y0, double z0,
	       double x1, double y1, double z1) {
  add_point( points, x0, y0, z0 );
  add_point( points, x1, y1, z1 );
}

/*======== void draw_lines() ==========
Inputs:   struct matrix * points
         screen s
         color c
Returns:
Go through points 2 at a time and call draw_line to add that line
to the screen
====================*/
void draw_lines( struct matrix * points, screen s, color c) {

 if ( points->lastcol < 2 ) {
   printf("Need at least 2 points to draw a line!\n");
   return;
 }

 int point;
 for (point=0; point < points->lastcol-1; point+=2)
   draw_line( points->m[0][point],
	      points->m[1][point],
	      points->m[0][point+1],
	      points->m[1][point+1],
	      s, c);
}// end draw_lines









void draw_line(int x0, int y0, int x1, int y1, screen s, color c) {

  int x, y, d, A, B;
  //swap points if going right -> left
  int xt, yt;
  if (x0 > x1) {
    xt = x0;
    yt = y0;
    x0 = x1;
    y0 = y1;
    x1 = xt;
    y1 = yt;
  }

  x = x0;
  y = y0;
  A = 2 * (y1 - y0);
  B = -2 * (x1 - x0);

  //octants 1 and 8
  if ( abs(x1 - x0) >= abs(y1 - y0) ) {

    //octant 1
    if ( A > 0 ) {

      d = A + B/2;
      while ( x < x1 ) {
	plot( s, c, x, y );
	if ( d > 0 ) {
	  y+= 1;
	  d+= B;
	}
	x++;
	d+= A;
      } //end octant 1 while
      plot( s, c, x1, y1 );
    } //end octant 1

    //octant 8
    else {
      d = A - B/2;

      while ( x < x1 ) {
	//printf("(%d, %d)\n", x, y);
	plot( s, c, x, y );
	if ( d < 0 ) {
	  y-= 1;
	  d-= B;
	}
	x++;
	d+= A;
      } //end octant 8 while
      plot( s, c, x1, y1 );
    } //end octant 8
  }//end octants 1 and 8

  //octants 2 and 7
  else {

    //octant 2
    if ( A > 0 ) {
      d = A/2 + B;

      while ( y < y1 ) {
	plot( s, c, x, y );
	if ( d < 0 ) {
	  x+= 1;
	  d+= A;
	}
	y++;
	d+= B;
      } //end octant 2 while
      plot( s, c, x1, y1 );
    } //end octant 2

    //octant 7
    else {
      d = A/2 - B;

      while ( y > y1 ) {
	plot( s, c, x, y );
	if ( d > 0 ) {
	  x+= 1;
	  d+= A;
	}
	y--;
	d-= B;
      } //end octant 7 while
      plot( s, c, x1, y1 );
    } //end octant 7
  }//end octants 2 and 7
} //end draw_line
