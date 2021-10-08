/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the code you will need to replace for Lab 1.
 *
 * It contains two simple implementations that loop over the longest axis adding the gradient to the position on the other axis at each step.
 * The objective is for you to remove the use of floating point variables and work entirely with integers.
 * You should use Bresenhams Line Algorithm to achieve this.
 * 
 * e is the diff between 2 points, the midpoint between the 2 possible y points and the current point
 * 
 * dir means direction so 1 is positive gradient, -1 is negative gradient
 * 
 * sx and sy are the starting x and y points for the line - now x0 and y0
 * 
 * slope -> m (gradient)
 */

#include <iostream>
#include "linedrawer.h"

int draw_x_line(FrameBuffer *fb, int x0, int y0, int x1, int y1)
{
  int xdir = 1;
  if (x0 > x1) {
    xdir = -1;
  }

  int ydir = 1;
  if (y0 > y1) {
    ydir = -1;
  }

  int   x     = x0;
  int   wy    = y0;

  int   dy    = y1-y0;
  //dy    = dy * ydir;
  int   dx    = x1-x0;
  dx    = dx * xdir;
  int   fy    = dy/2;
  
  while (x != x1)
  {
    fb->plotPixel(x, (int)wy, 1.0f, 1.0f, 1.0f);
    x += xdir;

    fy += dy;

    if (fy >= dx) {
      wy += 1;
      fy -= dx;
    }

    
  }

}

int draw_y_line(FrameBuffer *fb, int x0, int y0, int x1, int y1)
{
  int ydir = 1;
  if (y0 > y1) {
    ydir = -1;
  }

  int xdir = 1;
  if (x0 > x1) {
    xdir = -1;
  }

  int   y     = y0;
  int   wx    = x0;

  int   dy    = y1-y0;
  dy    = dy * ydir;
  int   dx    = x1-x0;
  dx    = dx * xdir;
  int   fx    = dx/2;
  
  while (y != y1)
  {
    fb->plotPixel((int)wx, y, 1.0f, 1.0f, 1.0f);
    y += ydir;

    fx += dx;

    if (fx >= dy) {
      wx += 1;
      fx -= dy;
    }
  
  }

}


int draw_line(FrameBuffer *fb, int x0, int y0, int ex, int ey){
  if ((x0 == ex) && (y0==ey)) {
    // dot
    return fb->plotPixel(x0, y0, 1.0f, 1.0f, 1.0f);
    
  } else if (((ex-x0)* (ex-x0)) >= ((ey-y0)* (ey-y0))) { 
    // if  x gradient is greater than y gradient, then it's an x line
    return draw_x_line(fb, x0, y0, ex, ey);
    
  } else {
    return draw_y_line(fb, x0, y0, ex, ey);
  }
}
