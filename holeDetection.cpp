/*
HoleDetection, a program for detecting holes and computing their boundaries from the outer boundary detected Delaunay triangulation of a planar point set.
Copyright (C) 2017 Subhasree Methirumangalath, Shyam Sundar Kannan, Amal Dev Parakkat, and Ramanathan Muthuganapathy.
Advanced Geometric Computing Lab, Department of Engineering Design, IIT Madras, Chennai, India.
For comments or questions, please contact Subhasree Methirumangalath at subhasree.rajiv@gmail.com or Ramanathan Muthuganapathy at mraman@iitm.ac.in.

HoleDetection implements the algorithm described in the paper, "Hole Detection of a Planar Point Set: An Empty Disk Approach" by Subhasree Methirumangalath, Shyam Sundar Kannan, Amal Dev Parakkat, and Ramanathan Muthuganapathy, to appear at the Shape Modeling International (SMI-2017).
*/

/*C++ Header Files*/
#include <fstream>
#include <queue>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <limits>
#include <cstdlib>

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

/*CGAL Header Files*/
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangle_2.h>

/*CGAL type definitions*/
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3   Point;
typedef CGAL::Triangle_3<K> Triangle;
typedef std::list<Triangle> TriangleList;
typedef std::vector<Point> Points;

/*Global Variables*/
Delaunay dt;
float minX = std::numeric_limits<int>::max(), maxX = std::numeric_limits<int>::min();
float minY = std::numeric_limits<int>::max(), maxY = std::numeric_limits<int>::min();
Points inputPoints, exists;
int numberOfHoles = 0;

float v = 1; //parameter for outer bdry detection

int window;
float diagonalDistance;
GLdouble width, height;

class Edge
{
public:
  Point source;
  Point target;
};
typedef std::vector<Edge> Edges;
Edges shape, holeEdges;

/*Priority Queue*/
struct node
{
  float length;
  Delaunay::All_faces_iterator ffi;
  struct node *link;
};

class PriorityQueue
{
private:
  node *front;
public:
  PriorityQueue()
  {
    front = NULL;
  }

  void pqInsert(float length, Delaunay::All_faces_iterator ffi)
  {
    node *tmp, *q;
    tmp = new node;
    tmp->ffi = ffi;
    tmp->length = length;
    if (front == NULL || length > front->length)
    {
      tmp->link = front;
      front = tmp;
    }
    else
    {
      q = front;
      while (q->link != NULL && q->link->length >= length)
      {
        q = q->link;
      }
      tmp->link = q->link;
      q->link = tmp;
    }
  }

  Delaunay::All_faces_iterator pqDelete()
  {
    node *tmp;
    if(front == NULL)
    {
      std::cout << "Queue Underflow" << std::endl;
    }
    else
    {
      tmp = front;
      Delaunay::All_faces_iterator frontFace;
      frontFace = tmp->ffi;
      front = front->link;
      free(tmp);
      return frontFace;
    }
  }

  bool isEmpty()
  {
    if(front == NULL)
    return 1;
    return 0;
  }
};

PriorityQueue pqOuterBoundary; //priority queue for outer boundary detection
PriorityQueue pqInnerBoundary; //priority queue for inner boundary detection
/*End of Priority Queue*/

/*Basic Functions*/
float area(Point a, Point b, Point c) /*compute the area of the triangle formed by the points a, b, and c*/
{
  return (float)std::abs(0.5*(a.x()*(b.y()-c.y())+b.x()*(c.y()-a.y())+c.x()*(a.y()-b.y())));
}

float distance(Point a, Point b) /*compute the distance between two points*/
{
  return (float)(sqrt(std::abs(((a.x()-b.x())*(a.x()-b.x())))+std::abs(((a.y()-b.y())*(a.y()-b.y())))));
}

bool isInfinite(Delaunay::All_faces_iterator afi) /*Checking whether the face is infinite or not*/
{
  if(afi->vertex(0) == dt.infinite_vertex() || afi->vertex(1) == dt.infinite_vertex() || afi->vertex(2) == dt.infinite_vertex())
  return true;
  return false;
}

bool isFinite(Delaunay::All_faces_iterator afi) /*Checking whether the face is finite or not*/
{
  if(afi->vertex(0) != dt.infinite_vertex() && afi->vertex(1) != dt.infinite_vertex() && afi->vertex(2) != dt.infinite_vertex())
  return true;
  return false;
}

bool liesInsideScaled(Point a, Point b, Point c) /*To check whether p3 is inside scaled diametric circle(p1,p2)*/
{
  if(distance( Point((a.x() + b.x())/2, (a.y() + b.y())/2, 1), c) <= v * distance( Point((a.x() + b.x())/2, (a.y() + b.y())/2, 1), a))
  return true;
  return false;
}

bool liesInsideUnscaled(Point a, Point b, Point c) /*To check whether p3 is inside unscaled diametric circle(p1,p2)*/
{
  if(distance( Point((a.x() + b.x())/2, (a.y() + b.y())/2, 1), c) <= distance( Point((a.x() + b.x())/2, (a.y() + b.y())/2, 1), a))
  return true;
  return false;
}
/*End of Basic Functions*/

/*Supporting functions for outer boundary detection*/
bool liesWithinException(Point a, Point b, Point c, float radius)/*Checking whether mid point and chord circles are empty*/
{
  if(distance(a, b) >= radius)
  {
    Point p1 = Point(( (a.x()+b.x())/2 - (radius/2)), ((a.y()+b.y())/2), 1);
    Point p2 = Point(( (a.x()+b.x())/2 + (radius/2)), ((a.y()+b.y())/2), 1); /*(p1,p2) is the midpoint of edge (a,b)*/
    if(liesInsideScaled(p1, p2, c))
    return true;
    return false;
  }
  radius = (radius/2) * v;
  double d = sqrt(std::abs((radius*radius) - ((distance(a,b)/2)*(distance(a,b)/2))));
  float a1New = ((a.x()+b.x())/2)+(d/(sqrt(((a.y()-b.y())*(a.y()-b.y()))+((b.x()-a.x())*(b.x()-a.x())))))*(a.y()-b.y());
  float b1New = ((a.y()+b.y())/2)+(d/(sqrt(((a.y()-b.y())*(a.y()-b.y()))+((b.x()-a.x())*(b.x()-a.x())))))*(b.x()-a.x());
  float a2New = ((a.x()+b.x())/2)-(d/(sqrt(((a.y()-b.y())*(a.y()-b.y()))+((b.x()-a.x())*(b.x()-a.x())))))*(a.y()-b.y());
  float b2New = ((a.y()+b.y())/2)-(d/(sqrt(((a.y()-b.y())*(a.y()-b.y()))+((b.x()-a.x())*(b.x()-a.x())))))*(b.x()-a.x());

  if(liesInsideUnscaled( Point(a1New+radius, b1New, 1),Point(a1New-radius, b1New, 1), c))
  return true;
  if(liesInsideUnscaled( Point(a2New+radius, b2New,1),Point(a2New-radius, b2New, 1), c))
  return true;
  return false;
}

bool isInShape(Point a)/*Checking whether a point is already in shape*/
{
  for(int i = 0; i < exists.size(); i++)
  {
    if(exists[i] == a)
    return true;
  }
  return false;
}

void updateNeighbors(Delaunay::All_faces_iterator currentFace, Delaunay::All_faces_iterator neighbor1, Delaunay::All_faces_iterator neighbor2, float d1, float d2)
{
  if(neighbor1->neighbor(0) == currentFace)
  {
    neighbor1->set_neighbor(0,dt.infinite_face());
  }
  else
  {
    if(neighbor1->neighbor(1) == currentFace)
    {
      neighbor1->set_neighbor(1,dt.infinite_face());
    }
    else
    {
      if(neighbor1->neighbor(2) == currentFace)
      {
        neighbor1->set_neighbor(2,dt.infinite_face());
      }
    }
  }

  if(neighbor2->neighbor(0) == currentFace)
  {
    neighbor2->set_neighbor(0,dt.infinite_face());
  }
  else
  {
    if(neighbor2->neighbor(1) == currentFace)
    {
      neighbor2->set_neighbor(1,dt.infinite_face());
    }
    else
    {
      if(neighbor2->neighbor(2) == currentFace)
      {
        neighbor2->set_neighbor(2,dt.infinite_face());
      }
    }
  }


  dt.delete_face(currentFace);
  pqOuterBoundary.pqInsert(d1, neighbor1);
  pqOuterBoundary.pqInsert(d2, neighbor2);
}

void insertToShape(Point a, Point b)/*Inserting a new edge into ec-shape*/
{
  Edge e;
  e.source = a;
  e.target = b;
  shape.push_back(e);
}
/*End of Supporting functions for outer boundary detection*/


/*Supporting functions for hole detection*/
Delaunay::Face_handle getMaxAreaTriangleHandle(Delaunay::Finite_faces_iterator ffi) /*returns the largest triangle*/
{
  double biggestArea = 0.0;
  bool presentInShape = false;
  Delaunay::Face_handle highestAreaTriangleHandle;
  ffi = dt.finite_faces_begin();
  int fl = 0;

  for(Delaunay::Finite_faces_iterator ffi = dt.finite_faces_begin(); ffi != dt.finite_faces_end(); ffi++)
  {
    Point a = ffi->vertex(0)->point();
    Point b = ffi->vertex(1)->point();
    Point c = ffi->vertex(2)->point();

    for (int i = 0; i < shape.size(); i++)
    {
      if ((a.x() == shape[i].source.x() && a.y()==shape[i].source.y()) ||(a.x()==shape[i].target.x() && a.y()==shape[i].target.y())||
      (b.x() == shape[i].source.x() && b.y()==shape[i].source.y()) || (b.x()==shape[i].target.x() && b.y()==shape[i].target.y())||
      (c.x() == shape[i].source.x() && c.y()==shape[i].source.y()) || (c.x()==shape[i].target.x() && c.y()==shape[i].target.y()))
      presentInShape = true;
    }

    if (!presentInShape)
    {
      if (isFinite(ffi->neighbor(0)) && isFinite(ffi->neighbor(1)) && isFinite(ffi->neighbor(2)))
      {
        fl = 1;
        if(area(ffi->vertex(0)->point(), ffi->vertex(1)->point(), ffi->vertex(2)->point()) >= biggestArea)
        {
          biggestArea = area(ffi->vertex(0)->point(), ffi->vertex(1)->point(), ffi->vertex(2)->point());
          highestAreaTriangleHandle = ffi;
        }
      }
    }
    presentInShape = false;
  }
  if(fl == 0)
  {
    highestAreaTriangleHandle = dt.infinite_face();
  }
  return highestAreaTriangleHandle;
}

bool liesInsideHole(Point p1, Point p2, Point p3)
{
  if(distance(Point( (p1.x()+p2.x())/2, (p1.y()+p2.y())/2, 1), p3) <= distance(Point( (p1.x()+p2.x())/2, (p1.y()+p2.y())/2, 1), p1))
  {
    return true;
  }
  return 0;
}

bool liesWithinExceptionHole(Point a, Point b, Point c, float radius)/*Checking whether mid point and chord circles are empty*/
{
  if(distance(a, b) >= radius)
  {
    Point p1 = Point(((a.x()+b.x())/2 - (radius/2)), ((a.y()+b.y())/2), 1);
    Point p2 = Point(((a.x()+b.x())/2 + (radius/2)), ((a.y()+b.y())/2), 1);/*(p1,p2) is the midpoint of edge (a,b)*/
    if(liesInsideScaled(p1, p2, c))
    {
      return true;
    }
    return false;
  }
  radius = (radius/2);
  double d = sqrt(std::abs((radius*radius) - ((distance(a, b)/2)*(distance(a, b)/2))));
  float a1New = ((a.x()+b.x())/2)+(d/(sqrt(((a.y()-b.y())*(a.y()-b.y()))+((b.x()-a.x())*(b.x()-a.x())))))*(a.y()-b.y());
  float b1New = ((a.y()+b.y())/2)+(d/(sqrt(((a.y()-b.y())*(a.y()-b.y()))+((b.x()-a.x())*(b.x()-a.x())))))*(b.x()-a.x());
  float a2New = ((a.x()+b.x())/2)-(d/(sqrt(((a.y()-b.y())*(a.y()-b.y()))+((b.x()-a.x())*(b.x()-a.x())))))*(a.y()-b.y());
  float b2New = ((a.y()+b.y())/2)-(d/(sqrt(((a.y()-b.y())*(a.y()-b.y()))+((b.x()-a.x())*(b.x()-a.x())))))*(b.x()-a.x());
  if(liesInsideUnscaled(Point(a1New+radius, b1New, 1), Point(a1New-radius, b1New,1 ), c))
  {
    return true;
  }
  if(liesInsideUnscaled(Point(a2New+radius, b2New, 1), Point(a2New-radius, b2New, 1), c))
  {
    return true;
  }
  return false;
}

void insertToHoleBoundary( Point a, Point b)/*Inserting a new edge into ec-shape*/
{
  Edge e;
  e.source = a;
  e.target = b;
  holeEdges.push_back(e);
}

void updateHole(Delaunay::All_faces_iterator currentFace, Delaunay::All_faces_iterator neighbor1, Delaunay::All_faces_iterator neighbor2) /*updates the hole boundary*/
{
  if(neighbor1->neighbor(0) == currentFace)
  {
    neighbor1->set_neighbor(0, dt.infinite_face());
  }
  else
  {
    if(neighbor1->neighbor(1) == currentFace)
    {
      neighbor1->set_neighbor(1, dt.infinite_face());
    }
    else
    {
      if(neighbor1->neighbor(2) == currentFace)
      {
        neighbor1->set_neighbor(2, dt.infinite_face());
      }
    }
  }

  if(neighbor2->neighbor(0) == currentFace)
  {
    neighbor2->set_neighbor(0, dt.infinite_face());
  }
  else
  {
    if(neighbor2->neighbor(1) == currentFace)
    {
      neighbor2->set_neighbor(1, dt.infinite_face());
    }
    else
    {
      if(neighbor2->neighbor(2) == currentFace)
      {
        neighbor2->set_neighbor(2, dt.infinite_face());
      }
    }
  }

  dt.delete_face(currentFace);

  pqInnerBoundary.pqInsert(area(neighbor1->vertex(0)->point(), neighbor1->vertex(1)->point(), neighbor1->vertex(2)->point()), neighbor1);
  pqInnerBoundary.pqInsert(area(neighbor2->vertex(0)->point(), neighbor2->vertex(1)->point(), neighbor2->vertex(2)->point()), neighbor2);
}

bool inBoundary(Point a)/*Checking whether a point is already in shape*/
{
  for(int i = 0;  i < exists.size(); i++)
  {
    if(exists[i] == a)
    {
      return true;
    }
  }
  return false;
}
/*End of functions for hole detection*/

/*Outer Boundary detection condition*/
bool outerBoundaryCondition(Point a, Point b, Point c, Delaunay::Face_handle fh)
{
  if(liesInsideScaled(a, b, c)) /*Diametric circle is non-empty*/
  {
    return true;
  }

  /*Checking the status of mid point circle and chord circles*/
  if(isFinite(fh->neighbor(0)))
  {
    if(fh->neighbor(0)->vertex(0)->point() !=a && fh->neighbor(0)->vertex(0)->point() != b && fh->neighbor(0)->vertex(0)->point() !=c)
    {
      if(liesWithinException(fh->neighbor(0)->vertex(1)->point(), fh->neighbor(0)->vertex(2)->point(), fh->neighbor(0)->vertex(0)->point(), v*distance(a,b)))
      return true;
    }
    else
    {
      if(fh->neighbor(0)->vertex(1)->point() !=a && fh->neighbor(0)->vertex(1)->point() != b && fh->neighbor(0)->vertex(1)->point() != c)
      {
        if(liesWithinException(fh->neighbor(0)->vertex(0)->point(), fh->neighbor(0)->vertex(2)->point(), fh->neighbor(0)->vertex(1)->point(), v*distance(a,b)))
        return true;
      }
      else
      {
        if(fh->neighbor(0)->vertex(2)->point() !=a && fh->neighbor(0)->vertex(2)->point() !=b && fh->neighbor(0)->vertex(2)->point() != c)
        {
          if(liesWithinException(fh->neighbor(0)->vertex(1)->point(), fh->neighbor(0)->vertex(0)->point(), fh->neighbor(0)->vertex(2)->point(), v*distance(a,b)))
          return true;
        }
      }
    }
  }

  if(isFinite(fh->neighbor(1)))
  {
    if(fh->neighbor(1)->vertex(0)->point() !=a && fh->neighbor(1)->vertex(0)->point() != b && fh->neighbor(1)->vertex(0)->point() !=c)
    {
      if(liesWithinException(fh->neighbor(1)->vertex(1)->point(), fh->neighbor(1)->vertex(2)->point(), fh->neighbor(1)->vertex(0)->point(), v*distance(a,b)))
      return true;
    }
    else
    {
      if(fh->neighbor(1)->vertex(1)->point() != a && fh-> neighbor(1)->vertex(1)->point() !=b && fh->neighbor(1)->vertex(1)->point() != c)
      {
        if(liesWithinException(fh->neighbor(1)->vertex(2)->point(), fh->neighbor(1)->vertex(0)->point(), fh->neighbor(1)->vertex(1)->point(), v*distance(a,b)))
        return true;
      }
      else
      {
        if(fh->neighbor(1)->vertex(2)->point() != a && fh->neighbor(1)->vertex(2)->point() !=b && fh->neighbor(1)->vertex(2)->point() != c)
        {
          if(liesWithinException(fh->neighbor(1)->vertex(0)->point(), fh->neighbor(1)->vertex(1)->point(), fh->neighbor(1)->vertex(2)->point(), v*distance(a,b)))
          return true;
        }
      }
    }
  }

  if(isFinite(fh->neighbor(2)))
  {
    if(fh->neighbor(2)->vertex(0)->point() != a && fh->neighbor(2)->vertex(0)->point() != b && fh->neighbor(2)->vertex(0)->point() != c)
    {
      if(liesWithinException(fh->neighbor(2)->vertex(1)->point(), fh->neighbor(2)->vertex(2)->point(), fh->neighbor(2)->vertex(0)->point(), v*distance(a,b)))
      return true;
    }
    else
    {
      if(fh->neighbor(2)->vertex(1)->point() != a && fh->neighbor(2)->vertex(1)->point() != b && fh->neighbor(2)->vertex(1)->point() != c)
      {
        if(liesWithinException(fh->neighbor(2)->vertex(2)->point(), fh->neighbor(2)->vertex(0)->point(), fh->neighbor(2)->vertex(1)->point(), v*distance(a,b)))
        return true;
      }
      else
      {
        if(fh->neighbor(2)->vertex(2)->point() != a && fh->neighbor(2)->vertex(2)->point() != b && fh->neighbor(2)->vertex(2)->point() != c)
        {
          if(liesWithinException(fh->neighbor(2)->vertex(0)->point(),fh->neighbor(2)->vertex(1)->point(),fh->neighbor(2)->vertex(2)->point(),v*distance(a,b)))
          return true;
        }
      }
    }
  }

  return false;
}
/*End of Outer Boundary detection condition*/

/*Hole Detection condition*/
bool holeDetectionCondition(Point a, Point b, Point c, Delaunay::Face_handle fh)
{
  if(liesInsideHole(a, b, c))/*Diametric circle is non-empty*/
  {
    return true;
  }

  /*Checking the status of mid point circle and chord circles*/
  if(isFinite(fh->neighbor(0)))
  {
    if(fh->neighbor(0)->vertex(0)->point() !=a && fh->neighbor(0)->vertex(0)->point() !=b && fh->neighbor(0)->vertex(0)->point() != c)
    {
      if(liesWithinExceptionHole(fh->neighbor(0)->vertex(1)->point(), fh->neighbor(0)->vertex(2)->point(), fh->neighbor(0)->vertex(0)->point(), distance(a, b)))
      {
        return true;
      }
    }
    else
    {
      if(fh->neighbor(0)->vertex(1)->point() != a && fh->neighbor(0)->vertex(1)->point() != b && fh->neighbor(0)->vertex(1)->point() != c)
      {
        if(liesWithinExceptionHole(fh->neighbor(0)->vertex(0)->point(), fh->neighbor(0)->vertex(2)->point(), fh->neighbor(0)->vertex(1)->point(), distance(a, b)))
        {
          return true;
        }
      }
      else
      {
        if(fh->neighbor(0)->vertex(2)->point() != a && fh->neighbor(0)->vertex(2)->point() != b && fh->neighbor(0)->vertex(2)->point() != c)
        {
          if(liesWithinExceptionHole(fh->neighbor(0)->vertex(1)->point(), fh->neighbor(0)->vertex(0)->point(), fh->neighbor(0)->vertex(2)->point(), distance(a,b)))
          {
            return true;
          }
        }
      }
    }
  }

  if(isFinite(fh->neighbor(1)))
  {
    if(fh->neighbor(1)->vertex(0)->point() != a && fh->neighbor(1)->vertex(0)->point() != b && fh->neighbor(1)->vertex(0)->point() != c)
    {
      if(liesWithinExceptionHole(fh->neighbor(1)->vertex(1)->point(), fh->neighbor(1)->vertex(2)->point(), fh->neighbor(1)->vertex(0)->point(), distance(a, b)))
      {
        return true;
      }
    }
    else
    {
      if(fh->neighbor(1)->vertex(1)->point() != a && fh->neighbor(1)->vertex(1)->point() != b && fh->neighbor(1)->vertex(1)->point() != c)
      {
        if(liesWithinExceptionHole(fh->neighbor(1)->vertex(2)->point(), fh->neighbor(1)->vertex(0)->point(), fh->neighbor(1)->vertex(1)->point(), distance(a, b)))
        {
          return true;
        }
      }
      else
      {
        if(fh->neighbor(1)->vertex(2)->point() != a && fh->neighbor(1)->vertex(2)->point() != b && fh->neighbor(1)->vertex(2)->point() != c)
        {
          if(liesWithinExceptionHole(fh->neighbor(1)->vertex(0)->point(), fh->neighbor(1)->vertex(1)->point(), fh->neighbor(1)->vertex(2)->point(), distance(a, b)))
          {
            return true;
          }
        }
      }
    }
  }

  if(isFinite(fh->neighbor(2)))
  {
    if(fh->neighbor(2)->vertex(0)->point() != a && fh->neighbor(2)->vertex(0)->point() != b && fh->neighbor(2)->vertex(0)->point() != c)
    {
      if(liesWithinExceptionHole(fh->neighbor(2)->vertex(1)->point(), fh->neighbor(2)->vertex(2)->point(), fh->neighbor(2)->vertex(0)->point(), distance(a,b)))
      {
        return true;
      }
    }
    else
    {
      if(fh->neighbor(2)->vertex(1)->point() !=a && fh->neighbor(2)->vertex(1)->point() != b && fh->neighbor(2)->vertex(1)->point() != c)
      {
        if(liesWithinExceptionHole(fh->neighbor(2)->vertex(2)->point(), fh->neighbor(2)->vertex(0)->point(), fh->neighbor(2)->vertex(1)->point(), distance(a, b)))
        {
          return true;
        }
      }
      else
      {
        if(fh->neighbor(2)->vertex(2)->point() != a && fh->neighbor(2)->vertex(2)->point() != b && fh->neighbor(2)->vertex(2)->point() != c)
        {
          if(liesWithinExceptionHole(fh->neighbor(2)->vertex(0)->point(), fh->neighbor(2)->vertex(1)->point(), fh->neighbor(2)->vertex(2)->point(), distance(a,b)))
          {
            return true;
          }
        }
      }
    }
  }

  return false;
}
/*End of Hole Detection condition*/

/*Function to detect the outer boundary*/
void outerBoundary()
{
  for(Delaunay::Finite_faces_iterator ffi = dt.finite_faces_begin(); ffi != dt.finite_faces_end(); ffi++)
  {
    if(isInfinite(ffi->neighbor(0)) || isInfinite(ffi->neighbor(1)) || isInfinite(ffi->neighbor(2)))
    {
      if(isInfinite(ffi->neighbor(0)))
      {
        ffi->neighbor(0) = dt.infinite_face(); //set the 0th neighbor as infinite face
        exists.push_back(ffi->vertex(1)->point());
        exists.push_back(ffi->vertex(2)->point());
        pqOuterBoundary.pqInsert(distance(ffi->vertex(1)->point(), ffi->vertex(2)->point()), ffi);
      }
      else
      {
        if(isInfinite(ffi->neighbor(1)))
        {
          ffi->neighbor(1) = dt.infinite_face();
          exists.push_back(ffi->vertex(0)->point());
          exists.push_back(ffi->vertex(2)->point());
          pqOuterBoundary.pqInsert(distance(ffi->vertex(0)->point(), ffi->vertex(2)->point()), ffi);
        }
        else
        {
          if(isInfinite(ffi->neighbor(2)))
          {
            ffi->neighbor(2) = dt.infinite_face();
            exists.push_back(ffi->vertex(1)->point());
            exists.push_back(ffi->vertex(0)->point());
            pqOuterBoundary.pqInsert(distance(ffi->vertex(0)->point(), ffi->vertex(1)->point()), ffi);
          }
        }
      }
    }
  }

  while(!pqOuterBoundary.isEmpty())
  {
    Delaunay::All_faces_iterator currentFace = pqOuterBoundary.pqDelete();

    if(isInfinite(currentFace->neighbor(0)) && isFinite(currentFace->neighbor(1)) && isFinite(currentFace->neighbor(2)))
    {
      if(outerBoundaryCondition(currentFace->vertex(1)->point(), currentFace->vertex(2)->point(), currentFace->vertex(0)->point(), currentFace))/*If any of the circle is non-empty*/
      {
        if(!isInShape(currentFace->vertex(0)->point()))
        {
          exists.push_back(currentFace->vertex(0)->point());
          updateNeighbors(currentFace, currentFace->neighbor(1), currentFace->neighbor(2), distance(currentFace->vertex(0)->point(), currentFace->vertex(2)->point()), distance(currentFace->vertex(1)->point(), currentFace->vertex(0)->point()));
        }
        else
        {
          insertToShape(currentFace->vertex(1)->point(), currentFace->vertex(2)->point());
        }
      }
      else
      {
        insertToShape(currentFace->vertex(1)->point(), currentFace->vertex(2)->point());
      }
    }

    else
    {
      if(isInfinite(currentFace->neighbor(1)) && isFinite(currentFace->neighbor(0)) && isFinite(currentFace->neighbor(2)))
      {
        if(outerBoundaryCondition(currentFace->vertex(0)->point(), currentFace->vertex(2)->point(), currentFace->vertex(1)->point(), currentFace))
        {
          if(!isInShape(currentFace->vertex(1)->point()))
          {
            exists.push_back(currentFace->vertex(1)->point());
            updateNeighbors(currentFace, currentFace->neighbor(2), currentFace->neighbor(0), distance(currentFace->vertex(0)->point(), currentFace->vertex(1)->point()), distance(currentFace->vertex(1)->point(), currentFace->vertex(2)->point()));
          }
          else
          {
            insertToShape(currentFace->vertex(0)->point(), currentFace->vertex(2)->point());
          }
        }
        else
        {
          insertToShape(currentFace->vertex(0)->point(), currentFace->vertex(2)->point());
        }
      }
      else
      {
        if(isInfinite(currentFace->neighbor(2)) && isFinite(currentFace->neighbor(1)) && isFinite(currentFace->neighbor(0)))
        {
          if(outerBoundaryCondition(currentFace->vertex(1)->point(), currentFace->vertex(0)->point(), currentFace->vertex(2)->point(), currentFace))
          {
            if(!isInShape(currentFace->vertex(2)->point()))
            {
              exists.push_back(currentFace->vertex(2)->point());
              updateNeighbors(currentFace, currentFace->neighbor(1), currentFace->neighbor(0), distance(currentFace->vertex(0)->point(), currentFace->vertex(2)->point()), distance(currentFace->vertex(1)->point(), currentFace->vertex(2)->point()));
            }
            else
            {
              insertToShape(currentFace->vertex(0)->point(), currentFace->vertex(1)->point());
            }
          }
          else
          {
            insertToShape(currentFace->vertex(0)->point(), currentFace->vertex(1)->point());
          }
        }
        else
        {
          if(isInfinite(currentFace->neighbor(2)) && isInfinite(currentFace->neighbor(1)) && isFinite(currentFace->neighbor(0)))
          {
            insertToShape(currentFace->vertex(2)->point(), currentFace->vertex(0)->point());
            insertToShape(currentFace->vertex(0)->point(), currentFace->vertex(1)->point());
          }
          else
          {
            if(isInfinite(currentFace->neighbor(2)) && isFinite(currentFace->neighbor(1)) && isInfinite(currentFace->neighbor(0)))
            {
              insertToShape(currentFace->vertex(2)->point(), currentFace->vertex(1)->point());
              insertToShape(currentFace->vertex(1)->point(), currentFace->vertex(0)->point());
            }
            else
            {
              if(isFinite(currentFace->neighbor(2)) && isInfinite(currentFace->neighbor(1)) && isInfinite(currentFace->neighbor(0)))
              {
                insertToShape(currentFace->vertex(0)->point(), currentFace->vertex(2)->point());
                insertToShape(currentFace->vertex(2)->point(), currentFace->vertex(1)->point());
              }
            }
          }
        }
      }
    }
  }
}
/*End of outerBoundary function*/

/*Function to detect the hole boundary*/
void holeDetection()
{
  Delaunay::Finite_faces_iterator ffi = dt.finite_faces_begin();

  Delaunay::Face_handle highestAreaTriangleHandle;//the face handle for the highest area triangle
  Delaunay::Face_handle faceHandle0, faceHandle1, faceHandle2;//facehandles for neighbouring triangles
  bool inHole;
  int tempVar = 0;//variable for for multiple hole detection
  bool zerothPointFaceHandle0, firstPointFaceHandle0, secondPointFaceHandle0;
  bool zerothPointFaceHandle1, firstPointFaceHandle1, secondPointFaceHandle1;
  bool zerothPointFaceHandle2, firstPointFaceHandle2, secondPointFaceHandle2;
  bool inOuterBdry = false;

  //do while loop for multiple hole detection
  do //do of while ((count1-tempVar)>4);
  {
    tempVar = holeEdges.size(); //variable for multiple hole detection
    inHole = false;//variable to check whether already the edge is in OHL
    zerothPointFaceHandle0 = false, firstPointFaceHandle0 = false, secondPointFaceHandle0 = false;
    zerothPointFaceHandle1 = false, firstPointFaceHandle1 = false, secondPointFaceHandle1 = false;
    zerothPointFaceHandle2 = false, firstPointFaceHandle2 = false, secondPointFaceHandle2 = false;
    inOuterBdry = false;

    highestAreaTriangleHandle = getMaxAreaTriangleHandle(ffi);
    if(highestAreaTriangleHandle == dt.infinite_face())
    {
      if (holeEdges.size() == 0)
      {
        exit(0);
      }
      else
      {
        break;
      }
    }

    //checking any of the vertices of HAT are in the outerbdry, if so do not put that it to Q
    for(int i =0; i < holeEdges.size(); i++)
    {
      if ((highestAreaTriangleHandle->vertex(0)->point().x()==holeEdges[i].source.x()&& highestAreaTriangleHandle->vertex(0)->point().y()==holeEdges[i].source.y())||
      (highestAreaTriangleHandle->vertex(0)->point().x()==holeEdges[i].target.x()&& highestAreaTriangleHandle->vertex(0)->point().y()==holeEdges[i].target.y())||
      (highestAreaTriangleHandle->vertex(1)->point().x()==holeEdges[i].source.x()&&highestAreaTriangleHandle->vertex(1)->point().y()==holeEdges[i].source.y())||
      (highestAreaTriangleHandle->vertex(1)->point().x()==holeEdges[i].target.x()&&highestAreaTriangleHandle->vertex(1)->point().y()==holeEdges[i].target.y())||
      (highestAreaTriangleHandle->vertex(2)->point().x()==holeEdges[i].source.x()&&highestAreaTriangleHandle->vertex(2)->point().y()==holeEdges[i].source.y())||
      (highestAreaTriangleHandle->vertex(2)->point().x()==holeEdges[i].target.x()&&highestAreaTriangleHandle->vertex(2)->point().y()==holeEdges[i].target.y()))
      {
        inHole = true;
        break;
      }
    }

    //checking any of the vertices of HAT are in the outerbdry, if so do not put that it to Q
    for(int i = 0; i < shape.size(); i++)
    {
      if ((highestAreaTriangleHandle->vertex(0)->point().x()==shape[i].source.x()&& highestAreaTriangleHandle->vertex(0)->point().y()==shape[i].source.y())||
      (highestAreaTriangleHandle->vertex(0)->point().x()==shape[i].target.x()&& highestAreaTriangleHandle->vertex(0)->point().y()==shape[i].target.y())||
      (highestAreaTriangleHandle->vertex(1)->point().x()==shape[i].source.x()&& highestAreaTriangleHandle->vertex(1)->point().y()==shape[i].source.y())||
      (highestAreaTriangleHandle->vertex(1)->point().x()==shape[i].target.x()&& highestAreaTriangleHandle->vertex(1)->point().y()==shape[i].target.y())||
      (highestAreaTriangleHandle->vertex(2)->point().x()==shape[i].source.x()&& highestAreaTriangleHandle->vertex(2)->point().y()==shape[i].source.y())||
      (highestAreaTriangleHandle->vertex(2)->point().x()==shape[i].target.x()&& highestAreaTriangleHandle->vertex(2)->point().y()==shape[i].target.y()))
      {
        inOuterBdry = true;
        break;
      }
    }

    if (inOuterBdry)
    {
      std::cout << "highest area triangle is outer bdry triangle and it is an invalid input" << std::endl;
      if (holeEdges.size() == 0)
      {
        exit(0);
      }
      else
      {
        break;
      }
    }
    else if ((!inOuterBdry) && (!inHole))//if the hole is not along the outer bdry, initialize the queue1 with the highest area traingle
    {
      faceHandle0 = highestAreaTriangleHandle->neighbor(0);
      pqInnerBoundary.pqInsert(area(faceHandle0->vertex(0)->point(), faceHandle0->vertex(1)->point(), faceHandle0->vertex(2)->point()), faceHandle0);

      faceHandle1 = highestAreaTriangleHandle->neighbor(1);
      pqInnerBoundary.pqInsert(area(faceHandle1->vertex(0)->point(), faceHandle1->vertex(1)->point(), faceHandle1->vertex(2)->point()), faceHandle1);

      faceHandle2 = highestAreaTriangleHandle->neighbor(2);
      pqInnerBoundary.pqInsert(area(faceHandle2->vertex(0)->point(), faceHandle2->vertex(1)->point(), faceHandle2->vertex(2)->point()), faceHandle2);

      zerothPointFaceHandle2 = firstPointFaceHandle2 = secondPointFaceHandle2 = false;
      //the next three if loops are to find out the vertex number of the third vertex w.r.t facehandle2 when v0,v1 is processed. this is to use it in lies inside condition checking for diacircle
      //to get the unvisited vertex of the neighboring triangle of 2nd vertex for cheking lies inside condition

      if (((faceHandle2->vertex(0)->point().x()==highestAreaTriangleHandle->vertex(0)->point().x()) && (faceHandle2->vertex(0)->point().y()==highestAreaTriangleHandle->vertex(0)->point().y())) ||
      ((faceHandle2->vertex(0)->point().x()==highestAreaTriangleHandle->vertex(1)->point().x()) && (faceHandle2->vertex(0)->point().y()==highestAreaTriangleHandle->vertex(1)->point().y())) )
      {
        zerothPointFaceHandle2 = false;
      }
      else
      {
        zerothPointFaceHandle2 = 1;
      }

      if	((faceHandle2->vertex(1)->point().x()==highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle2->vertex(1)->point().y()==highestAreaTriangleHandle->vertex(0)->point().y()) ||
      (faceHandle2->vertex(1)->point().x()==highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle2->vertex(1)->point().y()==highestAreaTriangleHandle->vertex(1)->point().y()) )
      {
        firstPointFaceHandle2 = false;
      }
      else
      {
        firstPointFaceHandle2 = true;
      }

      if	((faceHandle2->vertex(2)->point().x()==highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle2->vertex(2)->point().y()==highestAreaTriangleHandle->vertex(0)->point().y()) ||
      (faceHandle2->vertex(2)->point().x()==highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle2->vertex(2)->point().y()==highestAreaTriangleHandle->vertex(1)->point().y()) )
      {
        secondPointFaceHandle2 = false;
      }
      else
      {
        secondPointFaceHandle2 = true;
      }

      if (zerothPointFaceHandle2)
      {
        faceHandle2->set_neighbor(0, dt.infinite_face());
      }
      else if (firstPointFaceHandle2)
      {
        faceHandle2->set_neighbor(1, dt.infinite_face());
      }


      else if (secondPointFaceHandle2)
      {
        faceHandle2->set_neighbor(2, dt.infinite_face());
      }

      zerothPointFaceHandle0 = firstPointFaceHandle0 = secondPointFaceHandle0 = false;
      if (((faceHandle0->vertex(0)->point().x()==highestAreaTriangleHandle->vertex(1)->point().x()) && (faceHandle0->vertex(0)->point().y()==highestAreaTriangleHandle->vertex(1)->point().y())) ||
      ((faceHandle0->vertex(0)->point().x()==highestAreaTriangleHandle->vertex(2)->point().x()) && (faceHandle0->vertex(0)->point().y()==highestAreaTriangleHandle->vertex(2)->point().y())) )
      {
        zerothPointFaceHandle0 = false;
      }
      else
      {
        zerothPointFaceHandle0 = true;
      }

      if	((faceHandle0->vertex(1)->point().x()==highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle0->vertex(1)->point().y()==highestAreaTriangleHandle->vertex(1)->point().y()) ||
      (faceHandle0->vertex(1)->point().x()==highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle0->vertex(1)->point().y()==highestAreaTriangleHandle->vertex(2)->point().y()) )
      {
        firstPointFaceHandle0 = false;
      }
      else
      {
        firstPointFaceHandle0 = true;
      }

      if	((faceHandle0->vertex(2)->point().x()==highestAreaTriangleHandle->vertex(1)->point().x() && faceHandle0->vertex(2)->point().y()==highestAreaTriangleHandle->vertex(1)->point().y()) ||
      (faceHandle0->vertex(2)->point().x()==highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle0->vertex(2)->point().y()==highestAreaTriangleHandle->vertex(2)->point().y()) )
      {
        secondPointFaceHandle0 = false;
      }
      else
      {
        secondPointFaceHandle0 = true;
      }

      if (zerothPointFaceHandle0)
      {
        faceHandle0->set_neighbor(0, dt.infinite_face());
      }
      else if (firstPointFaceHandle0)
      {
        faceHandle0->set_neighbor(1, dt.infinite_face());
      }
      else if (secondPointFaceHandle0)
      {
        faceHandle0->set_neighbor(2, dt.infinite_face());
      }

      zerothPointFaceHandle1 = firstPointFaceHandle1 = secondPointFaceHandle1 = false;

      if (((faceHandle1->vertex(0)->point().x()==highestAreaTriangleHandle->vertex(2)->point().x()) && (faceHandle1->vertex(0)->point().y()==highestAreaTriangleHandle->vertex(2)->point().y())) ||
      ((faceHandle1->vertex(0)->point().x()==highestAreaTriangleHandle->vertex(0)->point().x()) && (faceHandle1->vertex(0)->point().y()==highestAreaTriangleHandle->vertex(0)->point().y())) )
      {
        zerothPointFaceHandle1 = false;
      }
      else
      {
        zerothPointFaceHandle1 = true;
      }

      if	((faceHandle1->vertex(1)->point().x()==highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle1->vertex(1)->point().y()==highestAreaTriangleHandle->vertex(2)->point().y()) ||
      (faceHandle1->vertex(1)->point().x()==highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle1->vertex(1)->point().y()==highestAreaTriangleHandle->vertex(0)->point().y()) )
      {
        firstPointFaceHandle1 = false;
      }
      else
      {
        firstPointFaceHandle1 = true;
      }

      if	((faceHandle1->vertex(2)->point().x()==highestAreaTriangleHandle->vertex(2)->point().x() && faceHandle1->vertex(2)->point().y()==highestAreaTriangleHandle->vertex(2)->point().y()) ||
      (faceHandle1->vertex(2)->point().x()==highestAreaTriangleHandle->vertex(0)->point().x() && faceHandle1->vertex(2)->point().y()==highestAreaTriangleHandle->vertex(0)->point().y()) )
      {
        secondPointFaceHandle1 = false;
      }
      else
      {
        secondPointFaceHandle1 = true;
      }


      if (zerothPointFaceHandle1)
      {
        faceHandle1->set_neighbor(0, dt.infinite_face());
      }
      else if (firstPointFaceHandle1)
      {
        faceHandle1->set_neighbor(1, dt.infinite_face());
      }
      else if (secondPointFaceHandle1)
      {
        faceHandle1->set_neighbor(2, dt.infinite_face());
      }

      exists.push_back(highestAreaTriangleHandle->vertex(0)->point());
      exists.push_back(highestAreaTriangleHandle->vertex(1)->point());
      exists.push_back(highestAreaTriangleHandle->vertex(2)->point());
      dt.delete_face(highestAreaTriangleHandle);
    }

    while(!pqInnerBoundary.isEmpty())
    {
      highestAreaTriangleHandle = pqInnerBoundary.pqDelete();

      if(isInfinite(highestAreaTriangleHandle->neighbor(0)) && isFinite(highestAreaTriangleHandle->neighbor(1)) && isFinite(highestAreaTriangleHandle->neighbor(2)))
      {
        if(holeDetectionCondition(highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle))
        {
          if(!inBoundary(highestAreaTriangleHandle->vertex(0)->point()))
          {
            exists.push_back(highestAreaTriangleHandle->vertex(0)->point());
            updateHole(highestAreaTriangleHandle,highestAreaTriangleHandle->neighbor(1),highestAreaTriangleHandle->neighbor(2));
          }
          else
          {
            insertToHoleBoundary(highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(2)->point());
          }
        }
        else
        {
          insertToHoleBoundary(highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(2)->point());
        }
      }
      else if(isInfinite(highestAreaTriangleHandle->neighbor(1)) && isFinite(highestAreaTriangleHandle->neighbor(0)) && isFinite(highestAreaTriangleHandle->neighbor(2)))
      {
        if(holeDetectionCondition(highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle))
        {
          if(!inBoundary(highestAreaTriangleHandle->vertex(1)->point()))
          {
            exists.push_back(highestAreaTriangleHandle->vertex(1)->point());
            updateHole(highestAreaTriangleHandle,highestAreaTriangleHandle->neighbor(2), highestAreaTriangleHandle->neighbor(0));
          }
          else
          {
            insertToHoleBoundary(highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point());
          }
        }
        else
        {
          insertToHoleBoundary(highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point());
        }
      }
      else if(isInfinite(highestAreaTriangleHandle->neighbor(2)) && isFinite(highestAreaTriangleHandle->neighbor(1)) && isFinite(highestAreaTriangleHandle->neighbor(0)))
      {
        if(holeDetectionCondition(highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle))
        {
          if(!inBoundary(highestAreaTriangleHandle->vertex(2)->point()))
          {
            exists.push_back(highestAreaTriangleHandle->vertex(2)->point());
            updateHole(highestAreaTriangleHandle, highestAreaTriangleHandle->neighbor(1), highestAreaTriangleHandle->neighbor(0));
          }
          else
          {
            insertToHoleBoundary(highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(1)->point());
          }
        }
        else
        {
          insertToHoleBoundary(highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(1)->point());
        }
      }
      else
      {
        if(isInfinite(highestAreaTriangleHandle->neighbor(2)) && isInfinite(highestAreaTriangleHandle->neighbor(1)) && isFinite(highestAreaTriangleHandle->neighbor(0)))
        {
          insertToHoleBoundary(highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(0)->point());
          insertToHoleBoundary(highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(1)->point());
        }
        else if(isInfinite(highestAreaTriangleHandle->neighbor(2)) && isFinite(highestAreaTriangleHandle->neighbor(1)) && isInfinite(highestAreaTriangleHandle->neighbor(0)))
        {
          insertToHoleBoundary(highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(1)->point());
          insertToHoleBoundary(highestAreaTriangleHandle->vertex(1)->point(), highestAreaTriangleHandle->vertex(0)->point());
        }
        else if(isFinite(highestAreaTriangleHandle->neighbor(2))&&isInfinite(highestAreaTriangleHandle->neighbor(1))&&isInfinite(highestAreaTriangleHandle->neighbor(0)))
        {
          insertToHoleBoundary(highestAreaTriangleHandle->vertex(0)->point(), highestAreaTriangleHandle->vertex(2)->point());
          insertToHoleBoundary(highestAreaTriangleHandle->vertex(2)->point(), highestAreaTriangleHandle->vertex(1)->point());
        }
      }
    }
    numberOfHoles++;
  }while((holeEdges.size() - tempVar)>4);
}
/*End of holeDetection function*/

/*Drawing Functions*/
void drawFilledCircle(GLfloat x, GLfloat y, GLfloat radius)/* For displaying points as filled circles */
{
	int triangleAmount = 100;
	GLfloat twicePi = 2.0f * 3.14159265;
	glBegin(GL_TRIANGLE_FAN);
	glVertex2f(x, y);
	for(int i = 0; i <= triangleAmount; i++)
	glVertex2f(x + (radius * cos(i *  twicePi / triangleAmount)), y + (radius * sin(i * twicePi / triangleAmount)));
	glEnd();
}

void pointset(void)
{
  /***************Open GL commands to smoothen the line segments **************/
  glEnable (GL_LINE_SMOOTH);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  glPointSize( 6.0 );

  glLineWidth(5.0);
  glColor3f(0.0, 0.0, 0.0);
  for(int i = 0; i < holeEdges.size(); i++) /*For displaying edges of hole Boundary*/
  {
    glBegin(GL_LINES);
    glVertex2f(holeEdges[i].source.x(), holeEdges[i].source.y());
    glVertex2f(holeEdges[i].target.x(), holeEdges[i].target.y());
    glEnd();
  }

  glLineWidth(5.0);
  glColor3f(0.0, 0.0, 0.0);
  for(int i = 0; i < shape.size(); i++) /*For displaying edges of outerBoundary*/
  {
    glBegin(GL_LINES);
    glVertex2f(shape[i].source.x(),shape[i].source.y());
    glVertex2f(shape[i].target.x(),shape[i].target.y());
    glEnd();
  }

  glColor3f(1.0, 0.0, 0.0);
  float r1 = diagonalDistance/140;
  float r2 = diagonalDistance/260;

  for(int i = 0; i < inputPoints.size(); i++) /*For displaying point set*/
  {
    glColor3f(0.0, 0.0, 0.0);
    drawFilledCircle(inputPoints[i].x(), inputPoints[i].y(), r1);
    glColor3f(1.0, 0.0, 0.0);
    drawFilledCircle(inputPoints[i].x(), inputPoints[i].y(), r2);
  }
}
/*End of Drawing Functions*/

/*OpenGL Functions*/
void reshape(int w, int h)
{
	width = (GLdouble) w;
	height = (GLdouble) h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(minX-20.0, maxX+20.0, minY-20.0, maxY+20.0, -1.f, 1.f);
	return;
}

void key(unsigned char key, int x, int y)
{
	switch((char)key) {
		case 'q':
		case 27:
		{
      glutDestroyWindow(window);
  		exit(0);
    }
		default:
		{
      break;
    }
	}
	return;
}

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT);
  pointset();
  glFlush();
  return;
}
/*End of OpenGL Functions*/

int main(int argc, char **argv)
{
  std::string input = "";
  if (argc < 3)
  {
    std::cout << "Usage is ./holeDetection --in <filename> [--p <outer boundary parameter>]" << std::endl;
    std::cout << std::endl;
    std::cout << "Refer the ReadMe file for furthur details." << std::endl;
    exit(0);
  }
  else
  {
    for (int i = 1; i < argc; i++)
    {
      if (std::string(argv[i]) == "--in")
      {
        input = "data/" + std::string(argv[i + 1]) + ".xyz";
        std::cout << input << std::endl;
        i++;
      }
      else if (std::string(argv[i]) == "--p")
      {
        v = atof(argv[i+1]);
        if(v < 0 || v > 1)
        {
          v = 1;
        }
      }
    }
  }

  if(input == "")
  {
    std::cout << "Invalid Input.." << std::endl;
    exit(0);
  }

  std::ifstream inputFile(input.c_str());

  std::istream_iterator<Point> begin(inputFile);
  std::istream_iterator<Point> end;

  dt.insert(begin, end);

  for(Delaunay::Vertex_iterator vi = dt.vertices_begin(); vi != dt.vertices_end(); vi++)
  {
    if(vi->point().x() > maxX)
    {
      maxX = vi->point().x();
    }
    if(vi->point().y() > maxY)
    {
      maxY = vi->point().y();
    }
    if(vi->point().x() < minX)
    {
      minX = vi->point().x();
    }
    if(vi->point().y()<minY)
    {
      minY = vi->point().y();
    }

    inputPoints.push_back(vi->point());
  }

  diagonalDistance = distance(Point(minX, minY, 1), Point(maxX, maxY, 1));

  outerBoundary();
  holeDetection();

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
  float w = glutGet(GLUT_WINDOW_WIDTH);
  float h = glutGet(GLUT_WINDOW_HEIGHT);
  std::cout << w << " " << h << std::endl;
  //glutInitWindowSize(h, h);
  window = glutCreateWindow("Hole Detection");
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glutDisplayFunc(display);
  glutMainLoop();
  return 0;
}
