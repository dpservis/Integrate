#ifndef _2D_RiemannRule_hpp
#define _2D_RiemannRule_hpp

#include<functional>
#include<boost/optional.hpp>

namespace _2D {

/** @class 
  * @brief A class that implements Riemann sums.
  * @author C.D. Clark III
  */
template<typename T>
class RiemannRule
{
  public:
    RiemannRule (){};
    virtual ~RiemannRule (){};

    // This version will integrate a callable between four points
    template<typename F, typename X>
    T operator()( F f, X xa, X xb, size_t xN, X ya, X yb, size_t yN );

    // This version will integrate a set of discrete points
    template<typename X, typename Y, typename F>
    T operator()( X &x, Y &y, F &f );

  protected:
};


template<typename T>
template<typename F, typename X>
T RiemannRule<T>::operator()( F f, X xa, X xb, size_t xN, X ya, X yb, size_t yN )
{
  T sum = 0;
  X dx = (xb-xa)/xN; // make sure we don't get integer rounding
  X dy = (yb-ya)/yN; // make sure we don't get integer rounding
  X y = ya;
  for(int i = 0; i < xN; i++)
  {
    X x = xa;
    for(int j = 0; j < yN; j++)
    {
      sum += f(x,y);
      x += dx;
    }
    y += dy;
  }
  sum *= dx*dy;
  return sum;
}

template<typename T>
template<typename X, typename Y, typename F>
T RiemannRule<T>::operator()( X &x, Y &y, F &f )
{
  T sum = 0;
  for(int i = 0; i < x.size()-1; i++)
    for(int j = 0; j < y.size()-1; j++)
      sum += f[i][j]*(x[i+1]-x[i])*(y[j+1]-y[j]);

  return sum;
}

}


#endif // include protector