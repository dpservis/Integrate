#pragma once

#include<vector>

#include "./Detail.hpp"
#include "../integrate1d/RandomAccessLambda.hpp"
#include "../integrate2d/RandomAccessLambda.hpp"

namespace integrate2d {

/** @class 
  * @brief A class that does 2D integration on discretized functions using various rules.
  * @author C.D. Clark III
  */
template<typename Integrator>
class DiscretizedIntegratorWrapper
{
  Integrator integrate;
  template<typename T>
  struct getDataType {};

  template<template<typename> class V, typename T>
  struct getDataType<V<T>> {using type = T;};

  template<template<typename,std::size_t> class V, typename T, std::size_t N>
  struct getDataType<V<T,N>> {using type = T;};

  using DataType = typename getDataType<Integrator>::type;

  public:
  template<typename ...Args>
    DiscretizedIntegratorWrapper(Args&& ...args):integrate(std::forward<Args>(args)...) {}

    template<typename X, typename Y, typename F>
    auto operator()( const X &x, const Y &y, const F &f ) const -> decltype(integrate::getSize(x),integrate::getSize(y),integrate::getElement(x,0),integrate::getElement(y,0),integrate::getElement(f,0,0),DataType())
    {
      using integrate::getSize;
      using integrate::getElement;
      std::vector<DataType> sums(getSize(x));
      for(std::size_t i = 0; i < sums.size(); ++i)
      {
        sums[i] = integrate(y, [&f,i](std::size_t j){ return getElement(f,i,j); });
      }
      return integrate(x,sums);
    }

    template<typename F>
    auto operator()( const F &f, DataType dx, DataType dy ) const -> decltype(integrate::getSizeX(f),integrate::getSizeY(f),integrate::getElement(f,0,0),DataType())
    {
      using integrate::getSizeX;
      using integrate::getSizeY;
      using integrate::getElement;
      std::vector<DataType> sums(getSizeX(f));
      for(std::size_t i = 0; i < sums.size(); ++i)
      {
        sums[i] = integrate(  integrate1d::RandomAccessLambda( [&f,i](std::size_t j){return f[i][j];}, [&f](){return integrate::getSizeY(f);}), dy);
      }
      return integrate(sums,dx);
    }

    template<typename F>
    DataType operator()( F f, DataType xa, DataType xb, std::size_t xN, DataType ya, DataType yb, std::size_t yN ) const
    {
      // discretize the function with lambda functions
      // and call the discretized function integrators
      DataType dx = (xb-xa)/xN;
      DataType dy = (yb-ya)/yN;
      return this->operator()(
          integrate1d::RandomAccessLambda([&xa,&dx](int i){return xa + i*dx;},[&xN](){return xN+1;}),
          integrate1d::RandomAccessLambda([&ya,&dy](int j){return ya + j*dy;},[&yN](){return yN+1;}),
          [&xa,&ya,&dx,&dy,&f](int i, int j){ return f(xa+i*dx,ya+j*dy); }
      );
    }

  protected:
};



}
