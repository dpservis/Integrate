#pragma once

#include<cstddef>
#include "./DiscretizedIntegratorWrapper.hpp"
#include "../integrate1d/TrapezoidRule.hpp"

namespace integrate2d {

/** @class 
  * @brief A class that implements Trapezoid sums.
  * @author C.D. Clark III
  */
template<typename T>
class TrapezoidRule : public DiscretizedIntegratorWrapper<integrate1d::TrapezoidRule<T>>
{ 
  public:

    using BaseType = DiscretizedIntegratorWrapper<integrate1d::TrapezoidRule<T>>;
    using BaseType::operator();
};



}
