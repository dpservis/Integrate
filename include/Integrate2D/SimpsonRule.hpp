#pragma once

#include<cstddef>
#include "./DiscretizedIntegratorWrapper.hpp"
#include "../Integrate1d/SimpsonRule.hpp"

namespace integrate2d {

/** @class 
  * @brief A class that implements Simpson sums.
  * @author C.D. Clark III
  */
template<typename T>
class SimpsonRule : public DiscretizedIntegratorWrapper<integrate1d::SimpsonRule<T>>
{ 
  public:

    using BaseType = DiscretizedIntegratorWrapper<integrate1d::SimpsonRule<T>>;
    using BaseType::operator();
};



}
