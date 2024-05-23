#pragma once

#include<cstddef>
#include "./DiscretizedIntegratorWrapper.hpp"
#include "../Integrate1D/SimpsonRule.hpp"

namespace _2D {

/** @class 
  * @brief A class that implements Simpson sums.
  * @author C.D. Clark III
  */
template<typename T>
class SimpsonRule : public DiscretizedIntegratorWrapper<_1D::SimpsonRule<T>>
{ 
  public:

    using BaseType = DiscretizedIntegratorWrapper<_1D::SimpsonRule<T>>;
    using BaseType::operator();
};



}
