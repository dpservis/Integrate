#pragma once

#include<cstddef>
#include "./DiscretizedIntegratorWrapper.hpp"
#include "../integrate1d/RiemannRule.hpp"

namespace integrate2d {

/** @class 
  * @brief A class that implements Riemann sums.
  * @author C.D. Clark III
  */
template<typename T>
class RiemannRule : public DiscretizedIntegratorWrapper<integrate1d::RiemannRule<T>>
{ 
  public:

    using BaseType = DiscretizedIntegratorWrapper<integrate1d::RiemannRule<T>>;
    using BaseType::operator();


};



}
