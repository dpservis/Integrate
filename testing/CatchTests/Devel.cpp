#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Utils.hpp>

using namespace Catch;
struct A {
  void method(int) const {}
};
struct B {
  void method(std::size_t) const {}
};
struct C {
  void method(double) const {}
};

template<typename T>
bool takeInt(const T& t)
{
  return takeInt(t, integrate::priority<10>{});
}
template<typename T>
bool takeInt(const T& t, integrate::priority<1>)
{
  return false;
}
template<typename T>
auto takeInt(const T& t, integrate::priority<2>) -> decltype(t.method(integrate::IntOnly{}), true)
{
  return true;
}

TEST_CASE("Developement tests")
{
  using namespace integrate;

  CHECK(takeInt(A()));
  CHECK(takeInt(B()));
  CHECK(!takeInt(C()));
}
