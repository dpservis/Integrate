#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Integrate.hpp>

using namespace Catch;
TEST_CASE("Passing integer limits")
{
  // passing integer literals to integrate functions is causing
  // integrand to evaluate to zero. probably caused by integer rounding
  // that leads to dx being zero
  auto f = [](double x) { x = 1;return x; };
  SECTION("Gaussian Quadrature")
  {
    SECTION("8'th Order")
    {
      integrate1d::gq::GaussLegendreQuadrature<double, 8> integrate;
      CHECK(integrate(f, 0, 1) == Approx(1));
    }
    SECTION("16'th Order")
    {
      integrate1d::gq::GaussLegendreQuadrature<double, 16> integrate;
      CHECK(integrate(f, 0, 1) == Approx(1));
    }
    SECTION("32'nd Order")
    {
      integrate1d::gq::GaussLegendreQuadrature<double, 32> integrate;
      CHECK(integrate(f, 0, 1) == Approx(1));
    }
    SECTION("64'th Order")
    {
      integrate1d::gq::GaussLegendreQuadrature<double, 64> integrate;
      CHECK(integrate(f, 0, 1) == Approx(1));
    }
  }
  SECTION("Riemann Sum")
  {
    integrate1d::RiemannRule<double> integrate;
    CHECK(integrate(f, 0, 1, 4) == Approx(1));
  }
  SECTION("Riemann Sum Mixed")
  {
    integrate1d::RiemannRule<double> integrate;
    CHECK(integrate(f, 0, 1., 4) == Approx(1));
  }
  SECTION("Trapezoid Sum")
  {
    integrate1d::TrapezoidRule<double> integrate;
    CHECK(integrate(f, 0, 1, 4) == Approx(1));
  }
  SECTION("Trapezoid Sum Mixed")
  {
    integrate1d::TrapezoidRule<double> integrate;
    CHECK(integrate(f, 0, 1., 4) == Approx(1));
  }
  SECTION("Simpson's Rule")
  {
    integrate1d::SimpsonRule<double> integrate;
    CHECK(integrate(f, 0, 1, 4) == Approx(1));
  }
  SECTION("Simpson's Rule Mixed")
  {
    integrate1d::SimpsonRule<double> integrate;
    CHECK(integrate(f, 0., 1, 4) == Approx(1));
  }
}
