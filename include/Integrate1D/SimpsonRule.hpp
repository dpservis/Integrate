#pragma once
#include "./Detail.hpp"

#include <cstddef>
#include <type_traits>

namespace integrate1d
{
template <typename T>
concept VectorLike = requires(T x, size_t i) {
  
    {
        x[i]
    } -> std::same_as<double&>;
    {
        x.size()
    } -> std::convertible_to<size_t>;
};

template <typename T>
concept VectorLikeWithElement = requires(T x, size_t i) {
    {
        x[i].x
    } -> std::convertible_to<double>;
    {
        x[i].y
    } -> std::convertible_to<double>;
    {
        x.size()
    } -> std::convertible_to<size_t>;
};

template <typename T>
concept VectorLikeWithElementMethod = requires(T x, size_t i) {
    {
        x[i].x()
    } -> std::convertible_to<double>;
    {
        x[i].y()
    } -> std::convertible_to<double>;
    {
        x.size()
    } -> std::convertible_to<size_t>;
};

    template <typename T>
concept IsCallable = std::invocable<T, size_t>;

template <typename T>
concept IsCallableWithReturn = requires(T x, size_t i) {
    {
        x(i)
    } -> std::convertible_to<double>;
};

template <typename T>
concept IsCallableWithSize = requires(T x, size_t i) {
    {
        x(i)
    } -> std::convertible_to<double>;
    {
        x.size()
    } -> std::convertible_to<size_t>;
};

template <typename T>
concept ArrayLike = requires(T x, size_t i) {
    {
        x[i]
    } -> std::convertible_to<double>;
};

template <typename T>
concept IsCallableWithIndex = IsCallableWithReturn<T> or ArrayLike<T>;

template<typename T>
concept IsCallableWithIndexAndSize = IsCallableWithSize<T> or VectorLike<T>;

/** @class
 * @brief A class that implements Simpson's (1/3) rule.
 * @author C.D. Clark III
 */
template <typename T, std::size_t NN = 0>
class SimpsonRule
{
public:
    SimpsonRule() = default;

    // This version will integrate a callable between two points
    template <typename F, std::size_t NN_ = NN,
              typename SFINAE = typename std::enable_if<(NN_ == 0)>::type>
        requires std::invocable<F, T>
    T operator()(F f, T a, T b, std::size_t N) const;

    template <typename F, std::size_t NN_ = NN,
              typename SFINAE = typename std::enable_if<(NN_ > 0)>::type>
        requires std::invocable<F, T>
    T operator()(F f, T a, T b) const;

    /**
     * This version will integrate a discretized function with the x and y
     * values held in one container.
     *
     * @param x a vector-like container with operator[](int) and .size()
     * The container elements must provide .x and .y
     *
     * the container must contain 3 or more elements,
     * othersize the call is undefined behavior.
     */
    template <VectorLikeWithElement X>
    auto operator()(const X& x, long ai = 0, long bi = -1) const -> decltype(integrate::getSize(x), integrate::getElement(x, 0).x, T())
    {
        // TODO check cocnepts
        using integrate::getElement;
        using integrate::getSize;

        T sum = 0;

        auto N = getSize(x);
        if(N == 0)
            return sum;

        // support for negative indices.
        // interpret -n to mean the N-n index
        while(ai < 0)
            ai += N;
        while(bi < 0)
            bi += N;

        long i;
        for(i = ai; i < bi - 1; i += 2)
        {
            // Integrate segment using three points
            // \int_a^b f(x) dx = (b-a)/6 [ f(a) + 4*f(m) + f(b) ]
            // a = x[i]
            // b = x[i+2]
            // m = (a+b)/2
            //
            // f(a) = y[i]
            // f(b) = y[i+2]
            // and we need to interpolate f(m)
            // clang-format off
      auto x1 = getElement(x,i).x;
      auto x2 = getElement(x,i+1).x;
      auto x3 = getElement(x,i+2).x;
      auto y1 = getElement(x,i).y;
      auto y2 = getElement(x,i+1).y;
      auto y3 = getElement(x,i+2).y;
      T m = (x1 + x3)/2;
      T ym = y1*LagrangePolynomial(m, x2, x3, x1)
           + y2*LagrangePolynomial(m, x1, x3, x2)
           + y3*LagrangePolynomial(m, x1, x2, x3);

      sum += (x3-x1)/6 * (y1 + 4*ym + y3);
            // clang-format on
        }

        // if the number of elements in the range that was integrated
        // is even, then there will be one element at the end that has not
        // been integrated yet.
        if((bi + 1 - ai) % 2 == 0)
        {
            // there is one extra point at the end we need to handle
            // we will use the last *three* points to fit the polynomial
            // but then integrate between the last *two* points.
            i       = bi - 2;
            auto x1 = getElement(x, i).x;
            auto x2 = getElement(x, i + 1).x;
            auto x3 = getElement(x, i + 2).x;
            auto y1 = getElement(x, i).y;
            auto y2 = getElement(x, i + 1).y;
            auto y3 = getElement(x, i + 2).y;
            T    m  = (x2 + x3) / 2;
            T    ym = y1 * LagrangePolynomial(m, x2, x3, x1) + y2 * LagrangePolynomial(m, x1, x3, x2) + y3 * LagrangePolynomial(m, x1, x2, x3);

            sum += (x3 - x2) / 6 * (y2 + 4 * ym + y3);
        }

        return sum;
    }

    /**
     * This version will integrate a discretized function with the x and y
     * values held in one container.
     *
     * @param x a vector-like container with operator[](int) and .size()
     * The container elements must provide .x() and .y()
     *
     * the container must contain 3 or more elements,
     * othersize the call is undefined behavior.
     */
    template <VectorLikeWithElementMethod X>
    auto operator()(const X& x, long ai = 0, long bi = -1) const -> decltype(integrate::getSize(x), integrate::getElement(x, 0).x(), T())
    {
        // TODO check cocnepts
        using integrate::getElement;
        using integrate::getSize;

        T sum = 0;

        auto N = getSize(x);
        if(N == 0)
            return sum;

        // support for negative indices.
        // interpret -n to mean the N-n index
        while(ai < 0)
            ai += N;
        while(bi < 0)
            bi += N;

        long i;
        for(i = ai; i < bi - 1; i += 2)
        {
            // Integrate segment using three points
            // \int_a^b f(x) dx = (b-a)/6 [ f(a) + 4*f(m) + f(b) ]
            // a = x[i]
            // b = x[i+2]
            // m = (a+b)/2
            //
            // f(a) = y[i]
            // f(b) = y[i+2]
            // and we need to interpolate f(m)
            // clang-format off
      auto x1 = getElement(x,i).x();
      auto x2 = getElement(x,i+1).x();
      auto x3 = getElement(x,i+2).x();
      auto y1 = getElement(x,i).y();
      auto y2 = getElement(x,i+1).y();
      auto y3 = getElement(x,i+2).y();

      // Test if x1 and x2 are the same
      T m = (x1 + x3)/2;
      T ym = y1*LagrangePolynomial(m, x2, x3, x1)
           + y2*LagrangePolynomial(m, x1, x3, x2)
           + y3*LagrangePolynomial(m, x1, x2, x3);

      sum += (x3-x1)/6 * (y1 + 4*ym + y3);
            // clang-format on
        }

        // if the number of elements in the range that was integrated
        // is even, then there will be one element at the end that has not
        // been integrated yet.
        if((bi + 1 - ai) % 2 == 0)
        {
            // there is one extra point at the end we need to handle
            // we will use the last *three* points to fit the polynomial
            // but then integrate between the last *two* points.
            i       = bi - 2;
            auto x1 = getElement(x, i).x();
            auto x2 = getElement(x, i + 1).x();
            auto x3 = getElement(x, i + 2).x();
            auto y1 = getElement(x, i).y();
            auto y2 = getElement(x, i + 1).y();
            auto y3 = getElement(x, i + 2).y();
            T    m  = (x2 + x3) / 2;
            T    ym = y1 * LagrangePolynomial(m, x2, x3, x1) + y2 * LagrangePolynomial(m, x1, x3, x2) + y3 * LagrangePolynomial(m, x1, x2, x3);

            sum += (x3 - x2) / 6 * (y2 + 4 * ym + y3);
        }

        return sum;
    }

    /**
     * This version will integrate a discretized function with the x and y
     * values held in separate containers.
     *
     * @param x a vector-like container with operator[](int) and .size()
     * @param y a vector-like container with operator[](int)
     *
     * x and y must have the same size *and* contain 3 or more elements,
     * othersize the call is undefined behavior.
     */
    template <IsCallableWithIndexAndSize X, IsCallableWithIndex Y>
    auto operator()(const X& x, const Y& y, long ai = 0, long bi = -1) const -> decltype(integrate::getSize(x), integrate::getElement(x, 0), integrate::getElement(y, 0), T())
    {
        using integrate::getElement;
        using integrate::getSize;

        T sum = 0;

        auto N = getSize(x);
        if(N == 0)
            return sum;

        // support for negative indices.
        // interpret -n to mean the N-n index
        while(ai < 0)
            ai += N;
        while(bi < 0)
            bi += N;

        long i;
        for(i = ai; i < bi - 1; i += 2)
        {
            // Integrate segment using three points
            // \int_a^b f(x) dx = (b-a)/6 [ f(a) + 4*f(m) + f(b) ]
            // a = x[i]
            // b = x[i+2]
            // m = (a+b)/2
            //
            // f(a) = y[i]
            // f(b) = y[i+2]
            // and we need to interpolate f(m)
            // clang-format off
      auto x1 = getElement(x,i);
      auto x2 = getElement(x,i+1);
      auto x3 = getElement(x,i+2);
      auto y1 = getElement(y,i);
      auto y2 = getElement(y,i+1);
      auto y3 = getElement(y,i+2);
      T m = (x1 + x3)/2;
      T ym = y1*LagrangePolynomial(m, x2, x3, x1)
           + y2*LagrangePolynomial(m, x1, x3, x2)
           + y3*LagrangePolynomial(m, x1, x2, x3);

      sum += (x3-x1)/6 * (y1 + 4*ym + y3);
            // clang-format on
        }

        // if the number of elements in the range that was integrated
        // is even, then there will be one element at the end that has not
        // been integrated yet.
        if((bi + 1 - ai) % 2 == 0)
        {
            // there is one extra point at the end we need to handle
            // we will use the last *three* points to fit the polynomial
            // but then integrate between the last *two* points.
            i       = bi - 2;
            auto x1 = getElement(x, i);
            auto x2 = getElement(x, i + 1);
            auto x3 = getElement(x, i + 2);
            auto y1 = getElement(y, i);
            auto y2 = getElement(y, i + 1);
            auto y3 = getElement(y, i + 2);
            T    m  = (x2 + x3) / 2;
            T    ym = y1 * LagrangePolynomial(m, x2, x3, x1) + y2 * LagrangePolynomial(m, x1, x3, x2) + y3 * LagrangePolynomial(m, x1, x2, x3);

            sum += (x3 - x2) / 6 * (y2 + 4 * ym + y3);
        }

        return sum;
    }

    /**
     * This version will integrate a uniformly discretized function with y
     * values held in a container.
     *
     * @param y a vector-like container with operator[](int) and .size()
     * @param dx a the spacing between consecutive function values.
     *
     * y must contain 3 or more elements, othersize the call is undefined behavior.
     */
    template <IsCallableWithIndexAndSize Y>
    auto operator()(const Y& y, T dx = 1) const -> decltype(integrate::getSize(y), dx* integrate::getElement(y, 0), T())
    {
        using integrate::getElement;
        using integrate::getSize;
        T                    sum = 0;
        decltype(getSize(y)) i;

        for(i = 0; i < getSize(y) - 2; i += 2)
        {
            sum += getElement(y, i) + 4 * getElement(y, i + 1) + getElement(y, i + 2);
        }
        sum *= dx / 3.;

        // if the container size is even, there will be one extra interval we need to handle
        if(getSize(y) % 2 == 0)
        {
            // use the last three points to interpolate
            // the function to the midpoint of the last two
            // points. Then apply Simpson's rule
            // clang format off
            i    = getSize(y) - 3;
            T ym = getElement(y, i) * LagrangePolynomial(3. * dx / 2, dx, 2 * dx, 0) + getElement(y, i + 1) * LagrangePolynomial(3. * dx / 2, 0, 2 * dx, dx) + getElement(y, i + 2) * LagrangePolynomial(3. * dx / 2, 0, dx, 2 * dx);

            sum += dx / 6 * (getElement(y, i + 1) + 4 * ym + getElement(y, i + 2));
            // clang format on
        }

        return sum;
    }

protected:
    T LagrangePolynomial(T x, T A, T B, T C) const;
};

template <typename T, std::size_t NN>
template <typename F, std::size_t NN_, typename SFINAE>
    requires std::invocable<F, T>
T SimpsonRule<T, NN>::operator()(F f, T a, T b, std::size_t N) const
{
    T sum = 0;
    T dx  = static_cast<T>(b - a) / N; // size of each interval
    T x   = a;
    for(std::size_t i = 0; i < N; i++, x += dx)
    {
        sum += f(x) + 4 * f(x + dx / 2) + f(x + dx);
    }
    // note that 2*h = dx
    sum *= dx / 6;
    return sum;
}

template <typename T, std::size_t NN>
template <typename F, std::size_t NN_, typename SFINAE>
    requires std::invocable<F, T>
T SimpsonRule<T, NN>::operator()(F f, T a, T b) const
{
    T sum = 0;
    T dx  = static_cast<T>(b - a) / NN; // size of each interval
    T x   = a;
    for(std::size_t i = 0; i < NN; i++, x += dx)
    {
        sum += f(x) + 4 * f(x + dx / 2) + f(x + dx);
    }
    // note that 2*h = dx
    sum *= dx / 6;
    return sum;
}
template <typename T, std::size_t NN>
T SimpsonRule<T, NN>::LagrangePolynomial(T x, T A, T B, T C) const
{
    /// test
    auto coeff1 = (C - A);
    auto coeff2 = (C - B);

    return (x - A) * (x - B) / coeff1 / coeff2;
}

} // namespace integrate1d
