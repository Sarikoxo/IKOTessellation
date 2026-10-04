#include "tessellation/private/exact.h"

#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
std::uint64_t Next(std::uint64_t& state) {
  state ^= state >> 12;
  state ^= state << 25;
  state ^= state >> 27;
  return state * UINT64_C(2685821657736338717);
}

void Check(double first, double second, std::size_t index) {
  using iko::internal::exact::Rat;
  const Rat a(first);
  const Rat b(second);
  // Each CPU operation rounds the exact result once to a double. This
  // oracle does not use the SDK's integer or rational implementation.
  volatile double x = first;
  volatile double y = second;
  const std::array<double, 4> expected{x + y, x - y, x * y, x / y};
  const std::array<double, 4> actual{(a + b).Value(), (a - b).Value(),
                                     (a * b).Value(), (a / b).Value()};
  if (actual != expected) {
    std::cerr << std::setprecision(17) << "Pair " << index << ": " << first
              << ", " << second << '\n';
    for (std::size_t operation = 0; operation < actual.size(); ++operation) {
      std::cerr << "operation " << operation << ": got " << actual[operation]
                << ", expected " << expected[operation] << '\n';
    }
    throw std::runtime_error(
        "Exact arithmetic differs from the independent CPU result");
  }
}
}  // namespace

int main() {
  try {
    if (!std::numeric_limits<double>::is_iec559 ||
        std::fegetround() != FE_TONEAREST) {
      throw std::runtime_error(
          "Exact tests require IEEE round-to-nearest doubles");
    }
    std::size_t cases = 0;
    std::uint64_t state = 901527;
    for (int i = 0; i < 1000; ++i) {
      // Fixed integer operations avoid implementation-dependent distributions.
      // Mantissas retain the original [-1, 1) and [0.5, 1) ranges; coprime
      // strides visit every exponent in [-400, 400] for both operands.
      const double first =
          std::ldexp(static_cast<double>(Next(state) >> 11), -52) - 1.;
      const double second =
          .5 + std::ldexp(static_cast<double>(Next(state) >> 11), -54);
      Check(std::ldexp(first, -400 + (37 * i) % 801),
            std::ldexp(second, -400 + (137 * i + 400) % 801), cases++);
    }
    constexpr double tiny = std::numeric_limits<double>::denorm_min();
    const std::array<std::array<double, 2>, 7> boundaries{{
        {0., 1.},                      // Zero numerator.
        {-0., -1.},                    // Signed zero and a negative divisor.
        {tiny, 1.},                    // Smallest positive subnormal.
        {tiny, 3.},                    // Underflow below the halfway point.
        {tiny, 2.},                    // Halfway to zero: round to even.
        {3 * tiny, 2.},                // Halfway between subnormals.
        {1., std::ldexp(tiny, 1000)},  // Addition below the precision of one.
    }};
    for (const auto& pair : boundaries) {
      Check(pair[0], pair[1], cases++);
    }
    std::cout << cases
              << " exact arithmetic cases passed, four operations each\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
