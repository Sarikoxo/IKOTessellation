#ifndef IKO_SRC_EXACT_H_
#define IKO_SRC_EXACT_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>

// Exact arithmetic for the baseline's bounded Bernstein/Sturm admission
// predicates. Floating input is converted to its exact binary rational value.
namespace iko::internal::exact {
struct Nat {
  std::vector<std::uint32_t> words;
  Nat() = default;
  explicit Nat(std::uint64_t x) {
    if (x) {
      words.push_back(static_cast<std::uint32_t>(x));
      if (x >> 32) {
        words.push_back(static_cast<std::uint32_t>(x >> 32));
      }
    }
  }
  void Trim() {
    while (!words.empty() && words.back() == 0) {
      words.pop_back();
    }
  }
  bool Zero() const {
    return words.empty();
  }
  int Bits() const {
    if (Zero()) {
      return 0;
    }
    auto remaining_bits = words.back();
    int b = 0;
    while (remaining_bits) {
      ++b;
      remaining_bits >>= 1;
    }
    return 32 * (static_cast<int>(words.size()) - 1) + b;
  }
  int Trailing() const {
    int n = 0;
    for (auto word_value : words) {
      if (word_value == 0) {
        n += 32;
        continue;
      }
      while (!(word_value & 1)) {
        ++n;
        word_value >>= 1;
      }
      break;
    }
    return n;
  }
  int Compare(const Nat& b) const {
    if (words.size() != b.words.size()) {
      return words.size() < b.words.size() ? -1 : 1;
    }
    for (std::size_t word_index = words.size(); word_index-- > 0;) {
      if (words[word_index] != b.words[word_index]) {
        return words[word_index] < b.words[word_index] ? -1 : 1;
      }
    }
    return 0;
  }
  Nat Shift(int count) const {
    if (Zero() || count == 0) {
      return *this;
    }
    Nat out;
    if (count > 0) {
      int n = count / 32;
      int b = count % 32;
      out.words.assign(n, 0);
      std::uint64_t carry = 0;
      for (auto source_word : words) {
        std::uint64_t value =
            (static_cast<std::uint64_t>(source_word) << b) | carry;
        out.words.push_back(static_cast<std::uint32_t>(value));
        carry = value >> 32;
      }
      if (carry) {
        out.words.push_back(static_cast<std::uint32_t>(carry));
      }
    } else {
      int n = (-count) / 32;
      int b = (-count) % 32;
      if (n >= static_cast<int>(words.size())) {
        return out;
      }
      out.words.resize(words.size() - n);
      std::uint32_t carry = 0;
      for (std::size_t source_word_index = words.size();
           source_word_index-- > static_cast<std::size_t>(n);) {
        auto x = words[source_word_index];
        out.words[source_word_index - n] = b ? (x >> b) | carry : x;
        carry = b ? x << (32 - b) : 0;
      }
      out.Trim();
    }
    return out;
  }
  static Nat Add(const Nat& a, const Nat& b) {
    Nat out;
    std::size_t n = std::max(a.words.size(), b.words.size());
    out.words.resize(n);
    std::uint64_t carry = 0;
    for (std::size_t word_index = 0; word_index < n; ++word_index) {
      std::uint64_t value =
          carry + (word_index < a.words.size() ? a.words[word_index] : 0) +
          (word_index < b.words.size()
               ? static_cast<std::uint64_t>(b.words[word_index])
               : 0);
      out.words[word_index] = static_cast<std::uint32_t>(value);
      carry = value >> 32;
    }
    if (carry) {
      out.words.push_back(static_cast<std::uint32_t>(carry));
    }
    return out;
  }
  static Nat Subtract(const Nat& a, const Nat& b) {
    Nat out = a;
    std::uint64_t borrow = 0;
    for (std::size_t word_index = 0; word_index < a.words.size();
         ++word_index) {
      std::uint64_t value =
          (word_index < b.words.size()
               ? static_cast<std::uint64_t>(b.words[word_index])
               : 0) +
          borrow;
      borrow = static_cast<std::uint64_t>(a.words[word_index]) < value;
      out.words[word_index] = static_cast<std::uint32_t>(
          static_cast<std::uint64_t>(a.words[word_index]) - value);
    }
    out.Trim();
    return out;
  }
  static Nat Multiply(const Nat& a, const Nat& b) {
    Nat out;
    if (a.Zero() || b.Zero()) {
      return out;
    }
    out.words.resize(a.words.size() + b.words.size());
    for (std::size_t first_word_index = 0; first_word_index < a.words.size();
         ++first_word_index) {
      std::uint64_t carry = 0;
      for (std::size_t second_word_index = 0;
           second_word_index < b.words.size(); ++second_word_index) {
        std::uint64_t v =
            static_cast<std::uint64_t>(a.words[first_word_index]) *
                b.words[second_word_index] +
            out.words[first_word_index + second_word_index] + carry;
        out.words[first_word_index + second_word_index] =
            static_cast<std::uint32_t>(v);
        carry = v >> 32;
      }
      out.words[first_word_index + b.words.size()] =
          static_cast<std::uint32_t>(carry);
    }
    out.Trim();
    return out;
  }
  static std::pair<Nat, Nat> Divide(const Nat& a, const Nat& b) {
    if (b.Zero()) {
      throw std::runtime_error("Exact division by zero");
    }
    if (a.Compare(b) < 0) {
      return {Nat{}, a};
    }
    int shift = a.Bits() - b.Bits();
    Nat divisor = b.Shift(shift);
    Nat remainder = a;
    Nat quotient;
    quotient.words.resize(shift / 32 + 1);
    for (int quotient_bit_index = shift; quotient_bit_index >= 0;
         --quotient_bit_index) {
      if (remainder.Compare(divisor) >= 0) {
        remainder = Subtract(remainder, divisor);
        quotient.words[quotient_bit_index / 32] |= static_cast<std::uint32_t>(1)
                                                   << (quotient_bit_index % 32);
      }
      divisor = divisor.Shift(-1);
    }
    quotient.Trim();
    return {quotient, remainder};
  }
  static Nat Gcd(Nat first_operand, Nat second_operand) {
    if (first_operand.Zero()) {
      return second_operand;
    }
    if (second_operand.Zero()) {
      return first_operand;
    }
    int common = std::min(first_operand.Trailing(), second_operand.Trailing());
    first_operand = first_operand.Shift(-first_operand.Trailing());
    do {
      second_operand = second_operand.Shift(-second_operand.Trailing());
      if (first_operand.Compare(second_operand) > 0) {
        std::swap(first_operand, second_operand);
      }
      second_operand = Subtract(second_operand, first_operand);
    } while (!second_operand.Zero());
    return first_operand.Shift(common);
  }
  std::uint64_t Small() const {
    if (words.size() > 2) {
      throw std::runtime_error("Exact conversion overflow");
    }
    return (words.empty() ? 0 : words[0]) |
           (words.size() > 1 ? static_cast<std::uint64_t>(words[1]) << 32 : 0);
  }
};
struct Rat {
  int sign = 0;
  Nat n;
  Nat d{1};
  Rat() = default;
  explicit Rat(int value)
      : sign((value > 0) - (value < 0)),
        n(value < 0
              ? static_cast<std::uint64_t>(-static_cast<std::int64_t>(value))
              : static_cast<std::uint64_t>(value)) {}
  Rat(int numerator, int denominator) : Rat(numerator) {
    *this = *this / Rat(denominator);
  }
  explicit Rat(double value) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument("Nonfinite exact input");
    }
    std::uint64_t bits;
    std::memcpy(&bits, &value, 8);
    sign = value == 0 ? 0 : (bits >> 63 ? -1 : 1);
    std::uint64_t mantissa = bits & 0xfffffffffffffull;
    int exponent = static_cast<int>((bits >> 52) & 2047);
    if (exponent) {
      mantissa |= 0x10000000000000ull;
    }
    int shift = exponent ? exponent - 1023 - 52 : -1074;
    n = Nat(mantissa);
    if (shift >= 0) {
      n = n.Shift(shift);
    } else {
      d = d.Shift(-shift);
    }
    Normalize();
  }
  void Normalize() {
    if (n.Zero()) {
      sign = 0;
      d = Nat(1);
      return;
    }
    auto g = Nat::Gcd(n, d);
    n = Nat::Divide(n, g).first;
    d = Nat::Divide(d, g).first;
  }
  Rat operator-() const {
    Rat r = *this;
    r.sign = -r.sign;
    return r;
  }
  friend Rat operator+(const Rat& a, const Rat& b) {
    if (!a.sign) {
      return b;
    }
    if (!b.sign) {
      return a;
    }
    auto g = Nat::Gcd(a.d, b.d);
    auto ad = Nat::Divide(a.d, g).first;
    auto bd = Nat::Divide(b.d, g).first;
    auto x = Nat::Multiply(a.n, bd);
    auto y = Nat::Multiply(b.n, ad);
    Rat r;
    r.d = Nat::Multiply(ad, b.d);
    if (a.sign == b.sign) {
      r.sign = a.sign;
      r.n = Nat::Add(x, y);
    } else if (x.Compare(y) >= 0) {
      r.sign = a.sign;
      r.n = Nat::Subtract(x, y);
    } else {
      r.sign = b.sign;
      r.n = Nat::Subtract(y, x);
    }
    r.Normalize();
    return r;
  }
  friend Rat operator-(const Rat& a, const Rat& b) {
    return a + (-b);
  }
  friend Rat operator*(const Rat& a, const Rat& b) {
    if (!a.sign || !b.sign) {
      return {};
    }
    auto g1 = Nat::Gcd(a.n, b.d);
    auto g2 = Nat::Gcd(b.n, a.d);
    Rat r;
    r.sign = a.sign * b.sign;
    r.n = Nat::Multiply(Nat::Divide(a.n, g1).first, Nat::Divide(b.n, g2).first);
    r.d = Nat::Multiply(Nat::Divide(a.d, g2).first, Nat::Divide(b.d, g1).first);
    return r;
  }
  friend Rat operator/(const Rat& a, const Rat& b) {
    if (!b.sign) {
      throw std::runtime_error("Exact division by zero");
    }
    Rat inverse;
    inverse.sign = b.sign;
    inverse.n = b.d;
    inverse.d = b.n;
    return a * inverse;
  }
  friend int Compare(const Rat& a, const Rat& b) {
    if (a.sign != b.sign) {
      return a.sign < b.sign ? -1 : 1;
    }
    if (!a.sign) {
      return 0;
    }
    return a.sign * Nat::Multiply(a.n, b.d).Compare(Nat::Multiply(b.n, a.d));
  }
  friend bool operator<(const Rat& a, const Rat& b) {
    return Compare(a, b) < 0;
  }
  friend bool operator>(const Rat& a, const Rat& b) {
    return Compare(a, b) > 0;
  }
  friend bool operator<=(const Rat& a, const Rat& b) {
    return Compare(a, b) <= 0;
  }
  friend bool operator>=(const Rat& a, const Rat& b) {
    return Compare(a, b) >= 0;
  }
  friend bool operator==(const Rat& a, const Rat& b) {
    return Compare(a, b) == 0;
  }
  friend bool operator!=(const Rat& a, const Rat& b) {
    return Compare(a, b) != 0;
  }
  Rat& operator+=(const Rat& b) {
    return *this = *this + b;
  }
  Rat& operator-=(const Rat& b) {
    return *this = *this - b;
  }
  Rat& operator*=(const Rat& b) {
    return *this = *this * b;
  }
  double Value() const {
    if (!sign) {
      return 0.;
    }
    int e = n.Bits() - d.Bits();
    if (e >= 0 ? n.Compare(d.Shift(e)) < 0 : n.Shift(-e).Compare(d) < 0) {
      --e;
    }
    int shift = std::min(52 - e, 1074);
    auto a = shift >= 0 ? n.Shift(shift) : n;
    auto b = shift < 0 ? d.Shift(-shift) : d;
    auto qr = Nat::Divide(a, b);
    auto q = qr.first.Small();
    int round = qr.second.Shift(1).Compare(b);
    if (round > 0 || (round == 0 && (q & 1))) {
      ++q;
    }
    return sign * std::ldexp(static_cast<double>(q), -shift);
  }
};
inline Rat Abs(Rat r) {
  if (r.sign < 0) {
    r.sign = 1;
  }
  return r;
}
using Point = std::array<Rat, 2>;
using Range = std::array<Rat, 2>;
using Contact = std::array<Rat, 2>;
}  // namespace iko::internal::exact

#endif  // IKO_SRC_EXACT_H_
