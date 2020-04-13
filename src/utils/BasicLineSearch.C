#include "BasicLineSearch.h"

#include "MooseUtils.h"

namespace BasicLineSearch
{

Real
brents(const Real a0,
       const Real b0,
       const Real f_a0,
       const Real f_b0,
       const std::function<Real(const Real)> & f)
{
  mooseAssert(MooseUtils::absoluteFuzzyEqual(f_a0, f(a0)), "f_a0 incorrect");
  mooseAssert(MooseUtils::absoluteFuzzyEqual(f_b0, f(b0)), "f_b0 incorrect");

  if (f_a0 * f_b0 >= 0)
    bracketedError(a0, b0, f_a0, f_b0);

  Real a = a0;
  Real b = b0;
  Real f_a = f_a0;
  Real f_b = f_b0;
  Real temp;

  if (std::abs(f_a) < std::abs(f_b))
  {
    temp = a;
    a = b;
    b = temp;
    temp = f_a;
    f_a = f_b;
    f_b = temp;
  }

  Real s;
  Real f_s;

  Real c = a;
  Real f_c = f_a;

  Real d;
  bool mflag = true;

  do
  {
    if (!MooseUtils::absoluteFuzzyEqual(f_a, f_c) && !MooseUtils::absoluteFuzzyEqual(f_b, f_c))
    {
      s = a * f_b * f_c / ((f_a - f_b) * (f_a - f_c));
      s += b * f_a * f_c / ((f_b - f_a) * (f_b - f_c));
      s += c * f_a * f_b / ((f_c - f_a) * (f_c - f_b));
    }
    else
    {
      s = b - f_b * (b - a) / (f_b - f_a);
    }

    if (((s < 0.25 * (3 * a + b) && s > b) || (s < b && s > 0.25 * (3 * a + b))) ||
        (mflag && std::abs(s - b) >= 0.5 * std::abs(b - c)) ||
        (!mflag && std::abs(s - b) >= 0.5 * std::abs(c - d)) ||
        (mflag && std::abs(b - c) < TOLERANCE * TOLERANCE) ||
        (!mflag && std::abs(c - d) < TOLERANCE * TOLERANCE))
    {
      s = 0.5 * (a + b);
      mflag = true;
    }
    else
    {
      mflag = false;
    }

    f_s = f(s);
    d = c;
    c = b;
    f_c = f_b;

    if (f_a * f_s < 0)
    {
      b = s;
      f_b = f_s;
    }
    else
    {
      a = s;
      f_a = f_s;
    }

    if (std::abs(f_a) < std::abs(f_b))
    {
      temp = a;
      a = b;
      b = temp;
      temp = f_a;
      f_a = f_b;
      f_b = temp;
    }

    if (MooseUtils::absoluteFuzzyEqual(f_b, 0) || std::abs(b - a) < TOLERANCE * TOLERANCE)
      return b;
    if (MooseUtils::absoluteFuzzyEqual(f_s, 0))
      return s;
  } while (true);
}

void
bracketedError(const Real a0, const Real b0, const Real f_a0, const Real f_b0)
{
  mooseError("Root is not bracketed in line search:\n\n",
             "  [a0, b0] = [",
             a0,
             ", ",
             b0,
             "]\n  f(a0)    = ",
             f_a0,
             "\n  f(b0)    = ",
             f_b0);
}
} // namespace BasicLineSearch
