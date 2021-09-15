#pragma once

#include <gmp.h>
#include <iostream>

namespace doubleccd {

class
    Rational { // https://cs.nyu.edu/acsys/cvc3/releases/1.5/doc/rational-gmp_8cpp-source.html
public:
    mpq_t value;
    void canonicalize() { mpq_canonicalize(value); }
    int get_sign() const { return mpq_sgn(value); }
    std::string get_denominator_str()
    {
        mpz_t denominator;
        mpz_init(denominator);
        mpq_get_den(denominator, value);

        std::string v(mpz_get_str(NULL, 10, denominator));

        mpz_clear(denominator);
        return v;
    }
    std::string get_numerator_str()
    {
        mpz_t numerator;
        mpz_init(numerator);

        mpq_get_num(numerator, value);
        std::string v(mpz_get_str(NULL, 10, numerator));
        ;

        mpz_clear(numerator);
        return v;
    }
    double get_double(const std::string& num, const std::string& denom)
    {
        std::string tmp = num + "/" + denom;
        mpq_set_str(value, tmp.c_str(), 10);
        return mpq_get_d(value);
    }
    Rational()
    {
        mpq_init(value);
        mpq_set_d(value, 0);
    }

    Rational(double d)
    {
        mpq_init(value);
        mpq_set_d(value, d);
        //            canonicalize();
    }

    Rational(const mpq_t& v_)
    {
        mpq_init(value);
        mpq_set(value, v_);
        //            canonicalize();
    }

    Rational(const Rational& other)
    {
        mpq_init(value);
        mpq_set(value, other.value);
    }

    ~Rational() { mpq_clear(value); }

    friend Rational operator-(const Rational& v)
    {
        Rational r_out;
        mpq_neg(r_out.value, v.value);
        return r_out;
    }

    friend Rational operator+(const Rational& x, const Rational& y)
    {
        Rational r_out;
        mpq_add(r_out.value, x.value, y.value);
        return r_out;
    }

    friend Rational operator-(const Rational& x, const Rational& y)
    {
        Rational r_out;
        mpq_sub(r_out.value, x.value, y.value);
        return r_out;
    }

    friend Rational operator*(const Rational& x, const Rational& y)
    {
        Rational r_out;
        mpq_mul(r_out.value, x.value, y.value);
        return r_out;
    }

    friend Rational operator/(const Rational& x, const Rational& y)
    {
        Rational r_out;
        mpq_div(r_out.value, x.value, y.value);
        return r_out;
    }

    Rational& operator=(const Rational& x)
    {
        if (this == &x)
            return *this;
        mpq_set(value, x.value);
        return *this;
    }

    Rational& operator=(const double x)
    {
        mpq_set_d(value, x);
        //            canonicalize();
        return *this;
    }

    //> < ==
    friend bool operator<(const Rational& r, const Rational& r1)
    {
        return mpq_cmp(r.value, r1.value) < 0;
    }

    friend bool operator>(const Rational& r, const Rational& r1)
    {
        return mpq_cmp(r.value, r1.value) > 0;
    }

    friend bool operator<=(const Rational& r, const Rational& r1)
    {
        return mpq_cmp(r.value, r1.value) <= 0;
    }

    friend bool operator>=(const Rational& r, const Rational& r1)
    {
        return mpq_cmp(r.value, r1.value) >= 0;
    }

    friend bool operator==(const Rational& r, const Rational& r1)
    {
        return mpq_equal(r.value, r1.value);
    }

    friend bool operator!=(const Rational& r, const Rational& r1)
    {
        return !mpq_equal(r.value, r1.value);
    }

    // to double
    double to_double() { return mpq_get_d(value); }

    //<<
    friend std::ostream& operator<<(std::ostream& os, const Rational& r)
    {
        os << mpq_get_d(r.value);
        return os;
    }
};
} // namespace doubleccd
