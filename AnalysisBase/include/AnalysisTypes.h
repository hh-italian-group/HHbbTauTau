/*!
 * \file AnalysisMath.h
 * \brief Common simple types for analysis purposes.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-10-09 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <iostream>
#include <iomanip>
#include <map>
#include <cmath>

#include "TLorentzVector.h"

#include "exception.h"
#include "Tools.h"

namespace analysis {

enum class Channel { ETau = 0, MuTau = 1, TauTau = 2 };

enum class EventEnergyScale { Central = 0, TauUp = 1, TauDown = 2, JetUp = 3, JetDown = 4, BtagEfficiencyUp = 5,
                              BtagEfficiencyDown = 6 };

namespace detail {
const std::map<Channel, std::string> ChannelNameMap = {
    { Channel::ETau, "eTau" }, { Channel::MuTau, "muTau" }, { Channel::TauTau, "tauTau" }
};

const std::map<Channel, std::string> ChannelNameMapLatex = {
    { Channel::ETau, "e#tau" }, { Channel::MuTau, "#mu#tau" }, { Channel::TauTau, "#tau#tau" }
};

const std::map<EventEnergyScale, std::string> EventEnergyScaleNameMap = {
    { EventEnergyScale::Central, "Central" },
    { EventEnergyScale::TauUp, "TauUp" }, { EventEnergyScale::TauDown, "TauDown" },
    { EventEnergyScale::JetUp, "JetUp" }, { EventEnergyScale::JetDown, "JetDown" },
    { EventEnergyScale::BtagEfficiencyUp, "BtagEfficiencyUp" },
    { EventEnergyScale::BtagEfficiencyDown, "BtagEfficiencyDown" },
};

} // namespace detail

const std::set<EventEnergyScale> AllEventEnergyScales = tools::collect_map_keys(detail::EventEnergyScaleNameMap);

std::ostream& operator<< (std::ostream& s, const Channel& c)
{
    s << detail::ChannelNameMap.at(c);
    return s;
}

std::istream& operator>> (std::istream& s, Channel& c)
{
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::ChannelNameMap) {
        if(map_entry.second == name) {
            c = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown channel name '") << name << "'.";
}

std::ostream& operator<< (std::ostream& s, const EventEnergyScale& es)
{
    s << detail::EventEnergyScaleNameMap.at(es);
    return s;
}

template<typename T>
T sqr(const T& x) { return x * x; }

namespace detail {
template<typename char_type>
struct PhysicalValueErrorSeparator;

template<>
struct PhysicalValueErrorSeparator<char> {
    static std::string Get() { return " +/- "; }
};

template<>
struct PhysicalValueErrorSeparator<wchar_t> {
    static std::wstring Get() { return L" \u00B1 "; }
};

} // namespace detail

struct PhysicalValue {

    double value;
    double error;

    PhysicalValue() : value(0), error(0) {}
    PhysicalValue(double _value, double _error) : value(_value), error(_error)
    {
        if(error < 0)
            throw exception("Negative error = ") << error << ".";
    }

    PhysicalValue& operator+=(const PhysicalValue& other)
    {
        value += other.value;
        error = std::sqrt(sqr(error) + sqr(other.error));
        return *this;
    }

    PhysicalValue operator+(const PhysicalValue& other) const
    {
        PhysicalValue result(*this);
        result += other;
        return result;
    }

    PhysicalValue& operator-=(const PhysicalValue& other)
    {
        value -= other.value;
        error = std::sqrt(sqr(error) + sqr(other.error));
        return *this;
    }

    PhysicalValue operator-(const PhysicalValue& other) const
    {
        PhysicalValue result(*this);
        result -= other;
        return result;
    }

    PhysicalValue& operator*=(const PhysicalValue& other)
    {
        value *= other.value;
        error = std::sqrt(sqr(other.value * error) + sqr(value * other.error));
        return *this;
    }

    PhysicalValue operator*(const PhysicalValue& other) const
    {
        PhysicalValue result(*this);
        result *= other;
        return result;
    }

    PhysicalValue& operator/=(const PhysicalValue& other)
    {
        value /= other.value;
        error = std::sqrt(sqr(error) + sqr(value * other.error / sqr(other.value))) / std::abs(other.value);
        return *this;
    }

    PhysicalValue operator/(const PhysicalValue& other) const
    {
        PhysicalValue result(*this);
        result /= other;
        return result;
    }

    bool operator<(const PhysicalValue& other) const { return value < other.value; }

    bool IsCompatible(const PhysicalValue& other) const
    {
        const double delta = std::abs(value - other.value);
        const double error_sum = error + other.error;
        return delta < error_sum;
    }

    template<typename char_type>
    std::basic_string<char_type> ToString(bool print_error) const
    {
        static const int number_of_significant_digits_in_error = 2;
        const int precision = error ? std::floor(std::log10(error)) - number_of_significant_digits_in_error + 1
                                    : -15;
        const double ten_pow_p = std::pow(10.0, precision);
        const double error_rounded = std::ceil(error / ten_pow_p) * ten_pow_p;
        const double value_rounded = std::round(value / ten_pow_p) * ten_pow_p;
        const int decimals_to_print = std::max(0, -precision);
        std::basic_ostringstream<char_type> ss;
        ss << std::setprecision(decimals_to_print) << std::fixed << value_rounded;
        if(print_error)
            ss << detail::PhysicalValueErrorSeparator<char_type>::Get() << error_rounded;
        return ss.str();
    }
};

std::ostream& operator<<(std::ostream& s, const PhysicalValue& v)
{
    s << v.ToString<char>(true);
    return s;
}

std::wostream& operator<<(std::wostream& s, const PhysicalValue& v)
{
    s << v.ToString<wchar_t>(true);
    return s;
}

} // namespace analysis

std::ostream& operator<<(std::ostream& s, const TVector3& v){
    s << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return s;
}

std::ostream& operator<<(std::ostream& s, const TLorentzVector& v){
    s << "(pt=" << v.Pt() << ", eta=" << v.Eta() << ", phi=" << v.Phi() << ", E=" << v.E() << ")";
    return s;
}