/*!
 * \file FlatAnalyzerData.h
 * \brief Definition of histogram containers for flat tree analyzers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-12-02 created
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

#include "AnalysisBase/include/AnalyzerData.h"
#include "AnalysisBase/include/FlatEventInfo.h"
#include "Htautau_Summer13.h"
#include "AnalysisCategories.h"

namespace analysis {

class FlatAnalyzerData : public root_ext::AnalyzerData {
public:
    TH1D_ENTRY_CUSTOM_EX(m_sv, M_tt_Bins(), "M_{#tau#tau}[GeV]", "dN/dm_{#tau#tau}[1/GeV]", false, 1.5)
    TH1D_ENTRY_CUSTOM_EX(m_vis, M_tt_Bins(), "M_{vis}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb, M_ttbb_Bins(), "M_{#tau#taubb}[GeV]", "Events", false, 1.1)
    TH1D_ENTRY_CUSTOM_EX(m_ttbb_kinfit, M_ttbb_Bins(), "M_{#tau#taubb}[GeV]", "Events", false, 1.1)


    virtual const std::vector<double>& M_tt_Bins() const
    {
        static const std::vector<double> bins = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,
                                                  160, 170, 180, 190, 200, 225, 250, 275, 300, 325, 350 };
        return bins;
    }

    virtual const std::vector<double>& M_tt_bbSlice_Bins() const
    {
        static const size_t number_of_slices = 5;
        const static std::vector<double> bins = CreateMttSlicedBins(number_of_slices, M_tt_Bins());
        return bins;
    }

    virtual const std::vector<double>& M_ttbb_Bins() const
    {
        static const std::vector<double> bins = { 200, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 500, 550,
                                                  600, 650, 700 };
        return bins;
    }

    static std::string HistogramSuffix(EventSubCategory subCategory, EventEnergyScale eventEnergyScale)
    {
        std::ostringstream ss;
        if(subCategory != EventSubCategory::NoCuts)
            ss << subCategory;
        if(eventEnergyScale != EventEnergyScale::Central) {
            if(ss.str().size())
                ss << "_";
            ss << eventEnergyScale;
        }
        return ss.str();
    }

    static std::string FullHistogramName(const std::string& hist_name, EventSubCategory subCategory,
                                         EventEnergyScale eventEnergyScale)
    {
        std::ostringstream ss;
        ss << hist_name;
        const std::string suffix = HistogramSuffix(subCategory, eventEnergyScale);
        if(suffix.size())
            ss << "_" << suffix;
        return ss.str();
    }

    virtual void Fill(const FlatEventInfo& eventInfo, double weight, EventEnergyScale eventEnergyScale,
                      EventSubCategory subCategory)
    {
        const ntuple::Flat& event = *eventInfo.event;
        const std::string key = HistogramSuffix(subCategory, eventEnergyScale);

        const double mass_tautau = event.m_sv_MC;
        m_sv(key).Fill(mass_tautau, weight);
        if (eventInfo.fitResults.has_valid_mass)
            m_ttbb_kinfit(key).Fill(eventInfo.fitResults.mass, weight);

        if (eventEnergyScale != EventEnergyScale::Central) return;
        m_vis(key).Fill(eventInfo.Htt.M(),weight);
    }

    void FillSubCategories(const FlatEventInfo& eventInfo, double weight, EventEnergyScale eventEnergyScale)
    {
        using namespace cuts::massWindow;

        const double mass_tautau = eventInfo.event->m_sv_MC;
        Fill(eventInfo, weight, eventEnergyScale, EventSubCategory::NoCuts);
        const bool inside_mass_window = mass_tautau > m_tautau_low && mass_tautau < m_tautau_high
                && eventInfo.Hbb.M() > m_bb_low && eventInfo.Hbb.M() < m_bb_high;
        if(inside_mass_window)
            Fill(eventInfo, weight, eventEnergyScale, EventSubCategory::MassWindow);
        else
            Fill(eventInfo, weight, eventEnergyScale, EventSubCategory::OutsideMassWindow);
        if(eventInfo.fitResults.has_valid_mass)
            Fill(eventInfo, weight, eventEnergyScale, EventSubCategory::KinematicFitConverged);
        if(eventInfo.fitResults.has_valid_mass && inside_mass_window)
            Fill(eventInfo, weight, eventEnergyScale, EventSubCategory::KinematicFitConvergedWithMassWindow);

        if(eventEnergyScale == EventEnergyScale::Central) {
            FillSubCategories(eventInfo, weight, EventEnergyScale::BtagEfficiencyUp);
            FillSubCategories(eventInfo, weight, EventEnergyScale::BtagEfficiencyDown);
        }
    }

    void FillEnergyScales(const FlatEventInfo& eventInfo, double weight, const std::set<EventEnergyScale>& energyScales)
    {
        for (EventEnergyScale energyScale : energyScales)
            FillSubCategories(eventInfo, weight, energyScale);
    }

    void FillCentralAndJetEnergyScales(const FlatEventInfo& eventInfo, double weight)
    {
        static const std::set<EventEnergyScale> energyScales =
            { EventEnergyScale::Central, EventEnergyScale::JetUp, EventEnergyScale::JetDown };
        FillEnergyScales(eventInfo, weight, energyScales);
    }

    void FillAllEnergyScales(const FlatEventInfo& eventInfo, double weight)
    {
        FillEnergyScales(eventInfo, weight, AllEventEnergyScales);
    }

protected:
    static void FillSlice(TH1D& hist, double m_sv, double m_Hbb, double weight)
    {
        static const std::vector<double> slice_regions = { 60, 100, 140, 200, 600 };
        static const double slice_size = 350;

        if(m_sv < 0 || m_sv >= slice_size) return;
        const auto slice_region = std::find_if(slice_regions.begin(), slice_regions.end(),
                                               [&](double x) { return m_Hbb < x; });
        if(slice_region == slice_regions.end()) return;
        const std::ptrdiff_t slice_id = slice_region - slice_regions.begin();
        const double slice_shift = slice_size * slice_id;
        hist.Fill(m_sv + slice_shift, weight);
    }

    static std::vector<double> CreateMttSlicedBins(size_t number_of_slices, const std::vector<double>& mtt_bins)
    {
        if(!mtt_bins.size())
            throw exception("Invalid mtt bins.");
        if(!number_of_slices)
            throw exception("Invalid number of slices.");

        std::vector<double> sliced_bins;
        for(size_t n = 0; n < number_of_slices; ++n) {
            size_t k = n == 0 ? 0 : 1;
            const double shift = n * mtt_bins.back();
            for(; k < mtt_bins.size(); ++k) {
                const double value = mtt_bins.at(k) + shift;
                sliced_bins.push_back(value);
            }
        }
        return sliced_bins;
    }

};

} // namespace analysis
