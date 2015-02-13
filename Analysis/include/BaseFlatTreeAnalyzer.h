/*!
 * \file BaseFlatTreeAnalyzer.h
 * \brief Definition of BaseFlatTreeAnalyzer class, the base class for flat tree analyzers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-09-03 created
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

#ifndef __APPLE__
#define override
#endif

#include <iostream>
#include <cmath>
#include <set>
#include <list>
#include <locale>

#include <TColor.h>
#include <TLorentzVector.h>


#include "AnalysisBase/include/FlatEventInfo.h"
#include "AnalysisBase/include/AnalysisMath.h"
#include "AnalysisBase/include/AnalysisTypes.h"
#include "AnalysisBase/include/exception.h"
#include "AnalysisBase/include/Particles.h"
#include "PrintTools/include/RootPrintToPdf.h"

#include "MVASelections/include/MvaReader.h"

#include "Htautau_Summer13.h"
#include "AnalysisCategories.h"
#include "FlatAnalyzerDataCollection.h"

namespace analysis {

class BaseFlatTreeAnalyzer {
public:
    typedef std::map<EventRegion, PhysicalValue> PhysicalValueMap;

    virtual const EventCategorySet& EventCategoriesToProcess() const
    {
        static const EventCategorySet categories = {
            EventCategory::Inclusive, EventCategory::TwoJets_Inclusive, EventCategory::TwoJets_ZeroBtag,
            EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_OneLooseBtag,
            EventCategory::TwoJets_TwoLooseBtag
        };

        return categories;
    }

    const DataCategoryTypeSet& DataCategoryTypeToProcessForQCD() const { return dataCategoryTypeForQCD; }

    BaseFlatTreeAnalyzer(const DataCategoryCollection& _dataCategoryCollection, const std::string& _inputPath,
                         const std::string& _outputFileName, bool _applyPostFitCorrections = false)
        : inputPath(_inputPath), outputFileName(_outputFileName), dataCategoryCollection(_dataCategoryCollection),
          anaDataCollection(outputFileName + "_full.root", true),
          applyPostFitCorrections(_applyPostFitCorrections)
    {
        TH1::SetDefaultSumw2();
        gROOT->SetMustClean(kFALSE);
    }

    void Run()
    {
        std::cout << "Processing data categories... " << std::endl;
        for(const DataCategory* dataCategory : dataCategoryCollection.GetAllCategories()) {
            if(!dataCategory->sources_sf.size()) continue;
            std::cout << *dataCategory << std::endl;
            for(const auto& source_entry : dataCategory->sources_sf) {
                const std::string fullFileName = inputPath + "/" + source_entry.first;
                std::shared_ptr<TFile> file(new TFile(fullFileName.c_str(), "READ"));
                if(file->IsZombie())
                    throw exception("Input file '") << source_entry.first << "' not found.";
                std::shared_ptr<ntuple::FlatTree> tree(new ntuple::FlatTree(*file, "flatTree"));
                ProcessDataSource(*dataCategory, tree, source_entry.second);
            }
        }

        static const std::set< std::pair<std::string, EventSubCategory> > interesting_histograms = {
            { FlatAnalyzerData::m_sv_Name(), EventSubCategory::NoCuts },
            { FlatAnalyzerData::m_sv_Name(), EventSubCategory::MassWindow },
            { FlatAnalyzerData::m_ttbb_kinfit_Name(), EventSubCategory::KinematicFitConverged },
            { FlatAnalyzerData::m_ttbb_kinfit_Name(), EventSubCategory::KinematicFitConvergedWithMassWindow },
        };

        static const std::set<std::string> complete_histogram_names = {
            FlatAnalyzerData::m_sv_Name(), FlatAnalyzerData::m_ttbb_kinfit_Name()
        };

        for (const auto& hist_name : FlatAnalyzerData::GetOriginalHistogramNames<TH1D>()) {
            for(EventSubCategory subCategory : AllEventSubCategories) {
                for(EventEnergyScale energyScale : AllEventEnergyScales) {
                    if(energyScale != EventEnergyScale::Central && !complete_histogram_names.count(hist_name))
                        continue;
                    std::cout << "Processing '" << hist_name << "' in " << subCategory << "/" << energyScale
                              << "..." << std::endl;

                    std::ostringstream ss_debug;
                    std::ostream& s_out = interesting_histograms.count(std::make_pair(hist_name, subCategory))
                                        && energyScale == EventEnergyScale::Central ? std::cout : ss_debug;

                    for (EventCategory eventCategory : EventCategoriesToProcess()) {
                        const FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId(eventCategory, subCategory,
                                                                                   energyScale);
                        const auto ZTT_matched_yield = CalculateZTTmatchedYield(anaDataMetaId, hist_name, true);
                        for (const auto yield_entry : ZTT_matched_yield)
                            s_out << eventCategory << ": ZTT MC yield in " << yield_entry.first << " = "
                                  << yield_entry.second << ".\n";

                        CreateHistogramForZTT(anaDataMetaId, hist_name, ZTT_matched_yield, true);
                        CreateHistogramForZcategory(anaDataMetaId, hist_name);
                        CreateHistogramForVVcategory(anaDataMetaId, hist_name);
                    }

                    for (EventCategory eventCategory : EventCategoriesToProcess()) {
                        const FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId(eventCategory, subCategory,
                                                                                   energyScale);

                        const auto wjets_yields = CalculateWjetsYields(anaDataMetaId, hist_name, false);
                        for (const auto yield_entry : wjets_yields){
                            s_out << eventCategory << ": W+jets yield in " << yield_entry.first << " = "
                                  << yield_entry.second << ".\n";
                        }
                        EstimateWjets(anaDataMetaId, hist_name, wjets_yields);
                    }

                    for (EventCategory eventCategory : EventCategoriesToProcess()) {
                        const FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId(eventCategory, subCategory,
                                                                                   energyScale);

                        for (DataCategoryType dataCategoryType : DataCategoryTypeToProcessForQCD()) {
                            const auto qcd_yield = CalculateQCDYield(anaDataMetaId, hist_name, dataCategoryType, s_out);
                            s_out << eventCategory << ": QCD yield = " << qcd_yield << ".\n";
                            EstimateQCD(anaDataMetaId, hist_name, qcd_yield, dataCategoryType);
                        }
                        if(applyPostFitCorrections)
                            ApplyPostFitCorrections(anaDataMetaId, hist_name, false);
                        ProcessCompositDataCategories(anaDataMetaId, hist_name);
                        if(applyPostFitCorrections)
                            ApplyPostFitCorrections(anaDataMetaId, hist_name, true);
                    }
                    s_out << std::endl;
                }
            }
        }

        std::cout << "\nSaving tables... " << std::endl;
        PrintTables("comma", L",");
        PrintTables("semicolon", L";");

        std::cout << "Saving datacards... " << std::endl;
        ProduceFileForLimitsCalculation(FlatAnalyzerData::m_sv_Name(), EventSubCategory::NoCuts,
                                        &FlatAnalyzerData::m_sv);
        ProduceFileForLimitsCalculation(FlatAnalyzerData::m_ttbb_kinfit_Name(),
                                        EventSubCategory::KinematicFitConvergedWithMassWindow,
                                        &FlatAnalyzerData::m_ttbb_kinfit);

        std::cout << "Printing stacked plots... " << std::endl;
        PrintStackedPlots(EventRegion::OS_Isolated, true, true);
        PrintStackedPlots(EventRegion::OS_Isolated, false, true);

        std::cout << "Saving output file..." << std::endl;
    }

protected:
    virtual Channel ChannelId() const = 0;
    virtual EventRegion DetermineEventRegion(const ntuple::Flat& event, EventCategory eventCategory) = 0;
    virtual bool PassMvaCut(const FlatEventInfo& eventInfo, EventCategory eventCategory) { return true; }

    virtual PhysicalValue CalculateQCDYield(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                            const std::string& hist_name, DataCategoryType dataCategoryType,
                                            std::ostream& s_out) = 0;
    virtual void EstimateQCD(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId, const std::string& hist_name,
                             const PhysicalValue& scale_factor, DataCategoryType dataCategoryType) = 0;
    virtual PhysicalValueMap CalculateWjetsYields(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                                  const std::string& hist_name, bool fullEstimate) = 0;
    virtual void CreateHistogramForZTT(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                       const std::string& hist_name, const PhysicalValueMap& ztt_yield,
                                       bool useEmbedded) = 0;
    virtual void CreateHistogramForVVcategory(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                              const std::string& hist_name) = 0;

    PhysicalValueMap CalculateZTTmatchedYield(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                              const std::string& hist_name, bool useEmbedded)
    {
        const DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);

        PhysicalValueMap zttYield;
        if (useEmbedded){
            const DataCategory& embedded =  dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);
            const DataCategory& TTembedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::TT_Embedded);
            const FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_incl(EventCategory::Inclusive,
                                                                            anaDataMetaId.eventSubCategory,
                                                                            anaDataMetaId.eventEnergyScale);

            TH1D* hist_embedded_inclusive = GetSignalHistogram(anaDataMetaId_incl, embedded.name, hist_name);
            TH1D* hist_embedded_category = GetSignalHistogram(anaDataMetaId, embedded.name, hist_name);
            TH1D* hist_TTembedded_inclusive = GetSignalHistogram(anaDataMetaId_incl, TTembedded.name, hist_name);
            TH1D* hist_TTembedded_category = GetSignalHistogram(anaDataMetaId, TTembedded.name, hist_name);
            TH1D* hist_ztautau_inclusive = GetSignalHistogram(anaDataMetaId_incl, ZTT_MC.name, hist_name);
            if(!hist_embedded_inclusive)
                throw exception("embedded hist in inclusive category not found");
            if(!hist_embedded_category)
                throw exception("embedded hist not found in event category: ") << anaDataMetaId.eventCategory << "\n";
            if(!hist_TTembedded_inclusive)
                throw exception("TTembedded hist in inclusive category not found");
            if(!hist_TTembedded_category)
                throw exception("TTembedded hist not found in event category: ") << anaDataMetaId.eventCategory << "\n";
            if(!hist_ztautau_inclusive )
                throw exception("ztt hist in inclusive category not found");

            const PhysicalValue n_emb_inclusive =
                    Integral(*hist_embedded_inclusive, true) - Integral(*hist_TTembedded_inclusive, true);
            const PhysicalValue n_emb_category =
                    Integral(*hist_embedded_category, true) - Integral(*hist_TTembedded_category, true);
            const PhysicalValue n_ztautau_inclusive = Integral(*hist_ztautau_inclusive, true);
            const PhysicalValue embedded_eff = n_emb_category/n_emb_inclusive;
            zttYield[EventRegion::OS_Isolated] = n_ztautau_inclusive * embedded_eff;
        }

        for (EventRegion eventRegion : AllEventRegions){
            if (zttYield.count(eventRegion)) continue;
            TH1D* hist_ztautau = GetHistogram(anaDataMetaId, eventRegion, ZTT_MC.name, hist_name);

            if(!hist_ztautau ){
                if (eventRegion == EventRegion::OS_Isolated)
                    throw exception("ztt hist not found in event category") << anaDataMetaId.eventCategory;
                continue;
            }

            zttYield[eventRegion] = Integral(*hist_ztautau, true);
        }
        return zttYield;
    }

    virtual void CreateHistogramForZcategory(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                             const std::string& hist_name)
    {
        const std::map<DataCategoryType, DataCategoryType> z_type_category_map = {
            { DataCategoryType::ZL_MC, DataCategoryType::ZL }, { DataCategoryType::ZJ_MC, DataCategoryType::ZJ }
        };

        for (const auto& z_category : z_type_category_map){
            const DataCategory& originalZcategory = dataCategoryCollection.GetUniqueCategory(z_category.first);
            const DataCategory& newZcategory = dataCategoryCollection.GetUniqueCategory(z_category.second);

            PhysicalValueMap valueMap;

            for(EventRegion eventRegion : AllEventRegions) {
                auto z_hist_yield = GetHistogram(anaDataMetaId, eventRegion, originalZcategory.name, hist_name);
                if (z_hist_yield)
                    valueMap[eventRegion] = Integral(*z_hist_yield,true);
            }

            static const EventCategorySet categoriesToRelax = {
                EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag };
            const EventCategory shapeEventCategory = categoriesToRelax.count(anaDataMetaId.eventCategory)
                    ? MediumToLoose_EventCategoryMap.at(anaDataMetaId.eventCategory) : anaDataMetaId.eventCategory;
            const FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_shape(shapeEventCategory,
                                                                             anaDataMetaId.eventSubCategory,
                                                                             anaDataMetaId.eventEnergyScale);

            for(const auto& yield_iter : valueMap) {
                const EventRegion eventRegion = yield_iter.first;
                const PhysicalValue& yield = yield_iter.second;
                auto z_hist_shape = GetHistogram(anaDataMetaId_shape, eventRegion, originalZcategory.name, hist_name);
                if (z_hist_shape){
                    TH1D& z_hist = CloneHistogram(anaDataMetaId, eventRegion, newZcategory.name, *z_hist_shape);
                    RenormalizeHistogram(z_hist, yield, true);
                }
            }

        }
    }

    PhysicalValue CalculateYieldsForQCD(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                        EventRegion eventRegion, const std::string& hist_name, std::ostream& s_out)
    {
        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::Data);

        std::string bkg_yield_debug;
        const analysis::PhysicalValue bkg_yield =
                CalculateBackgroundIntegral(anaDataMetaId, eventRegion, hist_name, qcd.name, false, bkg_yield_debug);
        s_out << bkg_yield_debug;

        auto hist_data = GetHistogram(anaDataMetaId, eventRegion, data.name, hist_name);
        if(!hist_data)
            throw exception("Unable to find data histograms for QCD yield estimation.");
        const auto data_yield = Integral(*hist_data, true);
        const PhysicalValue yield = data_yield - bkg_yield;
        s_out << "Data yield = " << data_yield << "\nData-MC yield = " << yield << std::endl;
        if(yield.value < 0) {
            std::cout << bkg_yield_debug << "\nData yield = " << data_yield << std::endl;
            throw exception("Negative QCD yield for histogram '") << hist_name << "' in " << anaDataMetaId.eventCategory
                                                                  << " " << eventRegion << ".";
        }
        return yield;
    }

    virtual void EstimateWjets(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                               const std::string& hist_name, const PhysicalValueMap& yield_map)
    {
        static const EventCategorySet categoriesToRelax =
            { EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag };
        const EventCategory shapeEventCategory = categoriesToRelax.count(anaDataMetaId.eventCategory)
                ? MediumToLoose_EventCategoryMap.at(anaDataMetaId.eventCategory) : anaDataMetaId.eventCategory;
        return EstimateWjetsEx(anaDataMetaId, shapeEventCategory, hist_name, yield_map);
    }

    void EstimateWjetsEx(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId, EventCategory shapeEventCategory,
                         const std::string& hist_name, const PhysicalValueMap& yield_map)
    {
        const DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(DataCategoryType::WJets_MC);

        const FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_shape(shapeEventCategory,
                                                                         anaDataMetaId.eventSubCategory,
                                                                         anaDataMetaId.eventEnergyScale);

        for(const auto& yield_entry : yield_map) {
            const EventRegion eventRegion = yield_entry.first;
            const PhysicalValue& yield = yield_entry.second;
            TH1D* hist = nullptr;
            for (const DataCategory* wjets_category : wjets_mc_categories){
                if(auto hist_mc = GetHistogram(anaDataMetaId_shape, eventRegion, wjets_category->name, hist_name)) {
                    if (!hist)
                        hist = &CloneHistogram(anaDataMetaId, eventRegion, wjets.name, *hist_mc);
                    else
                        hist->Add(hist_mc);
                }
            }
            if (hist)
                RenormalizeHistogram(*hist, yield, true);
        }
    }

    const std::string& ChannelName() const { return detail::ChannelNameMap.at(ChannelId()); }
    const std::string& ChannelNameLatex() const { return detail::ChannelNameMapLatex.at(ChannelId()); }

    FlatAnalyzerData& GetAnaData(const FlatAnalyzerDataId& anaDataId)
    {
        return anaDataCollection.Get(anaDataId, ChannelId());
    }

    root_ext::SmartHistogram<TH1D>* GetHistogram(const FlatAnalyzerDataId& anaDataId, const std::string& histogramName)
    {
        return GetAnaData(anaDataId).GetPtr<TH1D>(histogramName);
    }

    root_ext::SmartHistogram<TH1D>* GetHistogram(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                                 EventRegion eventRegion, const std::string& dataCategoryName,
                                                 const std::string& histogramName)
    {
        return GetHistogram(anaDataMetaId.MakeId(eventRegion, dataCategoryName), histogramName);
    }

    root_ext::SmartHistogram<TH1D>* GetSignalHistogram(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                                       const std::string& dataCategoryName,
                                                       const std::string& histogramName)
    {
        return GetHistogram(anaDataMetaId, EventRegion::OS_Isolated, dataCategoryName, histogramName);
    }

    TH1D& CloneHistogram(const FlatAnalyzerDataId& anaDataId, const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return GetAnaData(anaDataId).Clone(originalHistogram);
    }

    TH1D& CloneHistogram(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId, EventRegion eventRegion,
                         const std::string& dataCategoryName, const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return CloneHistogram(anaDataMetaId.MakeId(eventRegion, dataCategoryName), originalHistogram);
    }

    TH1D& CloneSignalHistogram(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                               const std::string& dataCategoryName,
                               const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return CloneHistogram(anaDataMetaId, EventRegion::OS_Isolated, dataCategoryName, originalHistogram);
    }

    void ProcessDataSource(const DataCategory& dataCategory, std::shared_ptr<ntuple::FlatTree> tree,
                           double scale_factor)
    {

        static const bool applyMVAcut = false;

        const DataCategory& DYJets_incl = dataCategoryCollection.GetUniqueCategory(DataCategoryType::DYJets_incl);
        const DataCategory& DYJets_excl = dataCategoryCollection.GetUniqueCategory(DataCategoryType::DYJets_excl);
        const DataCategory& DY_Embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);

        for(Long64_t current_entry = 0; current_entry < tree->GetEntries(); ++current_entry) {
            tree->GetEntry(current_entry);
            const ntuple::Flat& event = tree->data;
            const bool useRetag = dataCategory.IsData() || dataCategory.name == DY_Embedded.name  ? false : true;
            const EventCategoryVector eventCategories = DetermineEventCategories(event.csv_Bjets,
                                                                                 event.nBjets_retagged,
                                                                                 cuts::Htautau_Summer13::btag::CSVL,
                                                                                 cuts::Htautau_Summer13::btag::CSVM,
                                                                                 useRetag);
            const bool useCustomSF = dataCategory.exclusive_sf.count(event.n_extraJets_MC);
            const double corrected_sf = useCustomSF
                    ? dataCategory.exclusive_sf.at(event.n_extraJets_MC) : scale_factor;

            if(std::isnan(event.weight)) {
                std::cerr << "ERROR: event " << event.evt << " will not be processed because event weight is NaN."
                          << std::endl;
                continue;
            }

            const double weight = event.weight * corrected_sf;
            std::shared_ptr<FlatEventInfo> eventInfo;
            for(auto eventCategory : eventCategories) {
                if (!EventCategoriesToProcess().count(eventCategory)) continue;
                const EventRegion eventRegion = DetermineEventRegion(event, eventCategory);
                if(eventRegion == EventRegion::Unknown) continue;

                if(!eventInfo)
                    eventInfo = std::shared_ptr<FlatEventInfo>(new FlatEventInfo(event, FlatEventInfo::BjetPair(0, 1),
                                                                                 false));

                UpdateMvaInfo(*eventInfo, eventCategory, false, false, false);
                if(applyMVAcut && !PassMvaCut(*eventInfo, eventCategory)) continue;

                if(dataCategory.name == DYJets_excl.name || dataCategory.name == DYJets_incl.name)
                    FillDYjetHistograms(*eventInfo, eventCategory, eventRegion, weight);

                const FlatAnalyzerDataMetaId_noSub_noES metaId_noSub_noES(eventCategory, eventRegion,
                                                                          dataCategory.name);
                const FlatAnalyzerDataMetaId_noSub metaId_noSub =
                        metaId_noSub_noES.MakeMetaId(eventInfo->eventEnergyScale);

                if (dataCategory.IsData())
                    anaDataCollection.FillAllEnergyScales(metaId_noSub_noES, ChannelId(), *eventInfo, weight);
                else if(dataCategory.name == DY_Embedded.name
                        && eventInfo->eventEnergyScale == EventEnergyScale::Central)
                    anaDataCollection.FillCentralAndJetRelatedScales(metaId_noSub_noES, ChannelId(),
                                                                     *eventInfo, weight);
                else
                    anaDataCollection.FillSubCategories(metaId_noSub, ChannelId(), *eventInfo, weight);
            }
        }
    }

    void FillDYjetHistograms(const FlatEventInfo& eventInfo, EventCategory eventCategory, EventRegion eventRegion,
                             double weight)
    {
        const DataCategory& ZL_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZL_MC);
        const DataCategory& ZJ_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZJ_MC);
        const DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);
        const DataCategory& ZTT_L = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_L);

        const std::map<ntuple::EventType, std::string> type_category_map = {
            { ntuple::EventType::ZL, ZL_MC.name }, { ntuple::EventType::ZJ, ZJ_MC.name },
            { ntuple::EventType::ZTT, ZTT_MC.name }, {ntuple::EventType::ZTT_L, ZTT_L.name}
        };

        if(type_category_map.count(eventInfo.eventType)) {
            const std::string& name = type_category_map.at(eventInfo.eventType);
            const FlatAnalyzerDataMetaId_noSub anaDataMetaId(eventCategory, eventRegion, eventInfo.eventEnergyScale,
                                                             name);
            anaDataCollection.FillSubCategories(anaDataMetaId, ChannelId(), eventInfo, weight);
        }
    }

    void UpdateMvaInfo(FlatEventInfo& eventInfo, EventCategory eventCategory, bool calc_BDT, bool calc_BDTD,
                       bool calc_BDTMitFisher)
    {
        static double const default_value = std::numeric_limits<double>::lowest();
        std::ostringstream category_name;
        category_name << eventCategory;

        auto getMVA = [&](bool calc_MVA, MVA_Selections::MvaMethod method) -> double {
            if(calc_MVA) {
                auto mvaReader = MVA_Selections::MvaReader::Get(ChannelName(), category_name.str(), method);
                if(mvaReader)
                    return mvaReader->GetMva(eventInfo.lepton_momentums.at(0), eventInfo.lepton_momentums.at(1),
                                             eventInfo.bjet_momentums.at(0), eventInfo.bjet_momentums.at(1),
                                             eventInfo.MET);
            }
            return default_value;
        };

        eventInfo.mva_BDT = getMVA(calc_BDT, MVA_Selections::BDT);
        eventInfo.mva_BDTD = getMVA(calc_BDTD, MVA_Selections::BDTD);
        eventInfo.mva_BDTMitFisher = getMVA(calc_BDTMitFisher, MVA_Selections::BDTMitFisher);
    }

    void PrintStackedPlots(EventRegion eventRegion, bool isBlind, bool drawRatio)
    {
        const std::string blindCondition = isBlind ? "_blind" : "_noBlind";
        const std::string ratioCondition = drawRatio ? "_ratio" : "_noRatio";
        std::ostringstream eventRegionName;
        eventRegionName << outputFileName << blindCondition << ratioCondition << "_" << eventRegion << ".pdf";
        root_ext::PdfPrinter printer(eventRegionName.str());

        for(EventCategory eventCategory : EventCategoriesToProcess()) {
            for (const auto& hist_name : FlatAnalyzerData::GetOriginalHistogramNames<TH1D>()) {
                for(EventSubCategory subCategory : AllEventSubCategories) {
                    const FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId(eventCategory, subCategory,
                                                                               EventEnergyScale::Central);
                    std::ostringstream ss_title;
                    ss_title << eventCategory;
                    if(subCategory != EventSubCategory::NoCuts)
                        ss_title << " " << subCategory;
                    ss_title << ": " << hist_name;
                    StackedPlotDescriptor stackDescriptor(ss_title.str(), false, ChannelNameLatex(), drawRatio);

                    for(const DataCategory* category : dataCategoryCollection.GetAllCategories()) {
                        if(!category->draw) continue;

                        const auto histogram = GetHistogram(anaDataMetaId, eventRegion, category->name, hist_name);
                        if(!histogram) continue;

                        if(category->IsSignal())
                            stackDescriptor.AddSignalHistogram(*histogram, category->title, category->color,
                                                               category->draw_sf);
                        else if(category->IsBackground())
                            stackDescriptor.AddBackgroundHistogram(*histogram, category->title, category->color);
                        else if(category->IsData())
                            stackDescriptor.AddDataHistogram(*histogram, category->title, isBlind,
                                                             GetBlindRegion(subCategory, hist_name));
                    }

                    printer.PrintStack(stackDescriptor);
                }
            }
        }
    }

    std::string FullDataCardName(const std::string& datacard_name, EventEnergyScale eventEnergyScale) const
    {
        if(eventEnergyScale == EventEnergyScale::Central)
            return datacard_name;

        std::string channel_name = ChannelName();
        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
        std::ostringstream full_name;
        full_name << datacard_name << "_CMS_scale_";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::TauDown)
            full_name << "t_" << channel_name;
        else if(eventEnergyScale == EventEnergyScale::JetUp || eventEnergyScale == EventEnergyScale::JetDown)
            full_name << "j";
        else if(eventEnergyScale == EventEnergyScale::BtagEfficiencyUp
                || eventEnergyScale == EventEnergyScale::BtagEfficiencyDown)
            full_name << "btagEff";
        else if(eventEnergyScale == EventEnergyScale::BtagFakeUp || eventEnergyScale == EventEnergyScale::BtagFakeDown)
            full_name << "btagFake";
        else
            throw exception("Unsupported event energy scale ") << eventEnergyScale;
        full_name << "_8TeV";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::JetUp ||
                eventEnergyScale == EventEnergyScale::BtagEfficiencyUp ||
                eventEnergyScale == EventEnergyScale::BtagFakeUp )
            full_name << "Up";
        else
            full_name << "Down";
        return full_name.str();
    }

    void ProduceFileForLimitsCalculation(const std::string& hist_name, EventSubCategory eventSubCategory,
                                         FlatAnalyzerData::HistogramAccessor histogramAccessor)
    {
        static const std::map<EventCategory, std::string> categoryToDirectoryNameSuffix = {
            { EventCategory::Inclusive, "inclusive" }, { EventCategory::TwoJets_ZeroBtag, "2jet0tag" },
            { EventCategory::TwoJets_OneBtag, "2jet1tag" }, { EventCategory::TwoJets_TwoBtag, "2jet2tag" }
        };

        static const std::map<std::string, std::string> channelNameForFolder = {
            { "eTau", "eleTau" }, { "muTau", "muTau" }, { "tauTau", "tauTau" }
        };

        static const double tiny_value = 1e-9;
        static const double tiny_value_error = tiny_value;

        std::ostringstream s_file_name;
        s_file_name << outputFileName << "_" << hist_name;
        if(eventSubCategory != EventSubCategory::NoCuts)
            s_file_name << "_" << eventSubCategory;
        s_file_name << ".root";
        const std::string file_name = s_file_name.str();

        std::shared_ptr<TFile> outputFile(new TFile(file_name.c_str(), "RECREATE"));
        outputFile->cd();
        for(EventCategory eventCategory : EventCategoriesToProcess()) {
            if(!categoryToDirectoryNameSuffix.count(eventCategory)) continue;
            const std::string directoryName = channelNameForFolder.at(ChannelName()) + "_"
                    + categoryToDirectoryNameSuffix.at(eventCategory);
            outputFile->mkdir(directoryName.c_str());
            outputFile->cd(directoryName.c_str());
            for(const DataCategory* dataCategory : dataCategoryCollection.GetCategories(DataCategoryType::Limits)) {
                if(!dataCategory->datacard.size())
                    throw exception("Empty datacard name for data category '") << dataCategory->name << "'.";
                for(const EventEnergyScale& eventEnergyScale : AllEventEnergyScales) {
                    const FlatAnalyzerDataMetaId_noRegion_noName meta_id(eventCategory, eventSubCategory,
                                                                         eventEnergyScale);
                    std::shared_ptr<TH1D> hist;
                    if(auto hist_orig = GetSignalHistogram(meta_id, dataCategory->name, hist_name))
                        hist = std::shared_ptr<TH1D>(static_cast<TH1D*>(hist_orig->Clone()));
                    else {
                        std::cout << "Warning - Datacard histogram '" << hist_name
                                  << "' not found for data category '" << dataCategory->name << "' in '"
                                  << eventCategory << "/" << eventSubCategory << "/" << eventEnergyScale
                                  << "'. Using histogram with a tiny yield in the central bin instead.\n";

                        FlatAnalyzerData& anaData = GetAnaData(meta_id.MakeId(EventRegion::OS_Isolated,
                                                                              dataCategory->name));
                        auto& new_hist = (anaData.*histogramAccessor)();
                        hist = std::shared_ptr<TH1D>(static_cast<TH1D*>(new_hist.Clone()));
                        const Int_t central_bin = hist->GetNbinsX() / 2;
                        hist->SetBinContent(central_bin, tiny_value);
                        hist->SetBinError(central_bin, tiny_value_error);
                    }
                    const std::string full_datacard_name = FullDataCardName(dataCategory->datacard, eventEnergyScale);
                    hist->Scale(dataCategory->limits_sf);
                    hist->Write(full_datacard_name.c_str());

                    if(eventEnergyScale == EventEnergyScale::Central && dataCategory->datacard == "ZL") {
                        std::string channel_name = ChannelName();
                        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
                        const std::string name_syst_prefix = dataCategory->datacard + "_CMS_htt_"
                                + dataCategory->datacard + "Scale_" + channel_name + "_8TeV";
                        const std::string name_syst_up = name_syst_prefix + "Up";
                        const std::string name_syst_down = name_syst_prefix + "Down";
                        std::shared_ptr<TH1D> hist_syst_up(static_cast<TH1D*>(hist->Clone()));
                        hist_syst_up->Scale(1.02);
                        hist_syst_up->Write(name_syst_up.c_str());
                        std::shared_ptr<TH1D> hist_syst_down(static_cast<TH1D*>(hist->Clone()));
                        hist_syst_down->Scale(0.98);
                        hist_syst_down->Write(name_syst_down.c_str());
                    }
                    //added shape systematics for QCD
                    if(eventEnergyScale == EventEnergyScale::Central && dataCategory->datacard == "QCD_alternative") {
                        std::string channel_name = ChannelName();
                        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
                        const std::string name_syst_prefix = "QCD_CMS_htt_" + dataCategory->datacard
                                + "Shape_" + channel_name + "_8TeV";
                        const std::string name_syst_up = name_syst_prefix + "Up";
                        const std::string name_syst_down = name_syst_prefix + "Down";
                        std::shared_ptr<TH1D> hist_syst_up(static_cast<TH1D*>(hist->Clone()));
                        hist_syst_up->Write(name_syst_up.c_str());
                        std::shared_ptr<TH1D> hist_syst_down(static_cast<TH1D*>(hist->Clone()));
                        hist_syst_down->Write(name_syst_down.c_str());
                    }
                }
            }
        }
        outputFile->Close();
    }

    void SubtractBackgroundHistograms(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                      EventRegion eventRegion, TH1D& histogram, const std::string& current_category,
                                      std::string& debug_info, std::string& negative_bins_info)
    {
        static const double correction_factor = 0.0000001;

        std::ostringstream ss_debug;

        ss_debug << "\nSubtracting background for '" << histogram.GetName() << "' in region " << eventRegion
                 << " for Event category '" << anaDataMetaId.eventCategory
                 << "' for data category '" << current_category
                 << "'.\nInitial integral: " << Integral(histogram, true) << ".\n";
        for (auto category : dataCategoryCollection.GetCategories(DataCategoryType::Background)) {
            if(category->IsComposit() || category->name == current_category || !category->isCategoryToSubtract)
                continue;

            ss_debug << "Sample '" << category->name << "': ";
            if(auto other_histogram = GetHistogram(anaDataMetaId, eventRegion, category->name, histogram.GetName())) {
                histogram.Add(other_histogram, -1);
                ss_debug << Integral(*other_histogram, true) << ".\n";
            } else
                ss_debug << "not found.\n";
        }

        const PhysicalValue original_Integral = Integral(histogram, true);
        ss_debug << "Integral after bkg subtraction: " << original_Integral.value << ".\n";
        debug_info = ss_debug.str();
        if (original_Integral.value < 0) {
            std::cout << debug_info << std::endl;
            throw exception("Integral after bkg subtraction is negative for histogram '")
                << histogram.GetName() << "' in event category " << anaDataMetaId.eventCategory
                << " for event region " << eventRegion << ".";
        }

        std::ostringstream ss_negative;

        for (Int_t n = 1; n <= histogram.GetNbinsX(); ++n) {
            if (histogram.GetBinContent(n) >= 0) continue;
            const std::string prefix = histogram.GetBinContent(n) + histogram.GetBinError(n) >= 0 ? "WARNING" : "ERROR";

            ss_negative << prefix << ": " << histogram.GetName() << " Bin " << n << ", content = "
                        << histogram.GetBinContent(n) << ", error = " << histogram.GetBinError(n)
                        << ", bin limits=[" << histogram.GetBinLowEdge(n) << "," << histogram.GetBinLowEdge(n+1)
                        << "].\n";
            const double error = correction_factor - histogram.GetBinContent(n);
            const double new_error = std::sqrt(sqr(error) + sqr(histogram.GetBinError(n)));
            histogram.SetBinContent(n, correction_factor);
            histogram.SetBinError(n, new_error);
        }
        const PhysicalValue new_Integral = Integral(histogram, true);
        const PhysicalValue correctionSF = original_Integral/new_Integral;
        histogram.Scale(correctionSF.value);
        negative_bins_info = ss_negative.str();
    }

    PhysicalValue CalculateBackgroundIntegral(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                              EventRegion eventRegion, const std::string& hist_name,
                                              const std::string& current_category,
                                              bool expect_at_least_one_contribution = false)
    {
        std::string debug_info;
        return CalculateBackgroundIntegral(anaDataMetaId, eventRegion, hist_name, current_category,
                                           expect_at_least_one_contribution, debug_info);
    }

    PhysicalValue CalculateBackgroundIntegral(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                              EventRegion eventRegion, const std::string& hist_name,
                                              const std::string& current_category,
                                              bool expect_at_least_one_contribution, std::string& debug_info)
    {
        DataCategoryPtrSet bkg_dataCategories;
        for (auto category : dataCategoryCollection.GetCategories(DataCategoryType::Background)) {
            if(category->IsComposit() || category->name == current_category || !category->isCategoryToSubtract )
                continue;
            bkg_dataCategories.insert(category);
        }

        return CalculateFullIntegral(anaDataMetaId, eventRegion, hist_name, bkg_dataCategories,
                                     expect_at_least_one_contribution, debug_info);
    }

    PhysicalValue CalculateFullIntegral(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                        EventRegion eventRegion, const std::string& hist_name,
                                        const DataCategoryPtrSet& dataCategories,
                                        bool expect_at_least_one_contribution = false)
    {
        std::string debug_info;
        return CalculateFullIntegral(anaDataMetaId, eventRegion, hist_name, dataCategories,
                                     expect_at_least_one_contribution, debug_info);
    }

    PhysicalValue CalculateFullIntegral(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                        EventRegion eventRegion, const std::string& hist_name,
                                        const DataCategoryPtrSet& dataCategories, bool expect_at_least_one_contribution,
                                        std::string& debug_info)
    {
        PhysicalValue integral;
        bool hist_found = false;

        std::ostringstream ss_debug;

        ss_debug << "\nCalculating full integral for '" << hist_name << "' in '"
                 << anaDataMetaId.MakeMetaId(eventRegion)
                 << "' considering data categories (" << dataCategories << ").\n";

        for(const auto& dataCategory : dataCategories) {
            ss_debug << "Sample '" << dataCategory->name << "': ";

            if(auto hist = GetHistogram(anaDataMetaId, eventRegion, dataCategory->name, hist_name)) {
                hist_found = true;
                const auto hist_integral = Integral(*hist, true);
                integral += hist_integral;
                ss_debug << hist_integral << ".\n";
            } else
                ss_debug << "not found.\n";
        }

        ss_debug << "Full integral: " << integral << ".\n";
        debug_info = ss_debug.str();

        if(!hist_found && expect_at_least_one_contribution) {
            std::cout << debug_info << std::endl;
            throw exception("No histogram with name '")
                << hist_name << "' was found in the given data category set (" << dataCategories
                << "), in eventCategory: '" << anaDataMetaId.MakeMetaId(eventRegion)
                << "' to calculate full integral.";
        }

        return integral;
    }

private:

    static const std::vector< std::pair<double, double> >& GetBlindRegion(EventSubCategory subCategory,
                                                                          const std::string& hist_name)
    {
        static const std::vector< std::vector< std::pair<double, double> > > blindingRegions = {
            { { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() } },
            { { 100, 150 } },
            { { 200, 400 } },
            { { 100, 150 }, { 450, 500 }, { 800, 850 }, { 1150, 1200 }, { 1500, 1550 } }
        };

        static const std::map<std::string, size_t> histogramsToBlind = {
            { FlatAnalyzerData::m_sv_Name(), 1 }, { FlatAnalyzerData::m_vis_Name(), 1 },
            { FlatAnalyzerData::m_bb_Name(), 1 }, { FlatAnalyzerData::m_ttbb_Name(), 2 },
            { FlatAnalyzerData::m_ttbb_kinfit_Name(), 2 }, { FlatAnalyzerData::m_bb_slice_Name(), 3 }
        };

        static const std::set<EventSubCategory> sidebandSubCategories = {
            EventSubCategory::OutsideMassWindow, EventSubCategory::KinematicFitConvergedOutsideMassWindow
        };

        const auto findRegionId = [&]() -> size_t {
            if(sidebandSubCategories.count(subCategory) || !histogramsToBlind.count(hist_name))
                return 0;
            return histogramsToBlind.at(hist_name);
        };

        const size_t regionId = findRegionId();
        if(regionId >= blindingRegions.size())
            throw exception("Bad blinding region index = ") << regionId;
        return blindingRegions.at(regionId);
    }

    void PrintTables(const std::string& name_suffix, const std::wstring& sep)
    {
        std::wofstream of(outputFileName + "_" + name_suffix + ".csv");

        static const std::set< std::pair<std::string, EventSubCategory> > interesting_histograms = {
            { FlatAnalyzerData::m_sv_Name(), EventSubCategory::NoCuts },
            { FlatAnalyzerData::m_sv_Name(), EventSubCategory::MassWindow },
            { FlatAnalyzerData::m_ttbb_kinfit_Name(), EventSubCategory::KinematicFitConverged },
            { FlatAnalyzerData::m_ttbb_kinfit_Name(), EventSubCategory::KinematicFitConvergedWithMassWindow }
        };

        for(const auto& hist_entry : interesting_histograms)
            PrintTables(of, sep, hist_entry.first, hist_entry.second, false, true);

        of.flush();
        of.close();
    }

    void PrintTables(std::wostream& of, const std::wstring& sep, const std::string& hist_name,
                     EventSubCategory subCategory, bool includeOverflow, bool includeError)
    {
        of << std::wstring(hist_name.begin(), hist_name.end());

        std::wstring table_name_suffix = L"";
        if(includeOverflow && includeError)
            table_name_suffix = L" with overflow and error";
        else if(includeOverflow && !includeError)
            table_name_suffix = L" with overflow";
        else if(!includeOverflow && includeError)
            table_name_suffix = L" with error";
        of << table_name_suffix << sep;

        for (EventCategory eventCategory : EventCategoriesToProcess())
            of << eventCategory << sep;
        of << std::endl;

        for (const DataCategory* dataCategory : dataCategoryCollection.GetAllCategories()) {
            of << std::wstring(dataCategory->title.begin(), dataCategory->title.end()) << sep;
            for (EventCategory eventCategory : EventCategoriesToProcess()) {
                const FlatAnalyzerDataMetaId_noRegion_noName meta_id(eventCategory, subCategory,
                                                                     EventEnergyScale::Central);

                if( TH1D* histogram = GetSignalHistogram(meta_id, dataCategory->name, hist_name) ) {
                    const PhysicalValue integral = Integral(*histogram, includeOverflow);
                    of << integral.ToString<wchar_t>(includeError) << sep;
                }
                else
                    of << "not found" << sep;
            }
            of << std::endl;
        }
        of << std::endl << std::endl;
    }

    void ProcessCompositDataCategories(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                       const std::string& hist_name)
    {
        for (analysis::EventRegion eventRegion : analysis::AllEventRegions) {
            for(const DataCategory* composit : dataCategoryCollection.GetCategories(DataCategoryType::Composit)) {
                for(const std::string& sub_name : composit->sub_categories) {
                    const DataCategory& sub_category = dataCategoryCollection.FindCategory(sub_name);
                    auto sub_hist = GetHistogram(anaDataMetaId, eventRegion, sub_category.name, hist_name);
                    if(!sub_hist) continue;
                    if(auto composit_hist = GetHistogram(anaDataMetaId, eventRegion, composit->name, hist_name))
                        composit_hist->Add(sub_hist);
                    else
                        CloneHistogram(anaDataMetaId, eventRegion, composit->name, *sub_hist);
                }
            }
        }
    }

    void ApplyPostFitCorrections(const FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                 const std::string& hist_name, bool compositFlag)
    {
        typedef std::map<std::string, double> ScaleFactorForDataCategoryMap;
        typedef std::map<EventCategory, ScaleFactorForDataCategoryMap> ScaleFactorForEventCategoryMap;
        typedef std::map<EventSubCategory, ScaleFactorForEventCategoryMap> ScaleFactorMap;
        typedef std::map<EventCategory, double> UncertaintyForEventCategoryMap;
        typedef std::map<EventSubCategory, UncertaintyForEventCategoryMap> UncertaintyMap;

        static UncertaintyMap uncertainties;
        static ScaleFactorMap scaleFactors;

        if(!scaleFactors.size()) {
            uncertainties[EventSubCategory::KinematicFitConverged] = {
                { EventCategory::TwoJets_ZeroBtag, 0.08235 },
                { EventCategory::TwoJets_OneBtag, 0.12865 },
                { EventCategory::TwoJets_TwoBtag, 0.17155 }
            };

            uncertainties[EventSubCategory::KinematicFitConvergedWithMassWindow] = {
                { EventCategory::TwoJets_ZeroBtag, 0.07615 },
                { EventCategory::TwoJets_OneBtag, 0.14487 },
                { EventCategory::TwoJets_TwoBtag, 0.18952 }
            };

            uncertainties[EventSubCategory::KinematicFitConverged][EventCategory::TwoJets_Inclusive] =
                    uncertainties[EventSubCategory::KinematicFitConverged][EventCategory::TwoJets_ZeroBtag];
            uncertainties[EventSubCategory::NoCuts] = uncertainties[EventSubCategory::KinematicFitConverged];


            scaleFactors[EventSubCategory::KinematicFitConverged][EventCategory::TwoJets_ZeroBtag] = {
                { "QCD", 0.98941 }, { "TT", 1.12182 }, { "VV", 1.12569 }, { "W", 1.01300 }, { "ZLL", 0.99556 },
                { "ZTT", 1.27644 }
            };

            scaleFactors[EventSubCategory::KinematicFitConverged][EventCategory::TwoJets_OneBtag] = {
                { "QCD", 1.03465 }, { "TT", 1.08095 }, { "VV", 1.08546 }, { "ZLL", 1.00888 }, { "ZTT", 1.19203 }
            };

            scaleFactors[EventSubCategory::KinematicFitConverged][EventCategory::TwoJets_TwoBtag] = {
                { "QCD", 0.97828 }, { "TT", 1.05903 }, { "VV", 1.04112 }, { "ZLL", 1.00152 }, { "ZTT", 1.04694 }
            };

            scaleFactors[EventSubCategory::KinematicFitConvergedWithMassWindow][EventCategory::TwoJets_ZeroBtag] = {
                { "QCD", 1.01377 }, { "TT", 1.05496 }, { "VV", 1.05401 }, { "W", 1.00883 }, { "ZLL", 1.01078 },
                { "ZTT", 1.10411 }
            };

            scaleFactors[EventSubCategory::KinematicFitConvergedWithMassWindow][EventCategory::TwoJets_OneBtag] = {
                { "QCD", 0.92723 }, { "TT", 0.98397 }, { "VV", 0.99815 }, { "ZLL", 0.99512 }, { "ZTT", 1.08041 }
            };

            scaleFactors[EventSubCategory::KinematicFitConvergedWithMassWindow][EventCategory::TwoJets_TwoBtag] = {
                { "QCD", 1.05262 }, { "TT", 1.23052 }, { "VV", 1.13869 }, { "ZLL", 0.99945 }, { "ZTT", 1.30434 }
            };

            scaleFactors[EventSubCategory::KinematicFitConverged][EventCategory::TwoJets_Inclusive] =
                    scaleFactors[EventSubCategory::KinematicFitConverged][EventCategory::TwoJets_ZeroBtag];
            scaleFactors[EventSubCategory::NoCuts] = scaleFactors[EventSubCategory::KinematicFitConverged];
        }

        const EventCategory eventCategory = anaDataMetaId.eventCategory;
        const EventSubCategory subCategory = anaDataMetaId.eventSubCategory;

        for (EventRegion eventRegion : AllEventRegions) {
            for(const DataCategory* dataCategory : dataCategoryCollection.GetCategories(DataCategoryType::Limits)) {
                if(dataCategory->IsComposit() != compositFlag) continue;
                if(!scaleFactors.count(subCategory) || !scaleFactors.at(subCategory).count(eventCategory)
                        || !scaleFactors.at(subCategory).at(eventCategory).count(dataCategory->datacard)) continue;
                const double uncertainty = uncertainties.at(subCategory).at(eventCategory);
                const double scaleFactor = scaleFactors.at(subCategory).at(eventCategory).at(dataCategory->datacard);

                if(auto hist = GetHistogram(anaDataMetaId, eventRegion, dataCategory->name, hist_name)) {
                    hist->Scale(scaleFactor);
                    for (Int_t n = 0; n <= hist->GetNbinsX() + 1; ++n) {
                        const double error = hist->GetBinError(n) * std::sqrt(1 + sqr(uncertainty));
                        hist->SetBinError(n, error);
                    }
                }
            }
        }
    }

protected:
    std::string inputPath;
    std::string outputFileName;
    DataCategoryCollection dataCategoryCollection;
    FlatAnalyzerDataCollection anaDataCollection;
    bool applyPostFitCorrections;
};

} // namespace analysis
