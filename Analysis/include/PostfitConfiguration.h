/*!
 * \file PostfitConfiguration.h
 * \brief Classes to define postfit uncertainties and scale-factors and read them from a configuration file.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2015-05-05 created
 *
 * Copyright 2015 Konstantin Androsov <konstantin.androsov@gmail.com>,
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

#include "AnalysisBase/include/ConfigReader.h"
#include "AnalysisCategories.h"

namespace analysis {

class PostfitCorrections {
public:
    typedef std::map<std::string, double> ScaleFactorMap;

    PostfitCorrections() : uncertainty(0) {}

    double GetUncertainty() const { return uncertainty; }
    void SetUncertainty(double _uncertainty) { uncertainty = _uncertainty; }

    bool HasScaleFactor(const std::string& data_name) const { return scale_factors.count(data_name); }

    double GetScaleFactor(const std::string& data_name) const
    {
        if(!HasScaleFactor(data_name))
            throw exception("Postfit SF for '") << data_name << "' not found.";
        return scale_factors.at(data_name);
    }

    void AddScaleFactor(const std::string& data_name, double sf)
    {
        if(HasScaleFactor(data_name))
            throw exception("Postfit SF for '") << data_name << "' already defined.";
        scale_factors[data_name] = sf;
    }

private:
    double uncertainty;
    ScaleFactorMap scale_factors;
};

class PostfitCorrectionsCollection {
public:
    typedef std::map<EventCategory, PostfitCorrections> CategoryCorrectionsMap;
    typedef std::map<EventSubCategory, CategoryCorrectionsMap> SubCategoryCorrectionsMap;
    typedef std::map<Channel, SubCategoryCorrectionsMap> ChannelCorrectionsMap;

    const bool HasCorrections(Channel channel, EventCategory category, EventSubCategory sub_category) const
    {
        return corrections_map.count(channel) && corrections_map.at(channel).count(sub_category)
                && corrections_map.at(channel).at(sub_category).count(category);
    }

    const PostfitCorrections& GetCorrections(Channel channel, EventCategory category,
                                             EventSubCategory sub_category) const
    {
        if(!HasCorrections(channel, category, sub_category))
            throw exception("Corrections for ") << channel << "/" << category << "/" << sub_category
                                                << " are not defined.";
        return corrections_map.at(channel).at(sub_category).at(category);
    }

    void AddCorrections(Channel channel, EventCategory category, EventSubCategory sub_category,
                        const PostfitCorrections& corrections)
    {
        if(HasCorrections(channel, category, sub_category))
            throw exception("Corrections for ") << channel << "/" << category << "/" << sub_category
                                                << " are already defined.";
        corrections_map[channel][sub_category][category] = corrections;
    }

    void CopyCorrections(Channel ref_channel, EventCategory ref_category, EventSubCategory ref_sub_category,
                         Channel channel, EventCategory category, EventSubCategory sub_category)
    {
        if(HasCorrections(channel, category, sub_category))
            throw exception("Corrections for ") << channel << "/" << category << "/" << sub_category
                                                << " are already defined.";
        if(!HasCorrections(ref_channel, ref_category, ref_sub_category))
            throw exception("Corrections for ") << ref_channel << "/" << ref_category << "/" << ref_sub_category
                                                << " are not defined.";

        corrections_map[channel][sub_category][category] =
                corrections_map.at(ref_channel).at(ref_sub_category).at(ref_category);
    }

    void CopyCorrections(Channel ref_channel, EventSubCategory ref_sub_category,
                         Channel channel, EventSubCategory sub_category)
    {
        if(corrections_map.count(channel) && corrections_map.at(channel).count(sub_category))
            throw exception("Corrections for ") << channel << "/*/" << sub_category << " are already defined.";
        if(!corrections_map.count(ref_channel) || !corrections_map.at(ref_channel).count(ref_sub_category))
            throw exception("Corrections for ") << ref_channel << "/*/" << ref_sub_category << " are not defined.";

        corrections_map[channel][sub_category] = corrections_map.at(ref_channel).at(ref_sub_category);
    }

private:
    ChannelCorrectionsMap corrections_map;
};

class PostfitCorrectionsCollectionReader : public ConfigEntryReader {
public:
    struct CorrectionsId {
        std::shared_ptr<Channel> channel;
        std::shared_ptr<EventSubCategory> sub_category;
        std::shared_ptr<EventCategory> category;

        CorrectionsId() {}

        CorrectionsId(Channel _channel, EventSubCategory _sub_category)
            : channel(new Channel(_channel)), sub_category(new EventSubCategory(_sub_category)) {}

        CorrectionsId(Channel _channel, EventSubCategory _sub_category, EventCategory _category)
            : channel(new Channel(_channel)), sub_category(new EventSubCategory(_sub_category)),
              category(new EventCategory(_category)) {}

        void SetChannel(Channel _channel) { channel = std::shared_ptr<Channel>(new Channel(_channel)); }

        void SetSubCategory(EventSubCategory _sub_category)
        {
            sub_category = std::shared_ptr<EventSubCategory>(new EventSubCategory(_sub_category));
        }

        void SetCategory(EventCategory _category)
        {
            category = std::shared_ptr<EventCategory>(new EventCategory(_category));
        }

        bool IsValid() const { return channel && sub_category; }
        bool IsFull() const { return IsValid() && category; }

        static CorrectionsId Parse(const std::string& path)
        {
            std::vector<std::string> params = ConfigReader::ParseOrderedParameterList(path, true, '/');
            if(params.size() < 2 || params.size() > 3)
                throw exception("Bad postfit corrections id.");
            std::istringstream ss_channel(params.at(0));
            Channel channel;
            ss_channel >> channel;
            std::istringstream ss_sub_category(params.at(1));
            EventSubCategory sub_category;
            ss_sub_category >> sub_category;
            if(params.size() < 3)
                return CorrectionsId(channel, sub_category);
            std::istringstream ss_category(params.at(2));
            EventCategory category;
            ss_category >> category;
            return CorrectionsId(channel, sub_category, category);
        }
    };

    PostfitCorrectionsCollectionReader(PostfitCorrectionsCollection& _output) : output(&_output) {}

    virtual void StartEntry(const std::string& /*name*/) override
    {
        current_id = CorrectionsId();
        reference_id = CorrectionsId();
        current_corrections = PostfitCorrections();
    }

    virtual void EndEntry() override
    {
        if(!current_id.IsValid())
            throw exception("Current postfit corrections id is not valid.");
        if(reference_id.IsValid()) {
            if(reference_id.IsFull()) {
                if(!current_id.IsFull())
                    throw exception("Trying to copy postfit corrections from full to not full id.");
                output->CopyCorrections(*reference_id.channel, *reference_id.category, *reference_id.sub_category,
                                        *current_id.channel, *current_id.category, *current_id.sub_category);
            } else {
                if(current_id.IsFull())
                    throw exception("Trying to copy postfit corrections from not full to full id.");
                output->CopyCorrections(*reference_id.channel, *reference_id.sub_category,
                                        *current_id.channel, *current_id.sub_category);
            }
        } else {
            if(!current_id.IsFull())
                throw exception("Unable to add not full postfit corrections id.");
            output->AddCorrections(*current_id.channel, *current_id.category, *current_id.sub_category,
                                   current_corrections);
        }
    }

    virtual void ReadParameter(const std::string& param_name, const std::string& param_value) override
    {
        std::istringstream ss(param_value);
        ss >> std::boolalpha;
        if(param_name == "channel") {
            Channel channel;
            ss >> channel;
            current_id.SetChannel(channel);
        } else if(param_name == "sub_category") {
            EventSubCategory sub_category;
            ss >> sub_category;
            current_id.SetSubCategory(sub_category);
        } else if(param_name == "category") {
            EventCategory category;
            ss >> category;
            current_id.SetCategory(category);
        } else if(param_name == "uncertainty") {
            double uncertainty;
            ss >> uncertainty;
            current_corrections.SetUncertainty(uncertainty);
        } else if(param_name == "scale_factors") {
            const std::vector<std::string> sf_vector = ConfigReader::ParseOrderedParameterList(param_value);
            for(const std::string& sf_entry : sf_vector) {
                const std::vector<std::string> sf_value = ConfigReader::ParseOrderedParameterList(sf_entry, false, ':');
                if(sf_value.size() != 2)
                    throw exception("Bad postfit sf definition '") << sf_entry << "'.";
                std::istringstream ss_sf(sf_value.at(1));
                double sf;
                ss_sf >> sf;
                current_corrections.AddScaleFactor(sf_value.at(0), sf);
            }
        } else if(param_name == "reference") {
            reference_id = CorrectionsId::Parse(param_value);
        } else
            throw exception("Unsupported parameter '") << param_name << "'.";

    }

private:
    PostfitCorrectionsCollection* output;
    CorrectionsId current_id;
    CorrectionsId reference_id;
    PostfitCorrections current_corrections;
};

} // namespace analysis
