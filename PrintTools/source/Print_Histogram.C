/*!
 * \file Print_Histogram.C
 * \brief Print histograms with specified name superimposing several files.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2014-04-04 created
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

#include "AnalysisBase/include/RootExt.h"
#include "../include/RootPrintToPdf.h"

class MyHistogramSource : public root_ext::SimpleHistogramSource<TH1D, Double_t> {
public:
    MyHistogramSource(unsigned _rebin)
        : rebinFactor(_rebin) {}
protected:
    virtual void Prepare(TH1D* histogram, const std::string& display_name,
                         const PlotOptions& plot_options) const
    {
        SimpleHistogramSource::Prepare(histogram,display_name,plot_options);
        histogram->Rebin(rebinFactor);
    }

private:
    unsigned rebinFactor;
};

class Print_Histogram {
public:
    typedef std::pair< std::string, std::string > FileTagPair;
    typedef root_ext::PdfPrinter Printer;
    //typedef root_ext::SimpleHistogramSource<TH1D, Double_t> MyHistogramSource;

    template<typename ...Args>
    Print_Histogram(const std::string& outputFileName, const std::string& _histogramName, const std::string& _title,
                    double _xMin, double _xMax, unsigned _rebin, bool _useLogX, bool _useLogY, const Args& ...args)
       : printer(outputFileName), histogramName(_histogramName), title(_title), xRange(_xMin, _xMax),
         useLogX(_useLogX), useLogY(_useLogY), source(_rebin)
    {
        Initialize(args...);
        for(const FileTagPair& fileTag : inputs) {
            auto file = root_ext::OpenRootFile(fileTag.first);
            source.Add(fileTag.second, file);
        }
    }

    void Run()
    {
        page.side.use_log_scaleX = useLogX;
        page.side.use_log_scaleY = useLogY;
        page.side.xRange = xRange;
        page.side.fit_range_x = false;
        page.side.layout.has_legend = false;
        page.side.axis_titleX = "P_{T} (#mu) [GeV]";
        page.side.axis_titleY = "N entries";
        Print(histogramName, title);
    }

private:
    template<typename ...Args>
    void Initialize(const std::string& inputName, const Args& ...args)
    {
        const size_t split_index = inputName.find_first_of(':');
        const std::string fileName = inputName.substr(0, split_index);
        const std::string tagName = inputName.substr(split_index + 1);
        inputs.push_back(FileTagPair(fileName, tagName));
        Initialize(args...);
    }

    void Initialize() {}

    void Print(const std::string& name, const std::string& title,
               std::string name_suffix = "", std::string title_suffix = "")
    {
        page.side.histogram_name = AddSuffix(name, name_suffix);
        page.title = page.side.histogram_title = AddSuffix(title, title_suffix, ". ");
        printer.Print(page, source);
    }

    static std::string AddSuffix(const std::string& name, const std::string& suffix, std::string separator = "_")
    {
        std::ostringstream full_name;
        full_name << name;
        if(suffix.size())
            full_name << separator << suffix;
        return full_name.str();
    }

private:
    std::vector<FileTagPair> inputs;
    root_ext::SingleSidedPage page;
    Printer printer;
    std::string histogramName, title;
    root_ext::Range xRange;
    bool useLogX, useLogY;
    MyHistogramSource source;
};

