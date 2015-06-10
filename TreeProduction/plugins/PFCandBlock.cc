/*!
 * \file PFCandBlock.cc
 * \author Original author: Subir Sarkar
 * \author Contributing author: Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Contributing author: Maria Teresa Grippo (University of Siena, INFN Pisa)
 *
 * Copyright 2011 Subir Sarkar
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


#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "Utilities/General/interface/FileInPath.h"

#include "HHbbTauTau/TreeProduction/interface/PFCand.h"
#include "HHbbTauTau/TreeProduction/interface/TriggerTools.h"

class PFCandBlock : public edm::EDAnalyzer {
public:
    explicit PFCandBlock(const edm::ParameterSet& iConfig) :
        _inputTag(iConfig.getParameter<edm::InputTag>("srcPFCandidates")),
        pfCandTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { pfCandTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    edm::InputTag _inputTag;
    ntuple::PFCandTree pfCandTree;
};

void PFCandBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    pfCandTree.RunId() = iEvent.id().run();
    pfCandTree.LumiBlock() = iEvent.id().luminosityBlock();
    pfCandTree.EventId() = iEvent.id().event();

    edm::Handle<reco::PFCandidateCollection> pfCandHandle;
    iEvent.getByLabel(_inputTag, pfCandHandle);
    const reco::PFCandidateCollection* pfCandidates = pfCandHandle.product();

    if(!pfCandidates)
        throw std::runtime_error("pfCandidates collection not found.");

    edm::LogInfo("PFCandBlock") << "Total # PFCandidates: " << pfCandidates->size();

    for (const reco::PFCandidate& PFCand : *pfCandidates) {
        // Store Tau variables
        pfCandTree.eta()    = PFCand.eta();
        pfCandTree.phi()    = PFCand.phi();
        pfCandTree.pt()     = PFCand.pt();
        pfCandTree.mass() = PFCand.mass();
        pfCandTree.charge() = PFCand.charge();

        const reco::TrackBase* track = nullptr;
        if(PFCand.trackRef().isNonnull()) track = PFCand.trackRef().get();
        else if(PFCand.gsfTrackRef().isNonnull()) track = PFCand.gsfTrackRef().get();

        pfCandTree.haveTrackInfo() = track != nullptr;
        if(track) {
            const reco::TrackBase::Point& trk_vertex = track->vertex();
            pfCandTree.trk_vx() = trk_vertex.x();
            pfCandTree.trk_vy() = trk_vertex.y();
            pfCandTree.trk_vz() = trk_vertex.z();
        }

        pfCandTree.Fill();
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFCandBlock);
