/*!
 * \file GenParticleBlock.cc
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


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HHbbTauTau/TreeProduction/interface/GenParticle.h"

#include <limits>

namespace {
size_t FindParticleIndex(const reco::GenParticleCollection& particles, const reco::Candidate* candidate)
    {
        const auto genCandidate = dynamic_cast<const reco::GenParticle*>(candidate);
        if(genCandidate)
        {
            for(size_t index = 0; index < particles.size(); ++index)
            {
                if(&particles[index] == genCandidate)
                    return index;
            }
        }
        throw std::logic_error("bad candidate");
    }
}

class GenParticleBlock : public edm::EDAnalyzer
{
public:
    explicit GenParticleBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _inputTag(iConfig.getParameter<edm::InputTag>("genParticleSrc")),
        genParticleTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { genParticleTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    int _verbosity;
    edm::InputTag _inputTag;
    ntuple::GenParticleTree genParticleTree;
};

void GenParticleBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //if (iEvent.isRealData()) return;

    genParticleTree.RunId() = iEvent.id().run();
    genParticleTree.LumiBlock() = iEvent.id().luminosityBlock();
    genParticleTree.EventId() = iEvent.id().event();

    edm::Handle<reco::GenParticleCollection> genParticles_handle;
    iEvent.getByLabel(_inputTag, genParticles_handle);
    
    if (!genParticles_handle.isValid()) {
        edm::LogError("GenParticleBlock") << "Error >> Failed to get GenParticleCollection for label: "
                                          << _inputTag;
        throw std::runtime_error("Failed to get GenParticleCollection for label");
    }

	const reco::GenParticleCollection& genParticles = *genParticles_handle.product();
    edm::LogInfo("GenParticleBlock") << "Total # GenParticles: " << genParticles.size();

    for(unsigned n = 0; n < genParticles.size(); ++n)
    {
        const reco::GenParticle& particle = genParticles.at(n);
        genParticleTree.PdgId() = particle.pdgId();
        genParticleTree.Status() = particle.status();
        genParticleTree.Charge() = particle.charge();
        genParticleTree.pt() = particle.pt();
        genParticleTree.eta() = particle.eta();
        genParticleTree.phi() = particle.phi();
        genParticleTree.mass() = particle.mass();
        genParticleTree.X() = particle.vx();
        genParticleTree.Y() = particle.vy();
        genParticleTree.Z() = particle.vz();
        genParticleTree.Mother_Indexes().size();
        for(size_t n = 0; n < particle.numberOfMothers(); ++n)
        {
            const size_t motherIndex = FindParticleIndex(genParticles, particle.mother(n));
            genParticleTree.Mother_Indexes().push_back(motherIndex);
        }
        genParticleTree.Fill();
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenParticleBlock);
