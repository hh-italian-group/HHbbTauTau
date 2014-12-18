## @package patElectrons
#  Configuration file that defines parameters related to PAT Electron objects.
#
#  \author Konstantin Androsov (Siena University, INFN Pisa)
#  \author Maria Teresa Grippo (Siena University, INFN Pisa)
#
#  Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
#                 Maria Teresa Grippo <grippomariateresa@gmail.com>
#
#  This file is part of X->HH->bbTauTau.
#
#  X->HH->bbTauTau is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  X->HH->bbTauTau is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.

import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso

def applyElectronParameters(process, isMC):
    process.electronIsoSequence = setupPFElectronIso(process,'gsfElectrons')

    process.patElectrons.isoDeposits = cms.PSet(
        pfAllParticles   = cms.InputTag("elPFIsoDepositPUPFIso"),             # all PU   CH+MU+E
        pfChargedHadrons = cms.InputTag("elPFIsoDepositChargedPFIso"),        # all noPU CH
        pfNeutralHadrons = cms.InputTag("elPFIsoDepositNeutralPFIso"),        # all NH
        pfPhotons        = cms.InputTag("elPFIsoDepositGammaPFIso"),          # all PH
        user = cms.VInputTag( cms.InputTag("elPFIsoDepositChargedAllPFIso") ) # all noPU CH+MU+E
        )

    process.patElectrons.isolationValues = cms.PSet(
        pfAllParticles   = cms.InputTag("elPFIsoValuePU04PFIdPFIso"),
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged04PFIdPFIso"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral04PFIdPFIso"),
        pfPhotons        = cms.InputTag("elPFIsoValueGamma04PFIdPFIso"),
        user = cms.VInputTag(
            cms.InputTag("elPFIsoValueChargedAll04PFIdPFIso"),
            cms.InputTag("elPFIsoValueChargedAll04NoPFIdPFIso"),
            cms.InputTag("elPFIsoValuePU04NoPFIdPFIso"),
            cms.InputTag("elPFIsoValueCharged04NoPFIdPFIso"),
            cms.InputTag("elPFIsoValueGamma04NoPFIdPFIso"),
            cms.InputTag("elPFIsoValueNeutral04NoPFIdPFIso")
            )
        )

    process.patElectronsWithEmbeddedVariables = cms.EDProducer('ElectronsUserEmbedder',
        electronTag = cms.InputTag("patElectronsTriggerMatch"),
        vertexTag = cms.InputTag("offlinePrimaryVerticesWithBS"),
        isMC = cms.bool(isMC),
        doMVAPOG = cms.bool(True),

        inputFileName0v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat1.weights.xml'),
        inputFileName1v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat2.weights.xml'),
        inputFileName2v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat3.weights.xml'),
        inputFileName3v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat4.weights.xml'),
        inputFileName4v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat5.weights.xml'),
        inputFileName5v2 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_TrigV0_Cat6.weights.xml'),

        inputFileName0v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml'),
        inputFileName1v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml'),
        inputFileName2v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml'),
        inputFileName3v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml'),
        inputFileName4v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml'),
        inputFileName5v3 = cms.FileInPath('EgammaAnalysis/ElectronTools/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml')
    )

    process.patElectrons.embedTrack = cms.bool(True)
    process.patElectrons.embedGsfTrack = cms.bool(True)

    process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
    process.patElectrons.electronIDSources = cms.PSet(mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"))

    return
