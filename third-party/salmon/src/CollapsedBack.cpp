#include <vector>
#include <unordered_map>
#include <atomic>
#include <random>

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/partitioner.h"

//#include "fastapprox.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/filesystem.hpp>

// C++ string formatting library
#include "spdlog/fmt/fmt.h"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"

#include "CollapsedGibbsSampler.hpp"
#include "Transcript.hpp"
#include "TranscriptGroup.hpp"
#include "SalmonMath.hpp"
#include "AlignmentLibrary.hpp"
#include "ReadPair.hpp"
#include "UnpairedRead.hpp"
#include "ReadExperiment.hpp"
#include "MultinomialSampler.hpp"
#include "BootstrapWriter.hpp"

using BlockedIndexRange =  tbb::blocked_range<size_t>;

// intelligently chosen value adopted from
// https://github.com/pachterlab/kallisto/blob/master/src/EMAlgorithm.h#L18
constexpr double minEQClassWeight = std::numeric_limits<double>::denorm_min();
constexpr double minWeight = std::numeric_limits<double>::denorm_min();


/**
 *  Populate the prior parameters for the VBEM
 *  Note: effLens *must* be valid before calling this function.
 */
std::vector<double> populatePrior_(
                                         std::vector<Transcript>& transcripts, // transcripts
                                         Eigen::VectorXd& effLens, // current effective length estimate
                                         double priorValue,        // the per-nucleotide prior value to use
                                         bool perTranscriptPrior   // true if prior is per-txp, else per-nucleotide
                                         ) {
    // start out with the per-txp prior
    std::vector<double> priorAlphas(transcripts.size(), priorValue);

    // If the prior is per-nucleotide (default, then we need a potentially different
    // value for each transcript based on its length).
    if (!perTranscriptPrior) {
        for (size_t i = 0; i < transcripts.size(); ++i) {
            priorAlphas[i] = priorValue * effLens(i);
        }
    }
    return priorAlphas;
}


void initCountMap_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<Transcript>& transcriptsIn,
        const std::vector<double>& priorAlphas,
        MultinomialSampler& msamp,
        std::vector<uint64_t>& countMap,
        std::vector<double>& probMap,
        Eigen::VectorXd& effLens,
        std::vector<int>& txpCounts) {

    size_t offset{0};
    for (auto& eqClass : eqVec) {
        uint64_t classCount = eqClass.second.count;

        // for each transcript in this class
        const TranscriptGroup& tgroup = eqClass.first;
        const size_t groupSize = tgroup.txps.size();
        if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const auto& auxs = eqClass.second.combinedWeights;

            double denom = 0.0;
            if (BOOST_LIKELY(groupSize > 1)) {

                for (size_t i = 0; i < groupSize; ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    denom += (priorAlphas[tid] + transcriptsIn[tid].mass(false)) * aux;
                    countMap[offset + i] = 0;
                }

		if (denom > ::minEQClassWeight) {
	   	   // Get the multinomial probabilities
		   double norm = 1.0 / denom;
		   for (size_t i = 0; i < groupSize; ++i) {
		     auto tid = txps[i];
		     auto aux = auxs[i];
		     probMap[offset + i] = norm *
                        ((priorAlphas[tid] + transcriptsIn[tid].mass(false)) * aux);
		    }

	   	    // re-sample
	            msamp(countMap.begin() + offset,
                      classCount,
                      groupSize,
                      probMap.begin() + offset);
		}
            } else {
                countMap[offset] = classCount;
            }


            for (size_t i = 0; i < groupSize; ++i) {
                auto tid = txps[i];
                txpCounts[tid] += countMap[offset + i];
            }

            offset += groupSize;
       } // valid group
    } // loop over all eq classes
}

void mmseqSampleRound_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<uint64_t>& countMap,
        std::vector<double>& probMap,
        Eigen::VectorXd& effLens,
        const std::vector<double>& priorAlphas,
        std::vector<int>& txpCount,
        MultinomialSampler& msamp) {
    // offset for 2d to 1d count map
    size_t offset{0};

    //retain original txp count
    std::vector<int> origTxpCount = txpCount;

    //reset txpCounts to zero
    std::fill(txpCount.begin(), txpCount.end(), 0);

    //generate norm. coeff for \mu from \alpha (countMap)
    double norm = 0.0;
    for (size_t i = 0; i < origTxpCount.size(); ++i){
        if (effLens(i)!=0){
            norm += (origTxpCount[i]+priorAlphas[i])/effLens(i);
        }
        else{
            std::cerr<<"eff length 0; something fishy";
        }
    }

    for (auto& eqClass : eqVec) {
        // get total number of reads for an equivalence class
        uint64_t classCount = eqClass.second.count;

        // for each transcript in this class
        const TranscriptGroup& tgroup = eqClass.first;
        const size_t groupSize = tgroup.txps.size();
        if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const auto& auxs = eqClass.second.combinedWeights;

            double denom = 0.0;
            // If this is a single-transcript group,
            // then it gets the full count --- otherwise,
            // sample!
            if (BOOST_LIKELY(groupSize > 1)) {

                std::vector<uint64_t> txpResamp(groupSize);
                std::vector<double> mu(groupSize);

                // For each transcript in the group
                for (size_t i = 0; i < groupSize; ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    mu[i] = ((origTxpCount[tid]+priorAlphas[tid])/effLens(tid)) / norm;
                    denom += (priorAlphas[tid] + origTxpCount[tid]) * aux;
                }

                //calculate prob vector
                double muSum = std::accumulate(mu.begin(), mu.end(), 0.0);
                for (size_t i = 0; i < groupSize; ++i) {
                    probMap[offset + i] = mu[i]/muSum;
                    txpResamp[i] = 0;
                }

                if (denom > ::minEQClassWeight) {
                    // re-sample
                    msamp(txpResamp.begin(),                    // count array to fill in
                            classCount,                         // multinomial n
                            groupSize,		                    // multinomial k
                            probMap.begin() + offset            // where to find multinomial probs
                         );

                    for (size_t i = 0; i < groupSize; ++i) {
                        auto tid = txps.at(i);
                        txpCount.at(tid) += txpResamp.at(i);
                        //txpCount.at(tid) -= countMap.at(offset + i);
                        //countMap.at(offset + i) = txpResamp.at(i);
                    }
                }//do nothing when denom less than minEQClassWeight
                else{
                    std::cerr<<"minEQClassWeight error";
                }
            }//do nothing if group size less than 2
            else{
                 auto tid = txps.at(0);
                 txpCount.at(tid) += countMap.at(offset);
            }
            offset += groupSize;
        }// valid group
    }// loop over all eq classes
//    std::cerr<<std::accumulate(txpCount.begin(), txpCount.end(), 0.0)<<"\n";
//    std::cerr<<std::accumulate(origTxpCount.begin(), origTxpCount.end(), 0.0)<<"\n";
}


void sampleRound_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<uint64_t>& countMap,
        std::vector<double>& probMap,
        Eigen::VectorXd& effLens,
        double priorAlpha,
        std::vector<int>& txpCount,
        MultinomialSampler& msamp) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.25, 0.75);
    size_t offset{0};
    // Choose a fraction of this class to re-sample

    // The count substracted from each transcript
    std::vector<uint64_t> txpResamp;

    for (auto& eqClass : eqVec) {
        uint64_t classCount = eqClass.second.count;
        double sampleFrac = dis(gen);

        // for each transcript in this class
        const TranscriptGroup& tgroup = eqClass.first;
        const size_t groupSize = tgroup.txps.size();
        if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const auto& auxs = eqClass.second.combinedWeights;

            double denom = 0.0;
            // If this is a single-transcript group,
            // then it gets the full count --- otherwise,
            // sample!
            if (BOOST_LIKELY(groupSize > 1)) {

                // Subtract some fraction of the current equivalence
                // class' contribution from each transcript.
                uint64_t numResampled{0};
                if (groupSize > txpResamp.size()) {
                    txpResamp.resize(groupSize, 0);
                }

                // For each transcript in the group
                for (size_t i = 0; i < groupSize; ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    auto currCount = countMap[offset + i];
                    uint64_t currResamp = std::round(sampleFrac * currCount);
                    numResampled += currResamp;
                    txpResamp[i] = currResamp;
                    txpCount[tid] -= currResamp;
                    countMap[offset + i] -= currResamp;
                    denom += (priorAlpha + txpCount[tid]) * aux;
                }

                if (denom > ::minEQClassWeight) {
                    // Get the multinomial probabilities
                    double norm = 1.0 / denom;
                    for (size_t i = 0; i < groupSize; ++i) {
                        auto tid = txps[i];
                        auto aux = auxs[i];
                        probMap[offset + i] = norm * ((priorAlpha + txpCount[tid]) * aux);
                    }

                    // re-sample
                    msamp(txpResamp.begin(),        // count array to fill in
                            numResampled,		// multinomial n
                            groupSize,		// multinomial k
                            probMap.begin() + offset  // where to find multinomial probs
                         );

                    for (size_t i = 0; i < groupSize; ++i) {
                        auto tid = txps[i];
                        countMap[offset + i] += txpResamp[i];
                        txpCount[tid] += txpResamp[i];
                    }

                } else { // We didn't sample
                    // add back to txp-count!
                    for (size_t i = 0; i < groupSize; ++i) {
                        auto tid = txps[i];
                        txpCount[tid] += txpResamp[i];
                        countMap[offset + i] += txpResamp[i];
                    }
                }
            }

            offset += groupSize;
        } // valid group
    } // loop over all eq classes

}

CollapsedGibbsSampler::CollapsedGibbsSampler() {}

class DistStats {
	public:
	DistStats() : meanVal(0.0), minVal(std::numeric_limits<double>::max()), maxVal(0.0) {}
	double meanVal;
	double minVal;
	double maxVal;
};

template <typename ExpT>
bool CollapsedGibbsSampler::sample(ExpT& readExp,
        SalmonOpts& sopt,
        std::function<bool(const std::vector<int>&)>& writeBootstrap,
        uint32_t numSamples) {

    namespace bfs = boost::filesystem;
    auto& jointLog = sopt.jointLog;
    tbb::task_scheduler_init tbbScheduler(sopt.numThreads);
    std::vector<Transcript>& transcripts = readExp.transcripts();

    // Fill in the effective length vector
    Eigen::VectorXd effLens(transcripts.size());

    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
        readExp.equivalenceClassBuilder().eqVec();

    using VecT = CollapsedGibbsSampler::VecType;

    std::vector<std::vector<int>> allSamples(numSamples,
                                        std::vector<int>(transcripts.size(),0));
    double priorAlpha = 1e-8;
    bool useScaledCounts = (!sopt.useQuasi and !sopt.allowOrphans);
    auto numMappedFragments = (useScaledCounts) ? readExp.upperBoundHits() : readExp.numMappedFragments();


    for (size_t i = 0; i < transcripts.size(); ++i) {
        auto& txp = transcripts[i];
//notsure this is going to effect or not
        txp.setMass(priorAlpha + (txp.mass(false) * numMappedFragments));
        effLens(i) = txp.EffectiveLength;
    }

    bool useVBEM{sopt.useVBOpt};
    bool perTranscriptPrior{sopt.perTranscriptPrior};
    double priorValue{sopt.vbPrior};

	// If we use VBEM, we'll need the prior parameters
 	std::vector<double> priorAlphas = populatePrior_(transcripts, effLens, priorValue, perTranscriptPrior);

    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(numSamples)),
                [&eqVec, &transcripts, priorAlphas, &effLens,
                 &allSamples, &writeBootstrap, useScaledCounts,
                 &jointLog, numMappedFragments]( const BlockedIndexRange& range) -> void {


                std::random_device rd;
                MultinomialSampler ms(rd);

                size_t countMapSize{0};
                for (size_t i = 0; i < eqVec.size(); ++i) {
                    if (eqVec[i].first.valid) {
                    countMapSize += eqVec[i].first.txps.size();
                    }
                }

                size_t numTranscripts{transcripts.size()};

                // will hold estimated counts
                std::vector<int> alphas(numTranscripts, 0.0);
                std::vector<uint64_t> countMap(countMapSize, 0);
                std::vector<double> probMap(countMapSize, 0.0);

                initCountMap_(eqVec, transcripts, priorAlphas, ms, countMap, probMap, effLens, allSamples[range.begin()]);

                // For each sample this thread should generate
                bool isFirstSample{true};
                int numInternalRounds = 50;
                for (auto sampleID : boost::irange(range.begin(), range.end())) {
                    if (sampleID % 100 == 0) {
                        std::cerr << "gibbs sampling " << sampleID << "\n";
                    }

                    if (!isFirstSample) {
                        // the counts start at what they were last round.
                        allSamples[sampleID] = allSamples[sampleID-1];
                    }

                    // Thin the chain by a factor of (numInternalRounds)
                    for (size_t i = 0; i < numInternalRounds; ++i){
                        mmseqSampleRound_(eqVec, countMap, probMap, effLens, priorAlphas,
                                allSamples[sampleID], ms);
                    }

                    // If we're scaling the counts, do it here.
                    if (useScaledCounts) {
                        double numMappedFrags = static_cast<double>(numMappedFragments);
                        double alphaSum = 0.0;
                        for (auto c : allSamples[sampleID]) { alphaSum += static_cast<double>(c); }
                        if (alphaSum > ::minWeight) {
                            double scaleFrac = 1.0 / alphaSum;
                            // scaleFrac converts alpha to nucleotide fraction,
                            // and multiplying by numMappedFrags scales by the total
                            // number of mapped fragments to provide an estimated count.
                            for (size_t tn = 0; tn < numTranscripts; ++tn) {
                                alphas[tn] = static_cast<int>(
                                        std::round(
                                            numMappedFrags *
                                            (static_cast<double>(allSamples[sampleID][tn]) * scaleFrac)));
                            }
                        } else { // This shouldn't happen!
                            jointLog->error("Gibbs sampler had insufficient number of fragments!"
                                    "Something is probably wrong; please check that you "
                                    "have run salmon correctly and report this to GitHub.");
                        }
                    } else { // otherwise, just copy over from the sampled counts
                        for (size_t tn = 0; tn < numTranscripts; ++tn) {
                            alphas[tn] = static_cast<int>(allSamples[sampleID][tn]);
                        }
                    }

                    writeBootstrap(alphas);
                    //bootstrapWriter->writeBootstrap(alphas);
                    isFirstSample = false;
                }
    });
    return true;
}

template
bool CollapsedGibbsSampler::sample<ReadExperiment>(ReadExperiment& readExp,
        SalmonOpts& sopt,
        std::function<bool(const std::vector<int>&)>& writeBootstrap,
        uint32_t maxIter);

template
bool CollapsedGibbsSampler::sample<AlignmentLibrary<UnpairedRead>>(
        AlignmentLibrary<UnpairedRead>& readExp,
        SalmonOpts& sopt,
        std::function<bool(const std::vector<int>&)>& writeBootstrap,
        uint32_t maxIter);


template
bool CollapsedGibbsSampler::sample<AlignmentLibrary<ReadPair>>(
        AlignmentLibrary<ReadPair>& readExp,
        SalmonOpts& sopt,
        std::function<bool(const std::vector<int>&)>& writeBootstrap,
        uint32_t maxIter);



/*
    // Deprecated Gibbs output code
    auto numTranscripts = transcripts.size();
    std::vector<DistStats> ds(numTranscripts);

    // get posterior means
    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(numTranscripts)),
                [&allSamples, &transcripts, &ds, numMappedFragments,
                 numSamples]( const BlockedIndexRange& range) -> void {

                // For each sample this thread should generate
                for (auto tid : boost::irange(range.begin(), range.end())) {
                    double meanNumReads = {0.0};
                    for (size_t i = 0; i < numSamples; ++i) {
                      auto val = allSamples[i][tid];
                      if (val < ds[tid].minVal) { ds[tid].minVal = val; }
                      if (val > ds[tid].maxVal) { ds[tid].maxVal = val; }
                      meanNumReads += (1.0 / numSamples) * val;
                    }
        		    ds[tid].meanVal = meanNumReads;
                    transcripts[tid].setMass(ds[tid].meanVal);
                }
    });

    bfs::path gibbsSampleFile = sopt.outputDirectory / "samples.txt";
    sopt.jointLog->info("Writing posterior samples to {}", gibbsSampleFile.string());

    std::ofstream statStream(gibbsSampleFile.string());
    statStream << "# txpName\tsample_1\tsample_2\t...\tsample_n\n";

    for (size_t i = 0; i < numTranscripts; ++i) {
	    statStream << transcripts[i].RefName;
        for (size_t s = 0; s < allSamples.size(); ++s) {
            statStream << '\t' << allSamples[s][i];
        }
        statStream << '\n';
    }
    statStream.close();
    sopt.jointLog->info("done writing posterior samples");

    double cutoff = priorAlpha + 1e-8;
    // Truncate tiny expression values
    double txpSumTrunc = 0.0;
    for (size_t i = 0; i < transcripts.size(); ++i) {
	// maybe use the first decile instead of the mean for the cutoff;
	// this could let too much through
        if (transcripts[i].mass(false) <= cutoff) { transcripts[i].setMass(0.0); }
        txpSumTrunc += transcripts[i].mass(false);
    }

    for (size_t i = 0; i < transcripts.size(); ++i) {
        // Set the mass to the normalized (after truncation)
        // relative abundance
        transcripts[i].setMass(transcripts[i].mass(false) / txpSumTrunc);
    }

*/
