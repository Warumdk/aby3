#pragma once

#include "Sh3Encryptor.h"
#include "Sh3Evaluator.h"
#include "Sh3Runtime.h"
#include <cryptoTools/Crypto/PRNG.h>

namespace aby3 {
    class Sh3Verifier {
    public:
        void init(u64 partyIdx, CommPkg& comm, block seed, u64 buffSize = 256) {
            mEncryptor.init(partyIdx, comm, seed, buffSize);
            mEvaluator.init(partyIdx, comm, seed, buffSize);
            mRuntime.init(partyIdx, comm);
            mPartyIdx = partyIdx;
        };

        bool verifyTripleWithOpening(CommPkg &comm, const std::array<si64, 3> &abc);

        Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<si64, 3> &abc, bool &dest);

        template<Decimal D>
        bool verifyTripleWithOpening(CommPkg &comm, const std::array<sf64<D>, 3> &abc);

        template<Decimal D>
        Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<sf64<D>, 3> &abc, bool &dest);

        bool verifyTripleWithOpening(CommPkg &comm, const std::array<sb64, 3> &abc);

        Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<sb64, 3> &abc, bool &dest);

        bool verifyTripleWithOpening(CommPkg &comm, const std::array<si64Matrix, 3> &abc);

        Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<si64Matrix, 3> &abc, bool &dest);

        template<Decimal D>
        bool verifyTripleWithOpening(CommPkg &comm, const std::array<sf64Matrix<D>, 3> &abc) {
            f64Matrix<D> a(abc[0].rows(), abc[0].cols());
            f64Matrix<D> b(abc[1].rows(), abc[1].cols());
            f64Matrix<D> c(abc[2].rows(), abc[2].cols());
            mEncryptor.revealAll(mRuntime.noDependencies(), abc[0], a).get();
            mEncryptor.revealAll(mRuntime.noDependencies(), abc[1], b).get();
            mEncryptor.revealAll(mRuntime.noDependencies(), abc[2], c).get();
            if (mPartyIdx == 0) {
                std::cout << c << std::endl <<  (a * b) << std::endl;
            }

            if (c == a * b) {
                return true;
            }

            comm.mNext.cancel();
            comm.mPrev.cancel();
            return false;
        }

        template<Decimal D>
        Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<sf64Matrix<D>, 3> &abc, bool &dest);

        bool verifyTripleWithOpening(CommPkg &comm, const std::array<sbMatrix, 3> &abc);

        Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<sbMatrix, 3> &abc, bool &dest);

        bool compareView(CommPkg &comm, i64 &x);
        bool compareView(CommPkg &comm, i64Matrix &x);

        Sh3Task compareView(Sh3Task dep, i64 &x, bool &dest);

        bool verifyTripleUsingAnother(CommPkg &comm, const std::array<si64, 3> &xyz, const std::array<si64, 3> &abc);

        bool verifyTripleUsingAnother(CommPkg &comm, const std::array<sb64, 3> &xyz, const std::array<sb64, 3> &abc);

        template<Decimal D>
        bool
        verifyTripleUsingAnother(CommPkg &comm, const std::array<sf64<D>, 3> &xyz, const std::array<sf64<D>, 3> &abc);

        bool verifyTripleUsingAnother(CommPkg &comm, const std::array<si64Matrix, 3> &xyz,
                                      const std::array<si64Matrix, 3> &abc);

        template<Decimal D>
        bool verifyTripleUsingAnother(CommPkg &comm, const std::array<sf64Matrix<D>, 3> &xyz,
                                                   const std::array<sf64Matrix<D>, 3> &abc) {

            return this->verifyTripleUsingAnother(comm, xyz, abc);
        }

        bool
        verifyTripleUsingAnother(CommPkg &comm, const std::array<sbMatrix, 3> &xyz, const std::array<si64Matrix, 3> &abc);


        Sh3Task verifyTripleUsingAnother(Sh3Task dep, const std::array<si64, 3> &xyz, const std::array<si64, 3> &abc,
                                         bool &dest);

        Sh3Task verifyTripleUsingAnother(Sh3Task dep, const std::array<sb64, 3> &xyz, const std::array<sb64, 3> &abc,
                                         bool &dest);

        template<Decimal D>
        Sh3Task
        verifyTripleUsingAnother(Sh3Task dep, const std::array<sf64<D>, 3> &xyz, const std::array<sf64<D>, 3> &abc,
                                 bool &dest);

        Sh3Task verifyTripleUsingAnother(Sh3Task dep, const std::array<si64Matrix, 3> &xyz,
                                         const std::array<si64Matrix, 3> &abc,
                                         bool &dest);

        template<Decimal D>
        Sh3Task verifyTripleUsingAnother(Sh3Task dep, const std::array<sf64Matrix<D>, 3> &xyz,
                                         const std::array<sf64Matrix<D>, 3> &abc,
                                         bool &dest);

        Sh3Task
        verifyTripleUsingAnother(Sh3Task dep, const std::array<sbMatrix, 3> &xyz, const std::array<si64Matrix, 3> &abc,
                                 bool &dest);

        template<typename T>
        std::vector<T> perm(CommPkg &comm, std::vector<T> &d) {
            std::vector<i64> indices = this->coin(comm, d.size());
            for (int j = d.size()-1; j >= 0; --j) {
                std::swap(d[j], d[indices[j] % d.size()]);
            }
            return d;
        }

        template<typename T>
        Sh3Task perm(Sh3Task dep, std::vector<T> &d, std::vector<T> &dest);

        std::vector<std::array<si64, 3>> generateTriples(CommPkg &comm, int N, int B, int C);

        std::vector<std::array<si64Matrix, 3>> generateTriples(CommPkg &comm, int N, int B, int C, int rowsA, int colsA, int rowsB, int colsB);

        template<Decimal D>
        std::vector<std::array<sf64Matrix<D>, 3>> generateTriples(CommPkg &comm, int N, int B, int C, int rowsA, int colsA, int rowsB, int colsB) {
            // 1) Generate Random Sharings
            int M = N*B+C;
            std::vector<std::array<sf64Matrix<D>, 3>> triples;
            for (int i = 0; i < M; ++i) {
                sf64Matrix<D> aMatrix(rowsA, colsA);
                sf64Matrix<D> bMatrix(rowsB, colsB);
                sf64Matrix<D> res;

                this->mEncryptor.rand(aMatrix.i64Cast());
                this->mEncryptor.rand(bMatrix.i64Cast());
                mEvaluator.asyncMul(mRuntime, aMatrix, bMatrix, res).get();
                std::array<sf64Matrix<D>, 3> triple = {aMatrix, bMatrix, res};
                triples.emplace_back(triple);
            }
            if (mPartyIdx == 0) {
                std::cout << "1" << std::endl;
            }

            // 3) Cut-and-bucket
            triples = this->perm(comm, triples);
            if (mPartyIdx == 0) {
                std::cout << "1.5" << std::endl;
            }
            for (int i = 0; i < C; ++i) {
                if (!this->verifyTripleWithOpening<D>(comm, triples.back())) {
                    return std::vector<std::array<sf64Matrix<D>, 3>>{};
                }
                triples.pop_back();
            }
            if (mPartyIdx == 0) {
                std::cout << "2" << std::endl;
            }
            std::vector<std::vector<std::array<sf64Matrix<D>, 3>>> buckets;
            for (int i = 0; i < N; ++i) {
                std::vector<std::array<sf64Matrix<D>, 3>> bucket;
                for (int j = 0; j < B; ++j) {
                    bucket.emplace_back(triples.back());
                    triples.pop_back();
                }
                buckets.emplace_back(bucket);
            }
            if (mPartyIdx == 0) {
                std::cout << "3" << std::endl;
            }

            // 4) Check Buckets
            std::vector<std::array<sf64Matrix<D>, 3>> out;
            for (int i = 0; i < N; ++i) {
                for (int j = 1; j < B; ++j) {
                    if(!this->verifyTripleUsingAnother(comm, buckets.at(i).at(0), buckets.at(i).at(j))) {
                        return std::vector<std::array<sf64Matrix<D>, 3>>{};
                    }
                }
                out.emplace_back(buckets.at(i).at(0));
            }
            if (mPartyIdx == 0) {
                std::cout << "4" << std::endl;
            }

            return out;
        }

        u64 mPartyIdx;
        Sh3Encryptor mEncryptor;
        Sh3Evaluator mEvaluator;
        Sh3Runtime mRuntime;
    private:
        std::vector<i64> coin(CommPkg &comm, int s);

        Sh3Task coin(Sh3Task dep, int s, std::vector<i64> &dest);
    };
}
