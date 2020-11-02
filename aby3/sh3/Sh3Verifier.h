#pragma once

#include "Sh3Encryptor.h"

namespace aby3 {
    class Sh3Verifier {
    public:
        void init(u64 partyIdx, block prevSeed, block nextSeed, u64 buffSize = 256) {
            mEncryptor.init(partyIdx, prevSeed, nextSeed, buffSize);
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
        bool verifyTripleWithOpening(CommPkg &comm, const std::array<sf64Matrix<D>, 3> &abc);

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
                                      const std::array<sf64Matrix<D>, 3> &abc);

        bool
        verifyTripleUsingAnother(CommPkg &comm, const std::array<sbMatrix, 3> &xyz, const std::array<sbMatrix, 3> &abc);


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
        verifyTripleUsingAnother(Sh3Task dep, const std::array<sbMatrix, 3> &xyz, const std::array<sbMatrix, 3> &abc,
                                 bool &dest);

        template<typename T>
        std::vector<T> perm(CommPkg &comm, std::vector<T> &d);

        u64 mPartyIdx;
        Sh3Encryptor mEncryptor;
    private:
        std::vector<i64> coin(CommPkg &comm, int s);

        Sh3Task coin(Sh3Task dep, int s, std::vector<i64> &dest);
    };
}
