#pragma once
#include "Sh3Encryptor.h"

namespace aby3
{
    class Sh3Verifier {
        public:
            void init(u64 partyIdx, block prevSeed, block nextSeed, u64 buffSize = 256) {mEncryptor.init(partyIdx, prevSeed, nextSeed, buffSize); mPartyIdx = partyIdx; }; 

            bool verifyTripleWithOpening(CommPkg& comm, const std::array<si64, 3>& abc);
            Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<si64, 3>& abc, bool& dest);
            
            template<Decimal D>
            bool verifyTripleWithOpening(CommPkg& comm, const std::array<sf64<D>, 3>& abc);
            template<Decimal D>
            Sh3Task verifyTripleWithOpening(Sh3Task dep,  const std::array<sf64<D>, 3>& abc, bool& dest);

            bool verifyTripleWithOpening(CommPkg& comm, const std::array<sb64, 3>& abc);
            Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<sb64, 3>& abc, bool& dest);

            bool verifyTripleWithOpening(CommPkg& comm, const std::array<si64Matrix, 3>& abc);
            Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<si64Matrix, 3>& abc, bool& dest);

            template<Decimal D>
            bool verifyTripleWithOpening(CommPkg& comm, const std::array<sf64Matrix<D>, 3>& abc);
            template<Decimal D>
            Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<sf64Matrix<D>, 3>& abc, bool& dest) ;

            bool verifyTripleWithOpening(CommPkg& comm, const std::array<sbMatrix, 3>& abc);
            Sh3Task verifyTripleWithOpening(Sh3Task dep, const std::array<sbMatrix, 3>& abc, bool& dest);

            bool compareView(CommPkg& comm, i64& x);

            bool verifyTripleUsingAnother(CommPkg& comm, const std::array<si64, 3>& xyz, const std::array<si64, 3>& abc);
            bool verifyTripleUsingAnother(CommPkg& comm, const std::array<sb64, 3>& xyz, const std::array<sb64, 3>& abc);

        u64 mPartyIdx;
        Sh3Encryptor mEncryptor;
    };
}
