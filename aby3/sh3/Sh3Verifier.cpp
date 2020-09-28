#include "Sh3Types.h"
#include "Sh3Verifier.h"

namespace aby3{
    bool Sh3Verifier::verifyTripleWithOpening(CommPkg& comm, const std::array<si64, 3>& abc){
        auto a = mEncryptor.revealAll(comm, abc[0]);
        auto b = mEncryptor.revealAll(comm, abc[1]);
        auto c = mEncryptor.revealAll(comm, abc[2]);
        if (c == a*b) {
            return true;
        }
        comm.mNext.cancel();
        comm.mPrev.cancel();
        return false;
    }

    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<si64, 3> &abc, bool &dest) {
        return dep.then([this, &abc, &dest](CommPkg& comm, Sh3Task& self){
            i64 a = this->mEncryptor.revealAll(comm, abc[0]);
            i64 b = this->mEncryptor.revealAll(comm, abc[1]);
            i64 c = this->mEncryptor.revealAll(comm, abc[2]);
            dest = (c == a*b);
        });
    }

    bool Sh3Verifier::verifyTripleWithOpening(CommPkg& comm, const std::array<sb64, 3>& abc) {
        auto a = mEncryptor.revealAll(comm, abc[0]);
        auto b = mEncryptor.revealAll(comm, abc[1]);
        auto c = mEncryptor.revealAll(comm, abc[2]);
        if (c == a * b) {
            return true;
        }
        comm.mNext.cancel();
        comm.mPrev.cancel();
        return false;
    }

    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<sb64, 3> &abc, bool &dest) {
        return dep.then([this, &abc, &dest](CommPkg& comm, Sh3Task& self){
            i64 a = this->mEncryptor.revealAll(comm, abc[0]);
            i64 b = this->mEncryptor.revealAll(comm, abc[1]);
            i64 c = this->mEncryptor.revealAll(comm, abc[2]);
            dest = (c == a*b);
        });
    }
}