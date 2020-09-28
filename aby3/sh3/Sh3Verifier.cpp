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
            auto a = this->mEncryptor.revealAll(comm, abc[0]);
            auto b = this->mEncryptor.revealAll(comm, abc[1]);
            auto c = this->mEncryptor.revealAll(comm, abc[2]);
            dest = (c == a*b);
        });
    }

    template<Decimal D>
    bool Sh3Verifier::verifyTripleWithOpening(CommPkg& comm, const std::array<sf64<D>, 3>& abc){
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

    template<Decimal D>
    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<sf64<D>, 3> &abc, bool &dest) {
        return dep.then([this, &abc, &dest](CommPkg& comm, Sh3Task& self){
            auto a = this->mEncryptor.revealAll(comm, abc[0]);
            auto b = this->mEncryptor.revealAll(comm, abc[1]);
            auto c = this->mEncryptor.revealAll(comm, abc[2]);
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
            auto a = this->mEncryptor.revealAll(comm, abc[0]);
            auto b = this->mEncryptor.revealAll(comm, abc[1]);
            auto c = this->mEncryptor.revealAll(comm, abc[2]);
            dest = (c == a*b);
        });
    }

    bool Sh3Verifier::verifyTripleWithOpening(CommPkg& comm, const std::array<si64Matrix, 3>& abc){
        i64Matrix a = i64Matrix(abc[0].rows(), abc[0].cols());
        i64Matrix b = i64Matrix(abc[1].rows(), abc[1].cols());
        i64Matrix c = i64Matrix(abc[2].rows(), abc[2].cols());
        mEncryptor.revealAll(comm, abc[0], a);
        mEncryptor.revealAll(comm, abc[1], b);
        mEncryptor.revealAll(comm, abc[2], c);
        if (c == a*b) {
            return true;
        }
        comm.mNext.cancel();
        comm.mPrev.cancel();
        return false;
    }

    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<si64Matrix, 3>& abc, bool& dest){
        return dep.then([this, &abc, &dest](CommPkg& comm, Sh3Task& self) {
            i64Matrix a = i64Matrix(abc[0].rows(), abc[0].cols());
            i64Matrix b = i64Matrix(abc[1].rows(), abc[1].cols());
            i64Matrix c = i64Matrix(abc[2].rows(), abc[2].cols());
            mEncryptor.revealAll(comm, abc[0], a);
            mEncryptor.revealAll(comm, abc[1], b);
            mEncryptor.revealAll(comm, abc[2], c);
            dest = (c == a*b);
        });
    }

    template<Decimal D>
    bool Sh3Verifier::verifyTripleWithOpening(CommPkg& comm, const std::array<sf64Matrix<D>, 3>& abc) {
        fpMatrix<i64, D> a = fpMatrix<i64, D>(abc[0].rows(), abc[0].cols());
        fpMatrix<i64, D> b = fpMatrix<i64, D>(abc[1].rows(), abc[1].cols());
        fpMatrix<i64, D> c = fpMatrix<i64, D>(abc[2].rows(), abc[2].cols());
        mEncryptor.revealAll(comm, abc[0], a);
        mEncryptor.revealAll(comm, abc[1], b);
        mEncryptor.revealAll(comm, abc[2], c);
        if (c == a*b) {
            return true;
        }
        comm.mNext.cancel();
        comm.mPrev.cancel();
        return false;
    }

    template<Decimal D>
    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<sf64Matrix<D>, 3>& abc, bool& dest) {
        return dep.then([this, &abc, &dest](CommPkg& comm, Sh3Task& self) {
            fpMatrix<i64, D> a = fpMatrix<i64, D>(abc[0].rows(), abc[0].cols());
            fpMatrix<i64, D> b = fpMatrix<i64, D>(abc[1].rows(), abc[1].cols());
            fpMatrix<i64, D> c = fpMatrix<i64, D>(abc[2].rows(), abc[2].cols());
            mEncryptor.revealAll(comm, abc[0], a);
            mEncryptor.revealAll(comm, abc[1], b);
            mEncryptor.revealAll(comm, abc[2], c);
            dest = (c == a*b);
        });
    }

    bool Sh3Verifier::verifyTripleWithOpening(CommPkg& comm, const std::array<sbMatrix, 3>& abc) {
        i64Matrix a = i64Matrix(abc[0].rows(), abc[0].i64Cols());
        i64Matrix b = i64Matrix(abc[1].rows(), abc[1].i64Cols());
        i64Matrix c = i64Matrix(abc[2].rows(), abc[2].i64Cols());
        mEncryptor.revealAll(comm, abc[0], a);
        mEncryptor.revealAll(comm, abc[1], b);
        mEncryptor.revealAll(comm, abc[2], c);
        if (c == a*b) {
            return true;
        }
        comm.mNext.cancel();
        comm.mPrev.cancel();
        return false;
    }

    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<sbMatrix, 3>& abc, bool& dest) {
        return dep.then([this, &abc, &dest](CommPkg& comm, Sh3Task& self) {
            i64Matrix a = i64Matrix(abc[0].rows(), abc[0].i64Cols());
            i64Matrix b = i64Matrix(abc[1].rows(), abc[1].i64Cols());
            i64Matrix c = i64Matrix(abc[2].rows(), abc[2].i64Cols());
            mEncryptor.revealAll(comm, abc[0], a);
            mEncryptor.revealAll(comm, abc[1], b);
            mEncryptor.revealAll(comm, abc[2], c);
            dest = (c == a*b);
        });
    }
}