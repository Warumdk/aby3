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

    bool Sh3Verifier::compareView(CommPkg &comm, i64& x) {
        // TODO: Move compare view to hash function optimization
        comm.mNext.asyncSendCopy(x);

        i64 xPrev;
        comm.mPrev.recv(xPrev);

        if (x != xPrev) {
            comm.mPrev.cancel();
            comm.mNext.cancel();
            return false;
        }
        return true;
    }

    Sh3Task Sh3Verifier::compareView(Sh3Task dep, i64& x, bool& dest) {
        // TODO: Move compare view to hash function optimization
        return dep.then([&x, &dest](CommPkg& comm, Sh3Task& self) {
            comm.mNext.asyncSendCopy(x);

            i64 xPrev;
            comm.mPrev.recv(xPrev);
            dest = (x == xPrev);
        });
    }

    bool Sh3Verifier::verifyTripleUsingAnother(CommPkg& comm, const std::array<si64, 3>& xyz, const std::array<si64, 3>& abc) {
        si64 sRho = xyz[0] - abc[0];
        si64 sSigma = xyz[1] - abc[1];
        i64 rho = mEncryptor.revealAll(comm, sRho);
        i64 sigma = mEncryptor.revealAll(comm, sSigma);

        if (!this->compareView(comm, rho)) {
            return  false;
        }

        if (!this->compareView(comm, sigma)) {
            return  false;
        }

        si64 sDelta = xyz[2] - abc[2] - (sigma * abc[0]) - (rho * abc[1]) - sigma * rho;

        i64 delta = mEncryptor.revealAll(comm, sDelta);

        if (delta != 0) {
            comm.mPrev.cancel();
            comm.mNext.cancel();
            return false;
        }

        return this->compareView(comm, delta);
    }

    bool Sh3Verifier::verifyTripleUsingAnother(CommPkg& comm, const std::array<sb64, 3>& xyz, const std::array<sb64, 3>& abc) {
        sb64 sRho = xyz[0] ^ abc[0];
        sb64 sSigma = xyz[1] ^ abc[1];
        i64 rho = mEncryptor.revealAll(comm, sRho);
        i64 sigma = mEncryptor.revealAll(comm, sSigma);

        if (!this->compareView(comm, rho)) {
            return  false;
        }

        if (!this->compareView(comm, sigma)) {
            return  false;
        }

        sb64 sDelta = xyz[2] ^ abc[2] ^ (sigma * abc[0]) ^ (rho * abc[1]) ^ sigma * rho;

        i64 delta = mEncryptor.revealAll(comm, sDelta);

        if (delta != 0) {
            comm.mPrev.cancel();
            comm.mNext.cancel();
            return false;
        }

        return this->compareView(comm, delta);
    }

    Sh3Task Sh3Verifier::verifyTripleUsingAnother(Sh3Task dep, const std::array<si64, 3>& xyz, const std::array<si64, 3>& abc, bool& dest) {
        return dep.then([this, &xyz, &abc, &dest](CommPkg& comm, Sh3Task& self){
            si64 sRho = xyz[0] - abc[0];
            si64 sSigma = xyz[1] - abc[1];
            i64 rho = this->mEncryptor.revealAll(comm, sRho);
            i64 sigma = this->mEncryptor.revealAll(comm, sSigma);

            if (!this->compareView(comm, rho)) {
                dest = false;
            }

            if (!this->compareView(comm, sigma)) {
                dest = false;
            }

            si64 sDelta = xyz[2] - abc[2] - (sigma * abc[0]) - (rho * abc[1]) - sigma * rho;

            i64 delta = mEncryptor.revealAll(comm, sDelta);

            if (delta != 0) {
                comm.mPrev.cancel();
                comm.mNext.cancel();
                return false;
            }

            dest = this->compareView(comm, delta);
        });
    }
}