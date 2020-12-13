#include "Sh3Types.h"
#include "Sh3Verifier.h"

namespace aby3 {
    bool Sh3Verifier::verifyTripleWithOpening(CommPkg &comm, const std::array<si64, 3> &abc) {
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

    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<si64, 3> &abc, bool &dest) {
        return dep.then([this, &abc, &dest](CommPkg &comm, Sh3Task &self) {
            auto a = this->mEncryptor.revealAll(comm, abc[0]);
            auto b = this->mEncryptor.revealAll(comm, abc[1]);
            auto c = this->mEncryptor.revealAll(comm, abc[2]);
            dest = (c == a * b);
        });
    }

    template<Decimal D>
    bool Sh3Verifier::verifyTripleWithOpening(CommPkg &comm, const std::array<sf64<D>, 3> &abc) {
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

    template<Decimal D>
    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<sf64<D>, 3> &abc, bool &dest) {
        return dep.then([this, &abc, &dest](CommPkg &comm, Sh3Task &self) {
            auto a = this->mEncryptor.revealAll(comm, abc[0]);
            auto b = this->mEncryptor.revealAll(comm, abc[1]);
            auto c = this->mEncryptor.revealAll(comm, abc[2]);
            dest = (c == a * b);
        });
    }

    bool Sh3Verifier::verifyTripleWithOpening(CommPkg &comm, const std::array<sb64, 3> &abc) {
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
        return dep.then([this, &abc, &dest](CommPkg &comm, Sh3Task &self) {
            auto a = this->mEncryptor.revealAll(comm, abc[0]);
            auto b = this->mEncryptor.revealAll(comm, abc[1]);
            auto c = this->mEncryptor.revealAll(comm, abc[2]);
            dest = (c == a * b);
        });
    }

    bool Sh3Verifier::verifyTripleWithOpening(CommPkg &comm, const std::array<si64Matrix, 3> &abc) {
        i64Matrix a = i64Matrix(abc[0].rows(), abc[0].cols());
        i64Matrix b = i64Matrix(abc[1].rows(), abc[1].cols());
        i64Matrix c = i64Matrix(abc[2].rows(), abc[2].cols());
        mEncryptor.revealAll(comm, abc[0], a);
        mEncryptor.revealAll(comm, abc[1], b);
        mEncryptor.revealAll(comm, abc[2], c);
        if (c == a * b) {
            return true;
        }
        comm.mNext.cancel();
        comm.mPrev.cancel();
        return false;
    }

    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<si64Matrix, 3> &abc, bool &dest) {
        return dep.then([this, &abc, &dest](CommPkg &comm, Sh3Task &self) {
            i64Matrix a = i64Matrix(abc[0].rows(), abc[0].cols());
            i64Matrix b = i64Matrix(abc[1].rows(), abc[1].cols());
            i64Matrix c = i64Matrix(abc[2].rows(), abc[2].cols());
            mEncryptor.revealAll(comm, abc[0], a);
            mEncryptor.revealAll(comm, abc[1], b);
            mEncryptor.revealAll(comm, abc[2], c);
            dest = (c == a * b);
        });
    }



    template<Decimal D>
    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<sf64Matrix<D>, 3> &abc, bool &dest) {
        return dep.then([this, &abc, &dest](CommPkg &comm, Sh3Task &self) {
            fpMatrix<i64, D> a = fpMatrix<i64, D>(abc[0].rows(), abc[0].cols());
            fpMatrix<i64, D> b = fpMatrix<i64, D>(abc[1].rows(), abc[1].cols());
            fpMatrix<i64, D> c = fpMatrix<i64, D>(abc[2].rows(), abc[2].cols());
            mEncryptor.revealAll(comm, abc[0], a);
            mEncryptor.revealAll(comm, abc[1], b);
            mEncryptor.revealAll(comm, abc[2], c);
            dest = (c == a * b);
        });
    }

    bool Sh3Verifier::verifyTripleWithOpening(CommPkg &comm, const std::array<sbMatrix, 3> &abc) {
        i64Matrix a = i64Matrix(abc[0].rows(), abc[0].i64Cols());
        i64Matrix b = i64Matrix(abc[1].rows(), abc[1].i64Cols());
        i64Matrix c = i64Matrix(abc[2].rows(), abc[2].i64Cols());
        mEncryptor.revealAll(comm, abc[0], a);
        mEncryptor.revealAll(comm, abc[1], b);
        mEncryptor.revealAll(comm, abc[2], c);
        if (c == a * b) {
            return true;
        }
        comm.mNext.cancel();
        comm.mPrev.cancel();
        return false;
    }

    Sh3Task Sh3Verifier::verifyTripleWithOpening(Sh3Task dep, const std::array<sbMatrix, 3> &abc, bool &dest) {
        return dep.then([this, &abc, &dest](CommPkg &comm, Sh3Task &self) {
            i64Matrix a = i64Matrix(abc[0].rows(), abc[0].i64Cols());
            i64Matrix b = i64Matrix(abc[1].rows(), abc[1].i64Cols());
            i64Matrix c = i64Matrix(abc[2].rows(), abc[2].i64Cols());
            mEncryptor.revealAll(comm, abc[0], a);
            mEncryptor.revealAll(comm, abc[1], b);
            mEncryptor.revealAll(comm, abc[2], c);
            dest = (c == a * b);
        });
    }

    bool Sh3Verifier::compareView(CommPkg &comm, i64 &x) {
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

    bool Sh3Verifier::compareView(CommPkg &comm, i64Matrix &x) {
        comm.mNext.asyncSendCopy(x.data(), x.size());

        i64Matrix xPrev = i64Matrix(x.rows(), x.cols());
        comm.mPrev.recv(xPrev.data(), xPrev.size());

        if (x != xPrev) {
            comm.mPrev.cancel();
            comm.mNext.cancel();
            return false;
        }
        return true;
    }

    Sh3Task Sh3Verifier::compareView(Sh3Task dep, i64 &x, bool &dest) {
        // TODO: Move compare view to hash function optimization
        return dep.then([&x, &dest](CommPkg &comm, Sh3Task &self) {
            comm.mNext.asyncSendCopy(x);

            i64 xPrev;
            comm.mPrev.recv(xPrev);
            dest = (x == xPrev);
        });
    }

    bool Sh3Verifier::verifyTripleUsingAnother(CommPkg &comm, const std::array<si64, 3> &xyz,
                                               const std::array<si64, 3> &abc) {
        si64 sRho = xyz[0] - abc[0];
        si64 sSigma = xyz[1] - abc[1];
        i64 rho = mEncryptor.revealAll(comm, sRho);
        i64 sigma = mEncryptor.revealAll(comm, sSigma);

        if (!this->compareView(comm, rho)) {
            return false;
        }

        if (!this->compareView(comm, sigma)) {
            return false;
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

    bool Sh3Verifier::verifyTripleUsingAnother(CommPkg &comm, const std::array<sb64, 3> &xyz,
                                               const std::array<sb64, 3> &abc) {
        sb64 sRho = xyz[0] ^abc[0];
        sb64 sSigma = xyz[1] ^abc[1];
        i64 rho = mEncryptor.revealAll(comm, sRho);
        i64 sigma = mEncryptor.revealAll(comm, sSigma);

        if (!this->compareView(comm, rho)) {
            return false;
        }

        if (!this->compareView(comm, sigma)) {
            return false;
        }

        sb64 sDelta = xyz[2] ^abc[2] ^(sigma * abc[0]) ^(rho * abc[1]) ^sigma * rho;

        i64 delta = mEncryptor.revealAll(comm, sDelta);

        if (delta != 0) {
            comm.mPrev.cancel();
            comm.mNext.cancel();
            return false;
        }

        return this->compareView(comm, delta);
    }


    template<Decimal D>
    bool Sh3Verifier::verifyTripleUsingAnother(CommPkg &comm, const std::array<sf64<D>, 3> &xyz,
                                               const std::array<sf64<D>, 3> &abc) {
        si64 sRho = xyz[0] - abc[0];
        si64 sSigma = xyz[1] - abc[1];
        i64 rho = mEncryptor.revealAll(comm, sRho);
        i64 sigma = mEncryptor.revealAll(comm, sSigma);

        if (!this->compareView(comm, rho)) {
            return false;
        }

        if (!this->compareView(comm, sigma)) {
            return false;
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

    bool Sh3Verifier::verifyTripleUsingAnother(CommPkg &comm, const std::array<si64Matrix, 3> &xyz,
                                               const std::array<si64Matrix, 3> &abc) {
        si64Matrix sRho = xyz[0] - abc[0];
        si64Matrix sSigma = xyz[1] - abc[1];
        i64Matrix rho = i64Matrix(abc[0].rows(), abc[0].cols());
        i64Matrix sigma = i64Matrix (abc[1].rows(), abc[1].cols());
        mEncryptor.revealAll(comm, sRho, rho);
        mEncryptor.revealAll(comm, sSigma, sigma);

        if (!this->compareView(comm, rho)) {
            return false;
        }

        if (!this->compareView(comm, sigma)) {
            return false;
        }

        si64Matrix sDelta = (xyz[2] - abc[2]) - (abc[0] * sigma) - (rho * abc[1]) - (rho * sigma);

        i64Matrix delta = i64Matrix(sDelta.rows(), sDelta.cols());
        mEncryptor.revealAll(comm, sDelta, delta);

        if (delta != i64Matrix::Zero(delta.rows(), delta.cols())) {
            if (mPartyIdx == 0) {
                std::stringstream msg;
                msg << mPartyIdx << std::endl << rho << std::endl << sigma << std::endl << delta << std::endl;
                std::cout << msg.str();
            }
            comm.mPrev.cancel();
            comm.mNext.cancel();
            return false;
        }

        return this->compareView(comm, delta);
    }



    bool Sh3Verifier::verifyTripleUsingAnother(CommPkg &comm, const std::array<sbMatrix, 3> &xyz,
                                               const std::array<si64Matrix, 3> &abc) {
        sbMatrix sRho = xyz[0] ^ abc[0];
        sbMatrix sSigma = xyz[1] ^ abc[1];
        i64Matrix rho = i64Matrix(abc[0].rows(), abc[0].cols());
        i64Matrix sigma = i64Matrix (abc[1].rows(), abc[1].cols());
        mEncryptor.revealAll(comm, sRho, rho);
        mEncryptor.revealAll(comm, sSigma, sigma);

        if (!this->compareView(comm, rho)) {
            return false;
        }

        if (!this->compareView(comm, sigma)) {
            return false;
        }

        sbMatrix sDelta = xyz[2] ^ abc[2] - (sigma * abc[0]) - (rho * abc[1]) - sigma * rho;
        i64Matrix delta = i64Matrix(sDelta.rows(), sDelta.i64Cols());
        mEncryptor.revealAll(comm, sDelta, delta);

        if (delta != i64Matrix(delta.rows(), delta.cols())) {
            comm.mPrev.cancel();
            comm.mNext.cancel();
            return false;
        }

        return this->compareView(comm, delta);
    }

    Sh3Task
    Sh3Verifier::verifyTripleUsingAnother(Sh3Task dep, const std::array<si64, 3> &xyz, const std::array<si64, 3> &abc,
                                          bool &dest) {
        return dep.then([this, &xyz, &abc, &dest](CommPkg &comm, Sh3Task &self) {
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
                dest = false;
            }

            dest = this->compareView(comm, delta);
        });
    }

    Sh3Task
    Sh3Verifier::verifyTripleUsingAnother(Sh3Task dep, const std::array<sb64, 3> &xyz, const std::array<sb64, 3> &abc,
                                          bool &dest) {
        return dep.then([this, &xyz, &abc, &dest](CommPkg &comm, Sh3Task self) {
            sb64 sRho = xyz[0] ^abc[0];
            sb64 sSigma = xyz[1] ^abc[1];
            i64 rho = mEncryptor.revealAll(comm, sRho);
            i64 sigma = mEncryptor.revealAll(comm, sSigma);

            if (!this->compareView(comm, rho)) {
                dest = false;
            }

            if (!this->compareView(comm, sigma)) {
                dest = false;
            }

            sb64 sDelta = xyz[2] ^abc[2] ^(sigma * abc[0]) ^(rho * abc[1]) ^sigma * rho;

            i64 delta = mEncryptor.revealAll(comm, sDelta);

            if (delta != 0) {
                comm.mPrev.cancel();
                comm.mNext.cancel();
                dest = false;
            }

            dest = this->compareView(comm, delta);
        });
    }

    template<Decimal D>
    Sh3Task Sh3Verifier::verifyTripleUsingAnother(Sh3Task dep, const std::array<sf64<D>, 3> &xyz,
                                                  const std::array<sf64<D>, 3> &abc, bool &dest) {
        return this->verifyTripleUsingAnother(dep, xyz, abc, dest);
    }

    Sh3Task Sh3Verifier::verifyTripleUsingAnother(Sh3Task dep, const std::array<si64Matrix, 3> &xyz,
                                                  const std::array<si64Matrix, 3> &abc, bool &dest) {
        return dep.then([this, &xyz, &abc, &dest](CommPkg &comm, Sh3Task self){
            si64Matrix sRho = xyz[0] - abc[0];
            si64Matrix sSigma = xyz[1] - abc[1];
            i64Matrix rho = i64Matrix(abc[0].rows(), abc[0].cols());
            i64Matrix sigma = i64Matrix (abc[1].rows(), abc[1].cols());
            mEncryptor.revealAll(comm, sRho, rho);
            mEncryptor.revealAll(comm, sSigma, sigma);

            if (!this->compareView(comm, rho)) {
                dest = false;
            }

            if (!this->compareView(comm, sigma)) {
                dest = false;
            }

            si64Matrix sDelta = xyz[2] - abc[2] - (sigma * abc[0]) - (rho * abc[1]) - sigma * rho;
            i64Matrix delta = i64Matrix(sDelta.rows(), sDelta.cols());
            mEncryptor.revealAll(comm, sDelta, delta);

            if (delta != i64Matrix(sDelta.rows(), sDelta.cols())) {
                comm.mPrev.cancel();
                comm.mNext.cancel();
                dest = false;
            }

            dest = this->compareView(comm, delta);
        });
    }

    template<Decimal D>
    Sh3Task Sh3Verifier::verifyTripleUsingAnother(Sh3Task dep, const std::array<sf64Matrix<D>, 3> &xyz,
                                                  const std::array<sf64Matrix<D>, 3> &abc, bool &dest) {
        return this->verifyTripleUsingAnother(dep, xyz, abc, dest);
    }

    Sh3Task Sh3Verifier::verifyTripleUsingAnother(Sh3Task dep, const std::array<sbMatrix, 3> &xyz,
                                                  const std::array<si64Matrix, 3> &abc, bool &dest) {
        return dep.then([this, &xyz, &abc, &dest](CommPkg &comm, Sh3Task self){
            sbMatrix sRho = xyz[0] ^ abc[0];
            sbMatrix sSigma = xyz[1] ^ abc[1];
            i64Matrix rho = i64Matrix(abc[0].rows(), abc[0].cols());
            i64Matrix sigma = i64Matrix (abc[1].rows(), abc[1].cols());
            mEncryptor.revealAll(comm, sRho, rho);
            mEncryptor.revealAll(comm, sSigma, sigma);

            if (!this->compareView(comm, rho)) {
                dest = false;
            }

            if (!this->compareView(comm, sigma)) {
                dest = false;
            }

            sbMatrix sDelta = xyz[2] ^ abc[2] - (sigma * abc[0]) - (rho * abc[1]) - sigma * rho;
            i64Matrix delta = i64Matrix(sDelta.rows(), sDelta.i64Cols());
            mEncryptor.revealAll(comm, sDelta, delta);

            if (delta != i64Matrix(delta.rows(), delta.cols())) {
                comm.mPrev.cancel();
                comm.mNext.cancel();
                dest = false;
            }

            dest = this->compareView(comm, delta);
        });
    }

    std::vector<i64> Sh3Verifier::coin(CommPkg &comm, int s) {
        std::vector<si64> shares;
        for (int i = 0; i < s; ++i) {
            shares.emplace_back(this->mEncryptor.mShareGen.getRandIntShare());
        }
        std::vector<i64> out;
        for (si64 share : shares) {
            out.emplace_back(this->mEncryptor.revealAll(comm, share));
        }
        for (i64 val : out) {
            bool cv = this->compareView(comm, val);
            if (!cv) {
                comm.mNext.close();
                comm.mPrev.close();
                return std::vector<i64>{};
            }
        }
        return out;
    }

    Sh3Task Sh3Verifier::coin(Sh3Task dep, int s, std::vector<i64> &dest) {
        return dep.then([this, s, &dest](CommPkg &comm, Sh3Task self){
            std::vector<si64> shares;
            for (int i = 0; i < s; ++i) {
                shares.emplace_back(this->mEncryptor.mShareGen.getRandIntShare());
            }
            std::vector<i64> out;
            for (si64 share : shares) {
                out.emplace_back(this->mEncryptor.revealAll(comm, share));
            }
            for (i64 val : out) {
                bool cv = this->compareView(comm, val);
                if (!cv) {
                    comm.mNext.close();
                    comm.mPrev.close();
                    dest = std::vector<i64>{};
                }
            }
            dest = out;
        });
    }


    template<typename T>
    Sh3Task Sh3Verifier::perm(Sh3Task dep, std::vector<T> &d, std::vector<T> &dest) {
        return dep.then([this, &d, &dest](CommPkg &comm, Sh3Task self){
            std::vector<i64> indices = this->coin(comm, d.size());
            for (int j = d.size()-1; j >= 0; --j) {
                std::swap(d[j], d[indices[j] % d.size()]);
            }
            dest = d;
        });
    }

    std::vector<std::array<si64, 3>>
    Sh3Verifier::generateTriples(CommPkg &comm, int N, int B, int C){
        // 1) Generate Random Sharings
        int M = N*B+C;
        si64 shares[2*M];
        for (int i = 0; i < 2*M; ++i) {
            shares[i] = mEncryptor.mShareGen.getRandIntShare();
        }

        // 2) Generate Multiplication Triples
        std::vector<std::array<si64, 3>> triples;
        for (int i = 0; i < (2*M)-1; i+=2) {
            si64 res{};
            mEvaluator.asyncMul(mRuntime, shares[i], shares[i+1], res).get();
            std::array<si64, 3> triple = {shares[i], shares[i+1], res};
            triples.emplace_back(triple);
        }

        // 3) Cut-and-bucket
        triples = this->perm(comm, triples);
        for (int i = 0; i < C; ++i) {
            if (!this->verifyTripleWithOpening(comm, triples.back())) {
                return std::vector<std::array<si64, 3>>{};
            }
            triples.pop_back();
        }

        std::vector<std::vector<std::array<si64, 3>>> buckets;
        for (int i = 0; i < N; ++i) {
            std::vector<std::array<si64, 3>> bucket;
            for (int j = 0; j < B; ++j) {
                bucket.emplace_back(triples.back());
                triples.pop_back();
            }
            buckets.emplace_back(bucket);
        }

        // 4) Check Buckets
        std::vector<std::array<si64, 3>> out;
        for (int i = 0; i < N; ++i) {
            for (int j = 1; j < B; ++j) {
                if(!this->verifyTripleUsingAnother(comm, buckets.at(i).at(0), buckets.at(i).at(j))) {
                    return std::vector<std::array<si64, 3>>{};
                }
            }
            out.emplace_back(buckets.at(i).at(0));
        }

        return out;
    }

    std::vector<std::array<si64Matrix, 3>> Sh3Verifier::generateTriples(CommPkg &comm, int N, int B, int C, int rowsA, int colsA, int rowsB, int colsB) {
        // 1) Generate Random Sharings
        int M = N*B+C;
        std::vector<std::array<si64Matrix, 3>> triples;
        for (int i = 0; i < M; ++i) {
            si64Matrix aMatrix(rowsA, colsA);
            si64Matrix bMatrix(rowsB, colsB);
            si64Matrix res;

            this->mEncryptor.rand(aMatrix);
            this->mEncryptor.rand(bMatrix);
            mEvaluator.asyncMul(mRuntime, aMatrix, bMatrix, res).get();
            std::array<si64Matrix, 3> triple = {aMatrix, bMatrix, res};
            triples.emplace_back(triple);
        }
        if (mPartyIdx == 0) {
            std::cout << "1" << std::endl;
        }

        // 3) Cut-and-bucket
        triples = this->perm(comm, triples);
        for (int i = 0; i < C; ++i) {
            if (!this->verifyTripleWithOpening(comm, triples.back())) {
                return std::vector<std::array<si64Matrix, 3>>{};
            }
            triples.pop_back();
        }
        if (mPartyIdx == 0) {
            std::cout << "2" << std::endl;
        }
        std::vector<std::vector<std::array<si64Matrix, 3>>> buckets;
        for (int i = 0; i < N; ++i) {
            std::vector<std::array<si64Matrix, 3>> bucket;
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
        std::vector<std::array<si64Matrix, 3>> out;
        for (int i = 0; i < N; ++i) {
            for (int j = 1; j < B; ++j) {
                if(!this->verifyTripleUsingAnother(comm, buckets.at(i).at(0), buckets.at(i).at(j))) {
                    return std::vector<std::array<si64Matrix, 3>>{};
                }
            }
            out.emplace_back(buckets.at(i).at(0));
        }
        if (mPartyIdx == 0) {
            std::cout << "4" << std::endl;
        }

        return out;
    }
}