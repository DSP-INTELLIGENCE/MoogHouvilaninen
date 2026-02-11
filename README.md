# MoogHouvilaninen
```c++
#ifndef GAMMA_MOOG_HUOVILAINEN_H_INC
#define GAMMA_MOOG_HUOVILAINEN_H_INC

#include "Gamma/Domain.h"
#include "Gamma/Types.h"
#include "Gamma/scl.h"
#include <cmath>
#include <cstdio>

namespace gam {

template<class Tv=gam::real, class Tp=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class MoogHuovilainen : public Td {
public:

    MoogHuovilainen(Tp cutoff=Tp(1000), Tp resonance=Tp(0.1))
    :   mCutoff(cutoff),
        mResonance(resonance),
        mThermal(Tp(0.000025)),
        mTune(Tp(0)),
        mAcr(Tp(1)),
        mResQuad(Tp(0))
    {
        reset();
        // DO NOT compute here
        set(1000.0f,1.0f);
        print();
    }

    // ------------------------------------------------------------
    // Public API
    // ------------------------------------------------------------

    void cutoff(Tp v){
        mCutoff = v;
        updateCoef();
    }

    void resonance(Tp v){
        mResonance = v;
        updateCoef();
    }

    void set(Tp c, Tp r){
        mCutoff = c;
        mResonance = r;
        updateCoef();
    }

    Tp cutoff() const { return mCutoff; }
    Tp resonance() const { return mResonance; }

    void thermal(Tp v){
        mThermal = (v > Tp(0)) ? v : Tp(0.000025);
        updateCoef();
    }

    void reset(){
        for(int i=0;i<4;++i) mStage[i] = Tv(0);
        for(int i=0;i<3;++i) mStageTanh[i] = Tv(0);
        for(int i=0;i<6;++i) mDelay[i] = Tv(0);
    }

    void onDomainChange(double){
        updateCoef();
    }

    // ------------------------------------------------------------
    // Processing
    // ------------------------------------------------------------

    Tv operator()(Tv in){

        Tv x = in;

        // 2x oversampling
        for(int j=0;j<2;++j){

            Tv inp = x - Tv(mResQuad) * mDelay[5];

            Tv t0 = Tv(std::tanh(Tp(inp) * mThermal));
            mDelay[0] = mDelay[0] + Tv(mTune) * (t0 - mStageTanh[0]);
            mStage[0] = mDelay[0];

            for(int k=1;k<4;++k){

                Tv s = mStage[k-1];
                Tv t = Tv(std::tanh(Tp(s) * mThermal));
                mStageTanh[k-1] = t;

                Tv sub = (k!=3)
                        ? mStageTanh[k]
                        : Tv(std::tanh(Tp(mDelay[k]) * mThermal));

                Tv y = mDelay[k] + Tv(mTune) * (t - sub);

                mStage[k] = y;
                mDelay[k] = y;
            }

            mDelay[5] = Tv(0.5) * (mStage[3] + mDelay[4]);
            mDelay[4] = mStage[3];
        }

        return mDelay[5];
    }

    // ------------------------------------------------------------
    // Debug
    // ------------------------------------------------------------

    void print() const {
        std::printf("SR:        %g\n", double(Td::spu()));
        std::printf("Cut:       %g\n", double(mCutoff));
        std::printf("Res:       %g\n", double(mResonance));
        std::printf("Tune:      %g\n", double(mTune));
        std::printf("ACR:       %g\n", double(mAcr));
        std::printf("ResQuad:   %g\n", double(mResQuad));
        std::printf("Thermal:   %g\n", double(mThermal));
    }

private:

    // state
    Tv mStage[4];
    Tv mStageTanh[3];
    Tv mDelay[6];

    // parameters
    Tp mCutoff;
    Tp mResonance;
    Tp mThermal;

    // coefficients
    Tp mTune;
    Tp mAcr;
    Tp mResQuad;

    static inline Tp pi(){
        return Tp(3.1415926535897932384626433832795);
    }

    inline Tp clampCutoff(Tp c, Tp sr){
        return std::clamp(c, Tp(1e-6), Tp(0.49)*sr);
    }

    inline Tp clampRes(Tp r){
        return std::clamp(r, Tp(0), Tp(4));
    }

    void updateCoef(){

        Tp sr = Tp(Td::spu());
        
        // clamp ONLY into locals
        Tp cutoff = clampCutoff(mCutoff, sr);
        Tp resonance = clampRes(mResonance);

            
        Tp fc = cutoff / sr;
        Tp f  = fc * Tp(0.5); // oversampled

        Tp fc2 = fc * fc;
        Tp fc3 = fc2 * fc;

        Tp fcr = Tp(1.8730)*fc3
               + Tp(0.4955)*fc2
               - Tp(0.6490)*fc
               + Tp(0.9988);

        mAcr = Tp(-3.9364)*fc2
             + Tp(1.8409)*fc
             + Tp(0.9968);

        Tp expo = -Tp(2) * pi() * f * fcr;

        mTune = (Tp(1) - Tp(std::exp(double(expo)))) / mThermal;

        mResQuad = Tp(4) * resonance * mAcr;
    }
};

} // namespace gam

#endif

```
