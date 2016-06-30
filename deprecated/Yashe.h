#pragma once

#include <vector>
#include <iostream>
#include <random>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/vector.h>

class Yashe{
public:
    using Ciphertext_t = NTL::ZZ_pX;
    using Plainetxt_t =  NTL::ZZ_pX;

    NTL::ZZ _q;
    NTL::ZZ _t;
    NTL::ZZ _delta;
    NTL::ZZ_pXModulus _PhiModulus;
    NTL::ZZ_pX _Phi;
    int64_t _sigma = 1;
    uint64_t _base = 1;
    uint64_t _q_bitsize = 0;
    uint64_t _base_bitsize = 0;
    uint64_t _logqw = 0;
    int64_t _seed = 0;
    uint64_t _phi_degree = 0;
    uint64_t _nb_slots = 0;
    std::mt19937 gen;
    
    NTL::ZZ_pX _pk;
    NTL::ZZ_pX _sk;
    std::vector<NTL::ZZ_pX> _evk;
 
    // Sample polynomial with Gaussian coefficicents
    NTL::ZZ_pX gaussian_poly(const double center) const {
	    NTL::ZZ_pX P;
	    std::normal_distribution<> dist(center, _sigma);
	
	    P.SetLength(_phi_degree);
	    for(long i = 0 ; i < _phi_degree-1 ; ++i){
		    P[i] = std::lrint(dist(gen));
	    }

      P.normalize();
      return P;
    }
       
    // Represent an integer n in the basis _base
    std::vector<NTL::ZZ> word_decomp_scal(const NTL::ZZ & n) const {
	    NTL::ZZ quo, rem, nn;
	    nn = n;
        std::vector<NTL::ZZ> res(_logqw);
	    for(uint64_t i = 0 ; i < _logqw ; ++i){
		    NTL::DivRem(quo, rem, nn, NTL::ZZ(_base));
            res[i] = rem;
		    nn = quo;
	    }
        return res;
    }
    
    NTL::ZZ power_of_scal(const std::vector<NTL::ZZ> & n) const {
	    NTL::ZZ r(0);
	    for(uint64_t i = 0 ; i < _logqw ; ++i){
		    r += n[i]*(NTL::power(NTL::ZZ(_base), i));
	    }
	    return r;
    }
    
    std::vector<NTL::ZZ_pX> word_decomp_poly(const NTL::ZZ_pX & P) const {
	    std::vector<NTL::ZZ_pX> res(_logqw);
        for(uint64_t i = 0 ; i < _logqw ; ++i){
		    res[i].SetLength(NTL::deg(P)+1);
	    }
        auto dd = NTL::deg(P);
	    for(long i = 0 ; i <= dd ; ++i){
		    auto tmp = word_decomp_scal(rep(P[i]));
		    for(uint64_t j = 0 ; j < _logqw ; ++j){
			    NTL::conv(res[j][i], tmp[j]);
		    }	
	    }
        return res;
    }

    std::vector<NTL::ZZ_pX> power_of_poly(const NTL::ZZ_pX & P) const {
	    NTL::ZZ_pX tmp;
	    tmp.SetLength(deg(P)+1);
        std::vector<NTL::ZZ_pX> res(_logqw);
	    for(uint64_t i = 0 ; i < _logqw ; ++i){
		    NTL::conv(tmp[0], power(NTL::ZZ(_base), i));
		    res[i] = P*tmp; 
	    }
        return res;
    }
    
    NTL::ZZ_pX MulModPhi(const NTL::ZZ_pX & a, const NTL::ZZ_pX & b) const {
        NTL::ZZ_pX r;
        r = a*b;
        auto dd = NTL::deg(r);       
        for(long i  = _phi_degree-1 ; i <= dd ; ++i){
            r[i-(_phi_degree-1)] -= r[i];
             r[i] = 0;
        }
        r.normalize();
        return r;
    }
    
    NTL::ZZ_pX dot_prod(const std::vector<NTL::ZZ_pX> & a, const std::vector<NTL::ZZ_pX> & b) const {
    	NTL::ZZ_pX r(0);
        uint64_t end = a.size();       
    	for(uint64_t i = 0 ; i < end ; ++i){
	    	// r += MulMod(a[i], b[i], _PhiModulus);
            r += MulModPhi(a[i], b[i]);
	    }
	    return r;
    }
 
    uint64_t phi_degree() const {
        return _phi_degree;
    }
     
    std::ostream& printPoly(std::ostream & out, const Ciphertext_t & P) const {
        out << P[0] << " + ";
        for(long i = 1 ; i < NTL::deg(P) ; ++i){
            out << P[i] << "*x^" << i << " + ";
        }
        out << P[NTL::deg(P)] << "*x^" << NTL::deg(P);
        out << std::endl;
        return out;
    }
    
    void number_of_slots(const NTL::ZZX & p, const long verbose = 0) {
        NTL::ZZ_pPush push(_t);
        auto r = NTL::CanZass(NTL::conv<NTL::ZZ_pX>(p), verbose);
        _nb_slots = r.length();
        // std::cout << "Factors degree: ";
        // for(int64_t i = 0 ; i < r.length() ; ++i){
        //     std::cout << NTL::deg(r[i].a) << " ";
        // }
        // std::cout << std::endl;
    }
         
    Yashe() = delete;
    
    Yashe(const NTL::ZZ & q, const NTL::ZZ & t, const uint64_t & d, const uint64_t & sigma, const uint64_t base, const int64_t & seed = 0) 
        : _q(q), _t(t), _sigma(sigma), _base(base), _seed(seed){
            NTL::ZZ_p::init(_q);
            _Phi.SetLength(d+1);
            _Phi[0] = 1;
            _Phi[d] = 1;
            _PhiModulus = NTL::ZZ_pXModulus(_Phi);
            // number_of_slots(NTL::conv<NTL::ZZX>(_Phi));
            // std::cout << "#slots: " << _nb_slots << std::endl;
            _phi_degree = d + 1;
            _delta = q/t;
            _q_bitsize = NTL::NumBits(_q);
            _base_bitsize = NTL::NumBits(_base);
            _logqw = std::ceil((double)_q_bitsize/(double)_base_bitsize);
            if(!_seed){
                std::random_device rd;
                _seed = rd();
            }

            
            gen = std::mt19937(_seed);
        }
        
    Yashe(const Yashe&) = default;
    Yashe(Yashe&&) = default;
    
    // Sample polynomial with uniform coefficicents
    NTL::ZZ_pX uniform_poly(const int64_t min, const int64_t max) {
	    NTL::ZZ_pX P;
	    std::uniform_int_distribution<int64_t> dist(min, max);
	
	    P.SetLength(_phi_degree);
	    for(long i = 0 ; i < _phi_degree-1  ; ++i){
		    P[i] = dist(gen);
	    }
        P.normalize();
        return P;
    }
    
    void keygen(){
	    NTL::ZZ_pX f1, f, u, v, d;
	    NTL::ZZ_p tp;
	    NTL::conv(tp, _t);
	    do{
		    f1 = uniform_poly(-1, 1);
		    NTL::rem(f1, f1, _PhiModulus);
		    f = tp*f1+1;
		    NTL::XGCD(d, u, v, f, _Phi);
	    }while(NTL::deg(d) != 0);
	    _pk = u;
	    _sk = f;

      _pk = tp*_pk;
      _pk += 12341229385748;

      std::cout << "secret key: " << _sk << std::endl;

      NTL::ZZ_pX g = uniform_poly(-1, 1);
      NTL::rem(g, g, _PhiModulus);
      _pk = MulMod(_pk, g, _PhiModulus);

      _pk = uniform_poly(-1, 1);
      NTL::rem(_pk, _pk, _PhiModulus);

      std::cout << "multiplied: " << MulMod(_pk, _sk, _PhiModulus) << std::endl;

        //eval key
        _evk.resize(_logqw);
        
        auto tmp = power_of_poly(f);
        for(uint64_t i = 0 ; i < _logqw ; ++i){
            auto ee = gaussian_poly(0);
            auto ss = gaussian_poly(0);
            _evk[i] = tmp[i]+ee+MulModPhi(_pk, ss);
        }
    }
    
    // Suppose that m coefficient are modulo t
    Ciphertext_t encrypt(const Plainetxt_t & m) const {
        // Error polynomials
	    auto s = gaussian_poly(0);
	    auto e = gaussian_poly(0);
   
	    NTL::ZZ_p deltap;
	    NTL::conv(deltap, _delta);
	    Ciphertext_t c;
        c = deltap*m+e+MulModPhi(_pk, s);
	    return c;
    }
    
    Plainetxt_t decrypt(const Ciphertext_t & c) const {
      NTL::ZZ_pX tmp = MulMod(c, _sk, _PhiModulus);
	    NTL::ZZX tmpZZ;
	    NTL::conv(tmpZZ, tmp);
	    NTL::RR tRR, qRR, tmpRR;
	    NTL::conv(qRR, _q);

	    for(long i = 0 ; i <= NTL::deg(tmpZZ) ; ++i){
		    tmpZZ[i] *= _t;
		    NTL::conv(tmpRR, tmpZZ[i]);
		    tmpZZ[i] = NTL::RoundToZZ(tmpRR/qRR);	
	    }
	    for(long i = 0 ; i <= NTL::deg(tmpZZ) ; ++i){
		    tmpZZ[i] %= _t;
	    }
	    NTL::conv(tmp, tmpZZ);
	    return tmp;
    }
    
    Ciphertext_t relin(const Ciphertext_t & c) const {       
        // std::chrono::time_point<std::chrono::system_clock> start, end;
        // start = std::chrono::system_clock::now();
        
        auto res = word_decomp_poly(c);
        
        // end = std::chrono::system_clock::now();
 	    // std::chrono::duration<double, std::milli> elapsed_seconds = end-start;
 	    // std::cout << "time word_decomp_poly: " << elapsed_seconds.count() << "ms" << std::endl;

        // start = std::chrono::system_clock::now();               
        
        auto r = dot_prod(res, _evk);       
        // end = std::chrono::system_clock::now();
 	    // elapsed_seconds = end-start;
 	    // std::cout << "time dot_prod: " << elapsed_seconds.count() << "ms" << std::endl;
        
        return r;
    }
    
    Ciphertext_t mul(const Ciphertext_t & c1, const Ciphertext_t & c2) const {
        Ciphertext_t cmul;
        NTL::ZZX c1ZZ, c2ZZ, cmulZZ;
	    NTL::conv(c1ZZ, c1);
	    NTL::conv(c2ZZ, c2);
	    cmulZZ = c1ZZ * c2ZZ;
	
        /*
	    NTL::RR tmpRR, qRR;
	    NTL::conv(qRR, _q);
	    NTL::ZZ tmp;
	    for(long i = 0 ; i <= NTL::deg(cmulZZ) ; ++i){
		    tmp = cmulZZ[i]*_t;
		    NTL::conv(tmpRR, tmp);
		    tmpRR /= qRR;
		    NTL::round(tmpRR);
		    NTL::conv(cmulZZ[i], tmpRR);
	    }
        //*/

        NTL::ZZ tmp, quo;
        for(long i = 0 ; i <= NTL::deg(cmulZZ) ; ++i){
            tmp = cmulZZ[i]*_t;
            quo = tmp/_q;
            if((tmp<<1) >= _q*(quo+1))
                quo++;
            cmulZZ[i] = quo;
        }

        /*
        NTL::RR tRR, qRR, tqRR;
	    NTL::conv(qRR, _q);
        NTL::conv(tRR, _t);
        tqRR = tRR/qRR;
	    NTL::ZZ tmp;
        auto dd = NTL::deg(cmulZZ);
	    for(long i = 0 ; i <= dd ; ++i){
            NTL::RR tmpRR;
		    NTL::conv(tmpRR, cmulZZ[i]);
		    tmpRR *= tqRR;
		    NTL::round(tmpRR);
		    NTL::conv(cmulZZ[i], tmpRR);
	    }
        */
	    NTL::conv(cmul, cmulZZ);
        NTL::rem(cmul, cmul, _PhiModulus);
	    return cmul;
    }
    
    Ciphertext_t mulcst(const Ciphertext_t & c1, const int64_t & cst) const {
        return c1*cst;
    }
    
    Ciphertext_t add(const Ciphertext_t & c1, const Ciphertext_t & c2) const {
        Ciphertext_t r;
        NTL::add(r, c1, c2);
        return r;
    }
    
    Ciphertext_t sub(const Ciphertext_t & c1, const Ciphertext_t & c2) const {
        Ciphertext_t r;
        NTL::sub(r, c1, c2);
        return r;
    }
    
    void add(Ciphertext_t & r, const Ciphertext_t & c1, const Ciphertext_t & c2) const {
        return NTL::add(r, c1, c2);
    }
    
    void addin(Ciphertext_t & c1, const Ciphertext_t & c2) const {
        return NTL::add(c1, c1, c2);
    }
    
    std::ostream& printParameters(std::ostream & out){
		out << "==== Yashe Parameters ==== " << std::endl;
		out << "q            : " << _q << std::endl;
		out << "t            : " << _t << std::endl;
		out << "delta        : " << _delta << std::endl;
		out << "degree       : " << NTL::deg(_Phi) << std::endl;
		out << "sigma        : " << _sigma << std::endl;
		out << "base         : " << _base << std::endl;
		out << "q_bitsize    : " << _q_bitsize << std::endl;
		out << "base_bitsize : " << _base_bitsize << std::endl;
		out << "logqw        : " << _logqw << std::endl;
		out << "seed         : " << _seed << std::endl;
		out << "========================== " << std::endl;
    return out;
	}
    
};
