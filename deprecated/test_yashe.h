#pragma once

NTL::ZZ_pX encode(const Yashe & SHE, const uint64_t & n){
	NTL::ZZ_pX m;
	m.SetLength(SHE.phi_degree()-1);
	for(long i = 0 ; i < SHE.phi_degree() - 1 ; ++i){
		m[i] = (n & (1 << i)) >> i;
	}
	return m;
}

NTL::ZZ decode(const NTL::ZZ_pX & c){
	NTL::ZZ tmp, r(0);
	for(long i = 0 ; i <= NTL::deg(c) ; ++i){
		NTL::conv(tmp, c[i]);
		r += tmp*(NTL::ZZ(1) << i);
	}
	return r;
}

void test_add_int64(const Yashe & SHE){
	std::mt19937 gen(SHE._seed);
	std::uniform_int_distribution<int64_t> dist(0, (1L << 32));
	int64_t a = dist(gen);
	int64_t b = dist(gen);
	
	// cout << "a: " << a << endl;
	// cout << "b: " << b << endl;
	
	auto ma = encode(SHE, a);
	auto mb = encode(SHE, b);
	
	auto ca = SHE.encrypt(ma);
	auto cb = SHE.encrypt(mb);
	
	std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
	
	auto cadd = SHE.add(ca, cb);
	
	end = std::chrono::system_clock::now();
 	std::chrono::duration<double, std::milli> elapsed_seconds = end-start;
 	std::cout << "time add: " << elapsed_seconds.count() << "ms" << std::endl;
	
	auto madd = SHE.decrypt(cadd);
	
	auto c = decode(madd);
	
	std::cout << "int64 add: " << ((c == a+b) ? "OK" : "ERROR") << std::endl;
}

void test_mul_int64(const Yashe & SHE){
	std::mt19937 gen(SHE._seed);
	std::uniform_int_distribution<int64_t> dist(0, (1L << 10));
	int64_t a = dist(gen);
	int64_t b = dist(gen);
	
	// cout << "a: " << a << endl;
	// cout << "b: " << b << endl;
	//b = 3;
	auto ma = encode(SHE, a);
	auto mb = encode(SHE, b);
	std::chrono::time_point<std::chrono::system_clock> start, endmul, endrelin;
	auto ca = SHE.encrypt(ma);
	auto cb = SHE.encrypt(mb);
	
    start = std::chrono::system_clock::now();
	
	auto cmul = SHE.mul(ca, cb);
	// auto cmul = NTL::ZZ_p(b)*ca;
	
	endmul = std::chrono::system_clock::now();
	std::chrono::duration<double, std::milli> elapsed_seconds_mul = endmul-start;
	start = std::chrono::system_clock::now();
	
	cmul = SHE.relin(cmul);
	
	endrelin = std::chrono::system_clock::now();
 	std::chrono::duration<double, std::milli> elapsed_seconds_relin = endrelin-start;
 	std::cout << "time mul: " << elapsed_seconds_mul.count() << "ms" << std::endl;
	std::cout << "time relin: " << elapsed_seconds_relin.count() << "ms" << std::endl;
	std::cout << "time total: " << elapsed_seconds_mul.count() + elapsed_seconds_relin.count()<< "ms" << std::endl;
	
	auto mmul = SHE.decrypt(cmul);
	
	auto c = decode(mmul);
	// cout << "a*b: " << a*b << endl;
	// cout << "c: " << c << endl;
	std::cout << "int64 mul: " << ((c == a*b) ? "OK" : "ERROR") << std::endl;
}

void test_mulcst_int64(const Yashe & SHE){
	std::mt19937 gen(SHE._seed);
	std::uniform_int_distribution<int64_t> dist(0, (1L << 10));
	int64_t a = dist(gen);
	int64_t b = dist(gen);
	
	// cout << "a: " << a << endl;
	// cout << "b: " << b << endl;
	//b = 3;
	auto ma = encode(SHE, a);
	
	auto ca = SHE.encrypt(ma);
	
	std::chrono::time_point<std::chrono::system_clock> start, endmul, endrelin;
    start = std::chrono::system_clock::now();
	
	auto cmul = SHE.mulcst(ca, b);
	// auto cmul = NTL::ZZ_p(b)*ca;
	
	endmul = std::chrono::system_clock::now();
	std::chrono::duration<double, std::milli> elapsed_seconds_mul = endmul-start;
	
 	std::cout << "time mulcst: " << elapsed_seconds_mul.count() << "ms" << std::endl;
	
	auto mmul = SHE.decrypt(cmul);
	
	auto c = decode(mmul);

	std::cout << "int64 mulcst: " << ((c == a*b) ? "OK" : "ERROR") << std::endl;
}

void test_wd_po(const Yashe & SHE){
	auto a = SHE.uniform_poly(-1000, 1000);
	auto b = SHE.uniform_poly(-1000, 1000);
	
	auto c = NTL::MulMod(a, b, SHE._PhiModulus);
	auto wda = SHE.word_decomp_poly(a);
	
	auto pob = SHE.power_of_poly(b);
	auto c2 = SHE.dot_prod(wda, pob);
	std::cout << "Dot prod: " << ((c == c2) ? "OK" : "ERROR" ) << std::endl; 
}

void test_decrypt(const Yashe & SHE){
	auto m = SHE.uniform_poly(-1000, 1000);
	
	auto c = SHE.decrypt(m);
}

void test_add(const Yashe & SHE){
	auto m1 = SHE.uniform_poly(0, 1);
	auto m2 = SHE.uniform_poly(0, 1);
	
	std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
	
 	auto c1 = SHE.encrypt(m1);
	auto c2 = SHE.encrypt(m2);

	
 	end = std::chrono::system_clock::now();
 	std::chrono::duration<double, std::milli> elapsed_seconds = end-start;
 	std::cout << "time encrypt: " << elapsed_seconds.count()/2 << "ms" << std::endl;
	
	start = std::chrono::system_clock::now();
	
	auto cadd = SHE.add(c1, c2);

  std::cout << cadd << std::endl;
	
 	end = std::chrono::system_clock::now();
 	elapsed_seconds = end-start;
 	std::cout << "time add: " << elapsed_seconds.count() << "ms" << std::endl;
	
	start = std::chrono::system_clock::now();
	
	auto madd = SHE.decrypt(cadd);
	
 	end = std::chrono::system_clock::now();
 	elapsed_seconds = end-start;
 	std::cout << "time decrypt: " << elapsed_seconds.count() << "ms" << std::endl;	

  if (m1+m2 == madd) {
    std::cout << "OK" << std::endl;
  } else {
    std::cout << "Error: expected " << m1 << " + " << m2 << " = " << m1+m2 << " but got " << madd << std::endl;
  }
}

void test_mul(const Yashe & SHE){
	auto m1 = SHE.uniform_poly(0, 1);
	auto m2 = SHE.uniform_poly(0, 1);
	
	std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
	
 	auto c1 = SHE.encrypt(m1);
	auto c2 = SHE.encrypt(m2);
	
 	end = std::chrono::system_clock::now();
 	std::chrono::duration<double, std::milli> elapsed_seconds = end-start;
 	std::cout << "time encrypt: " << elapsed_seconds.count()/2 << "ms" << std::endl;
	
	start = std::chrono::system_clock::now();
	
	auto cmul = SHE.mul(c1, c2);
	
 	end = std::chrono::system_clock::now();
 	elapsed_seconds = end-start;
 	std::cout << "time mul: " << elapsed_seconds.count() << "ms" << std::endl;
	
	start = std::chrono::system_clock::now();
	
	cmul = SHE.relin(cmul);
	
 	end = std::chrono::system_clock::now();
 	elapsed_seconds = end-start;
 	std::cout << "time relin: " << elapsed_seconds.count() << "ms" << std::endl;
	
	start = std::chrono::system_clock::now();
	
	auto mmul = SHE.decrypt(cmul);
	
 	end = std::chrono::system_clock::now();
 	elapsed_seconds = end-start;
 	std::cout << "time decrypt: " << elapsed_seconds.count() << "ms" << std::endl;
}

void test_MulModPhi(const Yashe & SHE){
	auto a = SHE.uniform_poly(0, 1);
	auto b = SHE.uniform_poly(0, 1);
	
	auto c = NTL::MulMod(a, b, SHE._PhiModulus);
	
	auto d = SHE.MulModPhi(a, b);
	
	std::cout << "MulModPhi: " << ((c == d) ? "OK" : "ERROR") << std::endl;
}
