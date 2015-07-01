/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* ffpack/ffpack_ludivine.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 *
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

	//---------------------------------------------------------------------
	// TURBO: rank computation algorithm
	//---------------------------------------------------------------------

	template <class Field>
	inline size_t
	TURBO (const Field& F, const size_t M, const size_t N,
	       typename Field::Element* A, const size_t lda, size_t * P, size_t * Q, const size_t cutoff)
	{

		size_t mo2 = (M>>1);
		size_t no2 = (N>>1);

		typename Field::Element * NW = A;
		typename Field::Element * NE = A + no2;
		typename Field::Element * SW = A + mo2*lda;
		typename Field::Element * SE = SW + no2;

		size_t ld1, ld2, ld3, ld4;
		ld1 = ld2 = ld3 = ld4 = lda;

		if ( !(M && N) ) return 0;

		// Column permutation
		size_t * P1 = new size_t[no2];
		size_t * P2 = new size_t[N-no2];
		// Row Permutation
		size_t * Q1 = new size_t[mo2];
		size_t * Q2 = new size_t[M-mo2];
		for (size_t i=0; i<mo2; ++i)
			Q1[i] = 0;
		for (size_t i=0; i<M-mo2; ++i)
			Q2[i] = 0;
		size_t q1,q2,q3,q3b,q4;
		/*q1=q2=q3=  */q3b=q4=0;


		// Step 1: NW = L1.Q1.U1.P1
		size_t mloc = mo2;
		size_t nloc ;
#if 0
		Timer tim;
		tim.clear();
		tim.start();
#endif
		q1 = LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mloc, no2, NW, ld1, P1, Q1, FFPACK::FfpackLQUP, cutoff);

#if 0
		tim.stop();
		cerr<<"LQUP1:"<<tim.realtime()<<std::endl;
		tim.start();
#endif
#if LB_DEBUG
		std::cerr<<"NW= L1.Q1.U1.P1"<<std::endl;
		write_field(F,std::cerr,NW,M,N,lda);
#endif
		// B1 = L^-1.NE
#ifdef LB_DEBUG
		std::cerr<<"avant B1 = L^-1.NE"<<std::endl;
		write_field(F,std::cerr,NE,mloc,N-no2,ld2);
#endif
		solveLB( F, FFLAS::FflasLeft, mo2, N-no2, q1, NW, ld1, Q1, NE, ld2);
#ifdef LB_DEBUG
		std::cerr<<"B1 = L^-1.NE"<<std::endl;
		write_field(F,std::cerr,NE,mloc,N-no2,ld2);
#endif

		// NE = Q^-1.NE

		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
				N-no2, 0,(int) mo2, NE, ld2, Q1);
#ifdef LB_DEBUG
		std::cerr<<"NE=Q^-1.NE"<<std::endl;
		write_field(F,std::cerr,NE,mloc,N-no2,ld2);
#endif

		// SW = SW.P1
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
				M-mo2, 0,(int) q1, SW, ld3, P1 );
#ifdef LB_DEBUG
		std::cerr<<"SW = SW.P1"<<std::endl;
		write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif

#if 0
		tim.stop();
		std::cerr<<"L^-1:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		// N1 = SW_{1,q1} . U1^-1
		ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, M-mo2, q1, F.one, NW, ld1 , SW, ld3 );
#ifdef LB_DEBUG
		std::cerr<<" N1 = SW_{1,q1} . U1^-1"<<std::endl;
		write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif

#if 0
		tim.stop();
		std::cerr<<"trsm:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		// I1 = SW_{q1+1,n} - N1.G1
		fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-mo2,  no2-q1, q1, F.mOne, SW, ld3, NW+q1, ld1, F.one, SW+q1, ld3);
#ifdef LB_DEBUG
		std::cerr<<" I1 = SW_{q1+1,n} - N1.G1"<<std::endl;
		write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif

#if 0
		tim.stop();
		std::cerr<<"fgemm1:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		// E1 = SE - N1.B1_{1,q1}
		fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-mo2, N-no2, q1, F.mOne, SW, ld3, NE, ld2, F.one, SE, ld4);
#ifdef LB_DEBUG
		std::cerr<<"  E1 = SE - N1.B1_{1,q1}"<<std::endl;
		write_field(F,std::cerr,SE,M-mo2,N-no2,ld4);
#endif

#if 0
		tim.stop();
		std::cerr<<"fgemm2:"<<tim.realtime()<<std::endl;
		tim.start();
#endif


		//Step 2: E1 = L2.Q2.U2.P2
		mloc = M-mo2;
		nloc = N-no2;
		q2 = LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mloc, nloc, SE, ld4, P2, Q2, FFPACK::FfpackLQUP, cutoff);
#ifdef LB_DEBUG
		std::cerr<<"  E1 = L2.Q2.U2.P2"<<std::endl;
		write_field(F,std::cerr,SE,M-mo2,N-no2,ld4);
#endif

#if 0
		tim.stop();
		std::cerr<<"LQUP2:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		// [I2;F2] = L2^-1.I1
		solveLB( F, FFLAS::FflasLeft, mloc, no2-q1, q2, SE, ld4, Q2, SW+q1, ld3);
#ifdef LB_DEBUG
		std::cerr<<"  [I2;F2] = L2^-1.I1"<<std::endl;
		write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif
		// I1 = Q2^-1.I1
		FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
				no2-q1, 0,(int) mloc, SW+q1, ld3, Q2 );
#ifdef LB_DEBUG
		std::cerr<<"I1 = Q2^-1.I1"<<std::endl;
		write_field(F,std::cerr,SW,mloc,no2,ld3);
#endif

		// B1 = B1.P2
		FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
				mo2, 0,(int) q2, NE, ld2, P2 );
#ifdef LB_DEBUG
		std::cerr<<"B1 = B1.P2"<<std::endl;
		write_field(F,std::cerr,NE,mo2,N-no2,ld2);
#endif
		// Updating P
#if 0
		for (size_t i=no2;i<N;++i)
			P[i] += no2;
		tim.stop();
		std::cerr<<"L2^-1:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		//alternative: de 0 a q2 avant
		// N2 = B1_{q1+1,mo2} . V2^-1
		ftrsm(F, FFLAS::FflasRight, FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit, mo2-q1, q2, F.one, SE, ld4, NE+q1*ld2,ld2);
#if 0
		tim.stop();
		std::cerr<<"trsm2:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		// H2 = B1_{q1+1,mo2;q2,N-no2} - N2.E2
		fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mo2-q1, N-no2-q2, q2, F.mOne, NE+q1*ld2, ld2, SE+q2, ld4, F.one, NE+q1*ld2+q2, ld2);

#if 0
		tim.stop();
		std::cerr<<"fgemm12:"<<tim.realtime()<<std::endl;
		tim.start();
		O2 = NW_{q1+1,mo2;q1+1,N-no2} = - N2.I2
		write_field (F,cerr<<"avant O2"<<endl, A, M, N, lda);
#endif

		fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mo2-q1, no2-q1, q2, F.mOne, NE+q1*ld2, ld2, SW+q1, ld3, F.zero,
		      NW+q1*(ld1+1), ld1);
		//	write_field (F,cerr<<"apres O2"<<endl, A, M, N, lda);
#if 0
		tim.stop();
		std::cerr<<"fgemm22:"<<tim.realtime()<<std::endl;
		tim.start();
#endif


		//Step 3: F2 = L3.Q3.U3.P3
		mloc = M-mo2-q2;
		nloc = no2-q1;
		q3 = LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mloc, nloc, SW+q2*ld3+q1, ld3, P1+q1, Q2+q2, FFPACK::FfpackLQUP, cutoff);

		// Updating P1,Q2
		for (size_t i=q1;i<no2;++i)
			P1[i] += q1;
		for (size_t i=q2;i<q2+q3;++i)
			Q2[i] += q2;

		//Step 3bis: H2 = L3b.Q3b.U3b.P3b
		mloc = mo2-q1;
		nloc = N-no2-q2;

		q3b = LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mloc, nloc, NE+q1*ld2+q2, ld2, P2+q2, Q1+q1, FFPACK::FfpackLQUP, cutoff);

		// Updating P2, Q1
		for (size_t i = q2; i < q2+q3b; ++i)
			P2[i] += q2;

#if 0
		tim.stop();
		std::cerr<<"LQUP3et3bis:"<<tim.realtime()<<std::endl;
		tim.start();
#endif

		if (( q3 < no2-q1) && (q3b<mo2-q1)){

			// [O3;_] = L3b^-1.O2
			if (q3b>0){
#if 0
				if ( mo2-q1 < N-no2-q2+q1)
					// L is expanded to a Lower triangular matrix
					solveLB( F, FFLAS::FflasLeft,mloc, no2-q1, q3b, NE+q1*ld2+q2 , ld2, rP3b, NW+q1*(ld1+1), ld1);
				else
#endif
					//std::cerr<<"USING SOLVELB2"<<std::endl;
					//no modification of L
					solveLB2( F, FFLAS::FflasLeft,mloc, no2-q1, q3b, NE+q1*ld2+q2 , ld2, Q1+q1, NW+q1*(ld1+1), ld1);
#ifdef LB_DEBUG
				std::cerr<<"O2 avant="<<std::endl;
				write_field(F,std::cerr,NW+q1*(ld1+1),mloc,no2-q1,ld1);
#endif

				// O2 = Q3b^-1.O2
				FFPACK::applyP( F, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
						no2-q1, 0,(int) mloc, NW+q1*(ld1+1), ld1, Q1+q1 );
#ifdef LB_DEBUG
				std::cerr<<"O2 apres="<<std::endl;
				write_field(F,std::cerr,NW+q1*(ld1+1),mloc,no2-q1,ld1);
#endif

				//updating Q
#if 0
				size_t tmp;
				for (size_t j=0;j<mo2-q1;++j)
					if (rP3b[j]!=j){
						//	std::cerr<<"(rP3b["<<j<<"]="<<rP3b[j]<<std::endl;
						tmp = Q[j+q1];
						Q[j+q1] = Q[rP3b[j]+q1];
						Q[rP3b[j]+q1] = tmp;
					}
#endif

				// X2 = X2.P3
				// Si plusieurs niveaux rec, remplacer X2 par [NW;I2]
				FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
						mo2-q1-q3b,(int) q1, (int)(q1+q3),
						NW/*+(q1+q3b)*ld1*/, ld1, P1);
				FFPACK::applyP( F, FFLAS::FflasRight, FFLAS::FflasTrans,
						q2,(int) q1, (int)(q1+q3),
						SW/*+(q1+q3b)*ld1*/, ld3, P1);


				// A faire si plusieurs niveaux recursifs
				// B2 = B2.P3b
				FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
						q1,(int) q2, (int)(q2+q3b),
						NW, ld2, P2);
				//flaswp(F,q1,NE,lda,no2+q2,no2+q2+q3b,P,1);
				// E2 = E2.P3b
				FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
						q2,(int) q2, (int)(q2+q3b),
						SE, ld4, P2);
				//flaswp(F,q2,SE+q2,lda,no2+q2,no2+q2+q3b,P,1);
			}

			// N3 = X2 . D3^-1
			ftrsm( F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit, mo2-q1-q3b, q3, F.one, SW+q2*ld3+q1, ld3 ,NW+(q1+q3b)*ld1+q1,ld1);

			// T2 = T2 - N3.F3
			fgemm( F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, mo2-q1-q3b, no2-q1-q3,q3, F.mOne, NW+(q1+q3b)*ld1+q1, ld1, SW+q2*ld3+q3+q1, ld3, F.one, NW+(q1+q3b)*ld1+q1+q3, ld1 );


			//Step 4: T2 = L4.Q4.U4.P4
			mloc = mo2-q1-q3b;
			nloc = no2-q1-q3;

#if 0
			size_t * rP4 = new size_t[mloc];
			for (size_t j=0;j<mo2-q1;++j)
				rP4[j]=0;
#endif
			q4 = LUdivine( F, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, mloc, nloc, NW+(q1+q3b)*ld1+q1+q3, ld1, P1+q1+q3, Q1+q1+q3b, FFPACK::FfpackLQUP, cutoff);

			// Updating P
			for (size_t i=q1+q3;i<q1+q3+q4;++i)
				P1[i] += q3;

#if 0
			size_t tmp;
			if (rP4[j]!=j){
				//	std::cerr<<"(rP3b["<<j<<"]="<<rP3b[j]<<std::endl;
				tmp = Q[j+q1+q3b];
				Q[j+q1+q3b] = Q[rP3b[j]+q1+q3b];
				Q[rP3b[j]+q1+q3b] = tmp;
			}
#endif

			// A faire si plusieurs niveaux recursifs
			// [G1;O3] = [G1;O3].P4
			FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
					q1+q3b, (int)(q1+q3), (int)(q1+q3+q4),
					NW, ld1, P1);
			//flaswp(F,q1+q3b,NE,lda,no2+q2,no2+q2+q3b,P,1);
			// [I2;F3] = [I2;F3].P4
			FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasTrans,
					q2+q3, (int)(q1+q3),(int) (q1+q3+q4),
					SW, ld3, P1);
			//flaswp(F,q2,SE+q2,lda,no2+q2,no2+q2+q3b,P,1);
		}
		//!!!!!! Attention a appliquer Q4, Q2, Q3, Q3b a gauche !!!!!!!

		//updating Q1
		for (size_t i = q1; i < q1+q3b; ++i)
			Q1[i] += q1;
		for (size_t i=q1+q3b;i<q1+q3b+q4;++i)
			Q1[i] += q1 + q3b;

		for (size_t i=0; i<q1; ++i)
			P[i] = P1[i];
		for (size_t i=q1; i<q1+q2; ++i)
			P[i] = P2[i-q1] + no2;
		for (size_t i=q1+q2; i<q1+q2+q3; ++i)
			P[i] = P1[i-q2];
		for (size_t i=q1+q2+q3; i<q1+q2+q3+q3b; ++i)
			P[i] = P2[i-q1-q3]+no2;
		for (size_t i=q1+q2+q3+q3b; i<q1+q2+q3+q3b+q4; ++i)
			P[i] = P1[i-q2-q3b];
		delete[] P1;
		delete[] P2;

		for (size_t i=0; i<q1; ++i)
			Q[i] = Q1[i];
		for (size_t i=q1; i<q1+q2; ++i)
			Q[i] = Q2[i-q1] + mo2;
		for (size_t i=q1+q2; i<q1+q2+q3; ++i)
			Q[i] = Q2[i-q1] + mo2;
		for (size_t i=q1+q2+q3; i<q1+q2+q3+q3b; ++i)
			Q[i] = Q1[i-q2-q3];
		for (size_t i=q1+q2+q3+q3b; i<q1+q2+q3+q3b+q4; ++i)
			P[i] = Q1[i-q2-q3];
		delete[] Q1;
		delete[] Q2;


		//write_field (F, cerr<<"avant reordonnancement"<<endl, A, M,N, lda)<<endl;
		typename Field::Element * R = new typename Field::Element[M*N];
		size_t ldr = N;
		// Copying first q1 cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, q1, NW+i*ld1,1, R+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, q1, SW+(i-q1)*ld3,1, R+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, q1, NW+(i-q2-q3)*ld1,1, R+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, q1, SW+(i-mo2)*ld3,1, R+i*ldr, 1);
		// Copying q1..q2 cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, q2, NE+i*ld2,1, R+q1+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, q2, SE+(i-q1)*ld4,1, R+q1+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, q2, NE+(i-q2-q3)*ld2,1, R+q1+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, q2, SE+(i-mo2)*ld4,1, R+q1+i*ldr, 1);
		// Copying q2..q3 cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, q3, NW+q1+i*ld1,1, R+q1+q2+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, q3, SW+q1+(i-q1)*ld3,1, R+q1+q2+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, q3, NW+q1+(i-q2-q3)*ld1,1, R+q1+q2+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, q3, SW+q1+(i-mo2)*ld3,1, R+q1+q2+i*ldr, 1);
		// Copying q3..q3b cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, q3b, NE+q2+i*ld2,1, R+q1+q2+q3+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, q3b, SE+q2+(i-q1)*ld4,1, R+q1+q2+q3+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, q3b, NE+q2+(i-q2-q3)*ld2,1, R+q1+q2+q3+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, q3b, SE+q2+(i-mo2)*ld4,1, R+q1+q2+q3+i*ldr, 1);
		// Copying q3b..q4 cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, q4, NW+q1+q3+i*ld1,1, R+q1+q2+q3+q3b+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, q4, SW+q1+q3+(i-q1)*ld3,1, R+q1+q2+q3+q3b+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, q4, NW+q1+q3+(i-q2-q3)*ld1,1, R+q1+q2+q3+q3b+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, q4, SW+q1+q3+(i-mo2)*ld3,1, R+q1+q2+q3+q3b+i*ldr, 1);
		// Copying the last cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, no2-q1-q3-q4, NW+q1+q3+q4+i*ld1,1, R+q1+q2+q3+q3b+q4+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, no2-q1-q3-q4, SW+q1+q3+q4+(i-q1)*ld3,1, R+q1+q2+q3+q3b+q4+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, no2-q1-q3-q4, NW+q1+q3+q4+(i-q2-q3)*ld1,1, R+q1+q2+q3+q3b+q4+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, no2-q1-q3-q4, SW+q1+q3+q4+(i-mo2)*ld3,1, R+q1+q2+q3+q3b+q4+i*ldr, 1);
		// Copying the last cols
		for (size_t i=0; i<q1; ++i)
			FFLAS::fcopy (F, N-no2-q2-q3b, NE+q2+q3b+i*ld2,1, R+no2+q2+q3b+i*ldr, 1);
		for (size_t i=q1; i<q1+q2+q3; ++i)
			FFLAS::fcopy (F, N-no2-q2-q3b, SE+q2+q3b+(i-q1)*ld4,1, R+no2+q2+q3b+i*ldr, 1);
		for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
			FFLAS::fcopy (F, N-no2-q2-q3b, NE+q2+q3b+(i-q2-q3)*ld2,1, R+no2+q2+q3b+i*ldr, 1);
		for (size_t i=q2+q3+mo2; i<M; ++i)
			FFLAS::fcopy (F, N-no2-q2-q3b, SE+q2+q3b+(i-mo2)*ld4,1, R+no2+q2+q3b+i*ldr, 1);

		// A=R : to be improved (avoid allocation of R). To be changed if rec data structure are used
		for (size_t i=0; i<M; ++i)
			FFLAS::fcopy (F, N, R+i*ldr,1, A+i*lda, 1);

		delete[] R;
		//delete[] Q;
		// Necessaire:
		// 1 traiter les flaswp manquants
		// Facultatif:
		// 2 permutations de lignes doivent etre coherentes
		// 3 effectuer les dernieres permutations lignes et colonnes
		//std::cerr<<q1<<" "<<q2<<" "<<q3<<" "<<q3b<<" "<<q4<<std::endl;
		return q1+q2+q3+q3b+q4;
	}
#if 0
		//---------------------------------------------------------------------
		// RectangleCopyTURBO: Copy A to T, with respect to the row permutation
		//                     defined by the lsp factorization of located in
		//                     A-dist2pivot
		//---------------------------------------------------------------------
		template <class Field>
		void
		RectangleCopyTURBO( const Field& F, const size_t M, const size_t N,
				    const size_t dist2pivot, const size_t rank,
				    typename Field::Element * T, const size_t ldt,
				    const typename Field::Element * A, const size_t lda )
		{

			const typename Field::Element * Ai = A;
			typename Field::Element * T1i = T, T2i = T + rank*ldt;
			size_t x = dist2pivot;
			for (; Ai<A+M*lda; Ai+=lda){
				while ( F.isZero(*(Ai-x)) ) { // test if the pivot is 0
					FFLAS::fcopy( F, N, Ai, 1, T2i, 1);
					Ai += lda;
					T2i += ldt;
				}
				FFLAS::fcopy( F, N, Ai, 1, T1i, 1);
				T1i += ldt;
				x--;
			}
		}
#endif
