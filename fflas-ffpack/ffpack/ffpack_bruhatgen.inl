void bruhat_gen (const Field& Fi,
           const size_t M, const size_t N,
           typename Field::Element_ptr A, const size_t lda,
           size_t BCThreshold)
    {
        if (std::max (M,N) == 1) Fi.assign(*(A+0), Fi.zero);
	

        FFLAS::FFLAS_DIAG OppDiag = (Diag == FFLAS::FflasUnit)? FFLAS::FflasNonUnit : FFLAS::FflasUnit;
        size_t M2 = M >> 1; //on divise par 2 M= Ligne
        size_t N2 = N >> 1; // N = colonnes
        size_t * P1 = FFLAS::fflas_new<size_t >(M2);
        size_t * Q1 = FFLAS::fflas_new<size_t >(N2);
        size_t r1,r2,r3,r4;

        // A1 = P1 [ L1 ] [ U1 V1 ] Q1
        //         [ M1 ]
        r1 = _PLUQ (Fi, Diag, M2, N2, A, lda, P1, Q1, BCThreshold);
        typename Field::Element_ptr A2 = A + N2;
        typename Field::Element_ptr A3 = A + M2*lda;
        typename Field::Element_ptr F = A2 + r1*lda;
        typename Field::Element_ptr G = A3 + r1;
        // [ B1 ] <- P1^T A2
        // [ B2 ]
#ifdef MONOTONIC_APPLYP
        MonotonicApplyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, size_t(0), M2, A2, lda, P1, r1);
        MonotonicApplyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, size_t(0), N2, A3, lda, Q1, r1);
#else
        applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasNoTrans, N-N2, size_t(0), M2, A2, lda, P1);
        // [ C1 C2 ] <- A3 Q1^T
        applyP (Fi, FFLAS::FflasRight, FFLAS::FflasTrans, M-M2, size_t(0), N2, A3, lda, Q1);
#endif
        // D <- L1^-1 B1
        ftrsm (Fi, FFLAS::FflasLeft, FFLAS::FflasLower, FFLAS::FflasNoTrans, OppDiag, r1, N-N2, Fi.one, A, lda, A2, lda);
        // E <- C1 U1^-1
        ftrsm (Fi, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, Diag, M-M2, r1, Fi.one, A, lda, A3, lda);
        // F <- B2 - M1 D
        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M2-r1, N-N2, r1, Fi.mOne, A + r1*lda, lda, A2, lda, Fi.one, A2+r1*lda, lda);
        // G <- C2 - E V1
        fgemm (Fi, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M-M2, N2-r1, r1, Fi.mOne, A3, lda, A+r1, lda, Fi.one, A3+r1, lda);
		
		//Definition de H et [D 0]
		float *H, *D;
		H= FFLAS::fflas_new(Fi, M2, N-N2);
		D= FFLAS::fflas_new(Fi, M2, N-N2); 
		typename Field::Element_ptr H2 = H + r1*(N-N2);
		FFLAS::fassign(Fi, M2-r1,N-N2, F, lda, H2, N-N2);
		FFLAS::fassign(Fi, r1,N-N2, A2, lda, D, N-N2);
		
		//On applique P
		applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, N-N2, size_t(0), M2, H, N-N2, P1);
		applyP (Fi, FFLAS::FflasLeft, FFLAS::FflasTrans, N-N2, size_t(0), M2, D, N-N2, P1);
		//On complete A2
		FFLAS::fadd(Fi, M2, N-N2, Left(D) ,N-N2, bruhat_gen(Fi,M2, N-N2, H,N-N2, size_t),N-N2, A2, lda);
		fflas_delete(H);
		fflas_delete(D);
		//Definition de I et [E 0]
		float *I, *E;
		I= FFLAS::fflas_new(Fi, M-M2, N2);
		E= FFLAS::fflas_new(Fi, M-M2,N2);
		typename Field::Element_ptr I2 = I + r1;
		FFLAS::fassign(Fi, M-M2,N2-r1, G, lda, I2, N2);
		FFLAS::fassign(Fi, M-M2, r1, A3, lda, E, N2);
		//On applique Q
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, M-M2, size_t(0), N2, I, N2, Q1);
		applyP (Fi, FFLAS::FflasRight, FFLAS::FflasNoTrans, M-M2, size_t(0), N2, E, N2, Q1);
		
		FFLAS::fadd(Fi, M-M2, N2, Left(E) ,N2, bruhat_gen(Fi,M-M2, N2, I,N2, size_t),N2, A3, lda);
		fflas_delete(I);
		fflas_delete(E);
		
		//Definition de R1 :
		float R1*;
		FFLAS::fflas_new(Fi, M2,N2);
		for (size_t i=0; i<r1; ++i) 
		{Fi.assign(*(R1+i+i*N2), Fi.one);}
		FFLAS::fadd(Fi,M2,N2, R1, N2, A1,lda,A1,lda);
		applyP(Fi,FFLAS::FflasLeft,FFLAS::FflasTrans, N2, size_t(0), M2, A1, lda, P1);
		applyP(Fi, FFLAS::FflasRight,FFLAS::FflasNoTrans, M2,size_t(0), N2, A1, lda, Q1);

    }

   