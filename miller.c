/*
 * miller.c
 *
 *  Created on: 2013-4-15
 *      Author: Suliangjian
 */


/*=======================================================
 *                     Parameters:
 * ******************************************************
 * z = -(2^62 + 2^55 + 1);
 * p := 36*z^4 + 36*z^3 + 24*z^2 + 6*z +1;
 * n := 36*z^4 + 36*z^3 + 18*z^2 + 6*z +1;
 * k = 12;
 * r = n;
 ********************************************************
 * Tower of extension
 * 		GF(p^2):  GF(p)[u]/u^2 - β where β = -1
 * 		GF(p^6):  GF(p^2)[v]/v^3 - ξ where  ξ = -1 - u
 * 		GF(p^12): GF(P^6)[w]/w^2 - γ where γ = v
 ********************************************************
 * BN Elliptic Curve
 * 			y^2 = x^3 + 2
 */


#include "miller.h"

un_bn r={{0xD, 0, 0, 0xA100, 0x10, 0, 0x8000, 0xFF9F, 0x7, 0, 0x4D80, 0xBA34, 0x1, 0x4000, 0x6482, 0x2523}, 0};;

/*
 * Jacobian coordinate(X, Y, Z) to affine coordinate(x, y)
 * x = X/Z^2
 * y = Y/Z^3
 */
point_ac jc_to_ac(point_jc a, un_bn p){
	point_ac res;
	un_bn t2,t3,t0,t1;

	t2 = mul_prime(a.z, a.z, p);		//z^2
	t3 = mul_prime(t2, a.z, p);			//z^3

	t0 = inverse_mod(t2, p);
	t1 = inverse_mod(t3, p);

	res.x = mul_prime(a.x, t0, p);		//x/z^2
	res.y = mul_prime(a.y, t1, p);		//y/z^3

	return res;
}


un_bn_ext12_t miller_function(point_jc P, point_ac_ext2 Q, un_bn p){
	un_bn_ext12_t res,t0;
	line_result t1,t2;
	point_jc T;
	point_ac P_t;
	int i,j,flag,ri;
	unint m[DIGIT_SIZE];

	flag = 0;
	//P_t = jc_to_ac(P, p);			//Convert to affine coordinate
	un_bn_cpy(&P_t.x, P.x);
	un_bn_cpy(&P_t.y, P.y);

	for (i = SIZE - 1; i >= 0; i--) {

		un_bn_binary(m, r.bn[i]);

		for (j = DIGIT_SIZE - 1; j >= 0; j--) {
			ri = m[j] & 0x1;						//r = ( r(l-1),...ri,....,r0 )2

			if (ri == 0 && flag == 0)
				continue;
			else {
				if (flag == 0) {					//The first time  having '1'
					un_bn_cpy(&T.x, P.x);
					un_bn_cpy(&T.y, P.y);
					un_bn_cpy(&T.z, P.z);			//T = P

					memset(&res, 0, sizeof(un_bn_ext12_t));
					res.W0.V0.bn0.bn[0]=1;			//res = 1

					flag = 1;
				}else {

					//t1 = tangent_line(Q, T, p);		//l(Q)
					//T  = point_double(T, p);		//T = 2*T
					point_double_and_tangentline(&T, &t1, Q, T, p);
					t0 = square_ext12(res, p);		//res^2
					res = ext12_mul_lineresult(t0, t1, p);	//res = res^2 * l(Q)

					if(ri == 1 && (i!=0||j!=0)){

						//t2 = secant_line(Q, T, P_t, p);
						//T = point_addition(T, P_t, p);
						point_addition_and_secantline(&T, &t2, Q, T, P_t, p);

						res = ext12_mul_lineresult(res, t2, p);


					}
				}
			}
		}
	}

	return res;
}


un_bn_ext12_t miller_tate(point_jc P, point_ac_ext2 Q, un_bn p){
	un_bn_ext12_t res, t0, t1, t2;

	t0 = miller_function(P, Q, p);
	//printf_mul_counts();

	t1 = first_step(t0, p);				//Final exponentiation
	//printf_mul_counts();
	t2 = second_step(t1, p);
	//printf_mul_counts();
	res = hard_step(t2, p);
	//printf_mul_counts();


	return res;
}

