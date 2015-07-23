/*
 * final_exp.c
 *
 *  Created on: 2013-4-3
 *      Author: Suliangjian
 */

#include "final_exp.h"


//precomputed
un_bn r1_p2={{0x8, 0, 0, 0xCD80, 0x6, 0, 0, 0x4909, 0x2, 0, 0x6240, 0x49B3, 0, 0, 0, 0}, 0};
un_bn r2_p2={{0x7, 0, 0, 0xCD80, 0x6, 0, 0, 0x4909, 0x2, 0, 0x6240, 0x49B3, 0, 0, 0, 0}, 0};
un_bn r3_p2={{0x12, 0, 0, 0xA700, 0x13, 0, 0, 0x6121, 0x8, 0, 0x4D80, 0xBA34, 0x1, 0x4000, 0x6482, 0x2523}, 0};
un_bn r4_p2={{0xB, 0, 0, 0xD980, 0xC, 0, 0, 0x1818, 0x6, 0, 0xEB40, 0x7080, 0x1, 0x4000, 0x6482, 0x2523}, 0};
un_bn r5_p2={{0xC, 0, 0, 0xD980, 0xC, 0, 0, 0x1818, 0x6, 0, 0xEB40, 0x7080, 0x1, 0x4000, 0x6482, 0x2523}, 0};

un_bn_ext2 r1_p={
		{{0x922A, 0x90D5, 0x193F, 0xC582, 0x8850, 0xB2C0, 0x8B6D, 0xDC17, 0x6AC8, 0x57B9, 0xB22F, 0x03EA, 0x8375, 0x1ED1, 0xEE69, 0x09EB}, 0},
		{{0x6DE9, 0x6F2A, 0xE6C0, 0xE17D, 0x77C2, 0x4D3F, 0x7492, 0x8509, 0x953F, 0xA846, 0x9B50, 0xB649, 0x7C8C, 0x212E, 0x7619, 0x1B37}, 0}};
un_bn r2_p={{0x000B, 0x0000, 0x0000, 0xD980, 0x000C, 0x0000, 0x0000, 0x1818, 0x0006, 0x0000, 0xEB40, 0x7080, 0x0001, 0x4000, 0x6482, 0x2523}, 0};		//0+a0*u
un_bn_ext2 r3_p={
		{{0x7E4E, 0x1A74, 0x7115, 0x5BE4, 0x7354, 0x9D28, 0xC5F1, 0xB9ED, 0x5F92, 0x7B75, 0xC5D7, 0xF398, 0xB248, 0x9C60, 0x9AB0, 0x0143}, 0},
		{{0x7E4E, 0x1A74, 0x7115, 0x5BE4, 0x7354, 0x9D28, 0xC5F1, 0xB9ED, 0x5F92, 0x7B75, 0xC5D7, 0xF398, 0xB248, 0x9C60, 0x9AB0, 0x0143}, 0}};
un_bn r4_p={{0x000C, 0x0000, 0x0000, 0xD980, 0x000C, 0x0000, 0x0000, 0x1818, 0x0006, 0x0000, 0xEB40, 0x7080, 0x0001, 0x4000, 0x6482, 0x2523}, 0};		//a0+0*u
un_bn_ext2 r5_p={
		{{0x1078, 0xAB4A, 0x8A54, 0x2166, 0xFBA5, 0x4FE8, 0x515F, 0x9605, 0xCA5B, 0xD32E, 0x7806, 0xF783, 0x35BD, 0xBB32, 0x8919, 0x0B2F}, 0},
		{{0xEF9B, 0x54B5, 0x75AB, 0x8599, 0x046E, 0xB017, 0xAEA0, 0xCB1B, 0x35AC, 0x2CD1, 0xD579, 0xC2B0, 0xCA43, 0x84CD, 0xDB68, 0x19F3}, 0}};


//-z = 2^62 + 2^55 + 1
int neg_z[64]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1};


/*
 * f∈GF(p^12)
 *
 *       return f^p
 */
un_bn_ext12_t ext12_frobenius_p(un_bn_ext12_t f, un_bn p){
	un_bn_ext12_t res;
	un_bn_ext2 f0,f1,f2,f3,f4,f5;

	//fi^(p) = conjugation(fi)
	un_bn_cpy(&f0.bn0, f.W0.V0.bn0);
	f0.bn1 = neg_mod(f.W0.V0.bn1, p);			//f0
	un_bn_cpy(&f1.bn0, f.W1.V0.bn0);
	f1.bn1 = neg_mod(f.W1.V0.bn1, p);			//f1
	un_bn_cpy(&f2.bn0, f.W0.V1.bn0);
	f2.bn1 = neg_mod(f.W0.V1.bn1, p);			//f2
	un_bn_cpy(&f3.bn0, f.W1.V1.bn0);
	f3.bn1 = neg_mod(f.W1.V1.bn1, p);			//f3
	un_bn_cpy(&f4.bn0, f.W0.V2.bn0);
	f4.bn1 = neg_mod(f.W0.V2.bn1, p);			//f4
	un_bn_cpy(&f5.bn0, f.W1.V2.bn0);
	f5.bn1 = neg_mod(f.W1.V2.bn1, p);			//f5

	//f*ri_p
	//res.W0.V0 = mul_ext2(f0, r1_p, p);		//f0*r1_p
	ext2_cop(&res.W0.V0, f0);
	res.W1.V0 = mul_ext2(f1, r1_p, p);			//f1*r1_p

	//(f20+f21*u)*(0+r2_p*u) = (β*f21*r2_p) + (f20*r2_p)*u, β=-1
	res.W0.V1.bn0 = mul_prime(f2.bn1, r2_p, p);
	res.W0.V1.bn0 = neg_mod(res.W0.V1.bn0, p);
	res.W0.V1.bn1 = mul_prime(f2.bn0, r2_p, p);	//f2*r2_p

	res.W1.V1 = mul_ext2(f3, r3_p, p);			//f3*r3_p
	res.W0.V2 = ext2_mul_prime(f4, r4_p, p);	//f4*r4_p
	res.W1.V2 = mul_ext2(f5, r5_p, p);			//f5*r5_p

	return res;
}

/*
 * f∈GF(p^12)
 *
 *      return f^(p^2)
 */
un_bn_ext12_t ext12_frobenius_p2(un_bn_ext12_t f, un_bn p){
	un_bn_ext12_t res;

	//fi^(p^2) = fi
	ext2_cop(&res.W0.V0, f.W0.V0);					//f0

	res.W1.V0.bn0 = mul_prime(f.W1.V0.bn0, r1_p2, p);
	res.W1.V0.bn1 = mul_prime(f.W1.V0.bn1, r1_p2, p);	//f1*r1_p2

	res.W0.V1.bn0 = mul_prime(f.W0.V1.bn0, r2_p2, p);
	res.W0.V1.bn1 = mul_prime(f.W0.V1.bn1, r2_p2, p);	//f2*r2_p2

	res.W1.V1.bn0 = mul_prime(f.W1.V1.bn0, r3_p2, p);
	res.W1.V1.bn1 = mul_prime(f.W1.V1.bn1, r3_p2, p);	//f3*r3_p2

	res.W0.V2.bn0 = mul_prime(f.W0.V2.bn0, r4_p2, p);
	res.W0.V2.bn1 = mul_prime(f.W0.V2.bn1, r4_p2, p);	//f4*r4_p2

	res.W1.V2.bn0 = mul_prime(f.W1.V2.bn0, r5_p2, p);
	res.W1.V2.bn1 = mul_prime(f.W1.V2.bn1, r5_p2, p);	//f5*r5_p2

	return res;
}

/*
 * f∈GF(p^12)
 *
 *      return f^a
 */
un_bn_ext12_t ext12_power_ext(un_bn_ext12_t f, un_bn a, un_bn p){
	un_bn_ext12_t res;
	unint m[DIGIT_SIZE];
	int i,j,flag,ai;
	un_bn_ext12_t q1,q2,t;

	flag = 0;

	for(i=SIZE-1; i>=0; i--){
		un_bn_binary(m, a.bn[i]);
		for(j=DIGIT_SIZE-1; j>=0; j--){
			ai = m[j]&0x1;
			if(ai == 0 && flag == 0)
				continue;
			else{
				if(flag == 0){				//The first time to have '1'
					ext12_cop(&q1, f);
					t = square_ext12(f, p);
					ext12_cop(&q2, t);
					flag = 1;
				}else{
					if(ai == 1){
						q1 = mul_ext12(q1, q2, p);
						q2 = square_ext12(q2, p);
					}else{
						q2 = mul_ext12(q1, q2, p);
						q1 = square_ext12(q1, p);
					}
				}
			}
		}
	}
	ext12_cop(&res, q1);

	return res;
}

/*
 * f∈GF(p^12), e=e_0 + e_1*2 + e_2*(2^2)+ ... + e_(L-1)*(2^(L-1)), e_i={-1,0,1}
 *
 *    return f^e
 */
un_bn_ext12_t ext12_power_e(un_bn_ext12_t f, un_bn p){
	un_bn_ext12_t res,t;
	int i;

	ext12_cop(&res, f);
	ext6_cop(&t.W0, f.W0);
	t.W1 = neg_ext6(f.W1, p);

	for(i=61; i>=0; i--){
		res = square_ext12(res, p);

		if(neg_z[i] != 0){
			if(neg_z[i] == 1){
				res = mul_ext12(res, f, p);
			}else{
				res = mul_ext12(res, t, p);
			}
		}
	}

	return res;
}

/*
 * f∈GF(p^12)
 *
 *     return f^(p^6-1)
 */
un_bn_ext12_t first_step(un_bn_ext12_t f, un_bn p){
	un_bn_ext12_t a,t0,t1;
	un_bn_ext12_t res;

	ext12_cop(&a, f);
	ext6_cop(&t0.W0, a.W0);
	t0.W1 = neg_ext6(a.W1, p);		//t0 = f^(p^6), conjugation(f)

	t1 = inverse_ext12(f, p);		//t1 = f^-1

	res = mul_ext12(t0, t1, p);		//unitary

	return res;
}

/*
 * f∈GF(p^12)
 *
 *        return f^(p^2 + 1)
 */
un_bn_ext12_t second_step(un_bn_ext12_t f, un_bn p){
	un_bn_ext12_t res,t;

	//f^(p^2)
	t = ext12_frobenius_p2(f, p);
	res = mul_ext12(t, f, p);

	return res;
}

/*
 * f∈GF(p^12)
 *
 *      return f^((p^4 - p^2 + 1)/r)
 *
 * After first step, f is unitary.
 */
un_bn_ext12_t hard_step(un_bn_ext12_t f, un_bn p){
	un_bn_ext12_t res,t1,t2,t3,t0,t;
	un_bn_ext12_t y0,y1,y2,y3,y4,y5,y6;

	//z = -(2^62 + 2^55 + 1)
	t1 = ext12_power_e(f, p);		//t1 = f^(-z)
	//printf("1:");printf_mul_counts();
	t2 = ext12_power_e(t1, p);		//t2 = f^(z^2)
	//printf("2:");printf_mul_counts();
	t3 = ext12_power_e(t2, p);		//t3 = f^(-z^3)
	//printf("3:");printf_mul_counts();

	t = ext12_frobenius_p2(f, p);		//t = f^(p^2)
	//printf("4:");printf_mul_counts();
	t0 = ext12_frobenius_p(f, p);		//t1 = f^p
	//printf("5:");printf_mul_counts();
	y0 = ext12_frobenius_p2(t0, p);		//y0 = f^(p^3)
	//printf("6:");printf_mul_counts();
	y0 = mul_ext12(y0, t0, p);			//y0 = f^(p^3) * f^p
	//printf("7:");printf_mul_counts();
	y0 = mul_ext12(y0, t, p);			//y0 = f^(p^3) * f^p * f^(p^2)
	//printf("8:");printf_mul_counts();

	ext6_cop(&y1.W0, f.W0);
	y1.W1 = neg_ext6(f.W1, p);			//y1 = f^-1  conjugation(f)
	//printf("9:");printf_mul_counts();

	y2 = ext12_frobenius_p2(t2, p);		//y2 = (f^(z^2))^(p^2)
	//printf("10:");printf_mul_counts();

	y3 = ext12_frobenius_p(t1, p);		//y3 = (f^-z)^p
	//printf("11:");printf_mul_counts();

	y4 = ext12_frobenius_p(t2, p);		//y4 = (f^(z^2))^p
	ext6_cop(&y4.W0, y4.W0);
	y4.W1 = neg_ext6(y4.W1, p);			//y4 = 1 / (f^(z^2))^p
	y4 = mul_ext12(t1, y4, p);			//y4 = f^(-z) / (f^(z^2))^p
	//printf("12:");printf_mul_counts();

	ext6_cop(&y5.W0, t2.W0);
	y5.W1 = neg_ext6(t2.W1, p);			//y5 = 1 / f^(z^2)
	//printf("13:");printf_mul_counts();

	y6 = ext12_frobenius_p(t3, p);		//y6 = ( f^(-z^3) )^p
	y6 = mul_ext12(y6, t3, p);			//y6 = f^(-z^3) * ( f^(-z^3) )^p
	//printf("14:");printf_mul_counts();


	t0 = square_ext12(y6, p);			//t0 = y6^2
	//printf("15:");printf_mul_counts();
	t0 = mul_ext12(t0, y4, p);			//t0 = y4 * y6^2
	//printf("16:");printf_mul_counts();
	t0 = mul_ext12(t0, y5, p);			//t0 = y4 * y5 * y6^2
	t = mul_ext12(y3, y5, p);			//t = y3 * y5
	t = mul_ext12(t, t0, p);			//t = y3 * y4 * y5^2 * y6^2
	t0 = mul_ext12(t0, y2, p);			//t0 = y2 * y4 * y5 * y6^2
	t = square_ext12(t, p);				//t = y3^2 * y4^2 * y5^4 * y6^4
	t = mul_ext12(t, t0, p);			//t = y2 * y3 ^2 * y4^3 * y5^5 * y6^6
	t = square_ext12(t, p);				//t = y2^2 * y3 ^4 * y4^6 * y5^10 * y6^12
	t0 = mul_ext12(t, y1, p);			//t0 = y1 * y2^2 * y3 ^4 * y4^6 * y5^10 * y6^12
	t = mul_ext12(t, y0, p);			//t = y0 * y2^2 * y3 ^4 * y4^6 * y5^10 * y6^12
	t0 = square_ext12(t0, p);			//t0 = y1^2 * y2^4 * y3 ^8 * y4^12 * y5^20 * y6^24
	res = mul_ext12(t0, t, p);			//res = y0 * y1^2 * y2^6 * y3 ^12 * y4^18 * y5^30 * y6^36
	//printf("17:");printf_mul_counts();

	return res;
}













