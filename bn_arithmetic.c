#include "bn_arithmetic.h"

//(p-1)/6
un_bn p_sub1_div6 = {{0x3, 0, 0, 0x4680, 0x3, 0, 0x8000, 0x6585, 0x1, 0, 0x6240, 0x49B3, 0, 0x6000, 0x90C0, 0x630}, 0};



//Copy the values of orig_op to assig_op.
void un_bn_cpy(un_bn *assig_op, un_bn orig_op){
	memset(assig_op, 0, sizeof(un_bn));

	int i;
	for(i=0; i<SIZE; i++)
		assig_op->bn[i] = orig_op.bn[i];
}

/*
 * return BIGGER if a>b
 * return LESS if a<b
 * else return EQUAL
 */
compareres un_bn_cmp(un_bn a, un_bn b){
	int i;

	for(i=SIZE-1; i>=0; i--){
		if(a.bn[i]>b.bn[i]){
			return BIGGER;
		}else if(a.bn[i]<b.bn[i])
			return LESS;
	}
	return EQUAL;
}

compareres un_bn_equal1(un_bn a){
	int i;

	for(i=SIZE-1; i>=1; i--){
		if(a.bn[i] != 0)
			return NO;
	}
	if(a.bn[0] == 1)
		return YES;
	else
		return NO;
}

compareres un_bn_even(un_bn a){
	unint t;

	t = a.bn[0];

	if(t%2 == 0)
		return YES;
	else
		return NO;
}

un_bn un_bn_negation(un_bn op){
	int i;
	un_bn res;

	memset(&res, 0, sizeof(un_bn));

	for(i=0; i<SIZE; i++)
		res.bn[i] = (op.bn[i]^0x0ffff)&0x0ffff;

	return res;
}

void un_bn_binary(unint res[], unint op){
	int i;

	memset(res, 0, sizeof(unint));

	for(i=0; i<DIGIT_SIZE; i++){
		res[i] = (op>>i)&0x1;
	}
}

//Left Shift on op with n bits.
un_bn un_bn_sll(un_bn op, int n){
	un_bn res;
	int i,q,r;
	unint t1,t2,flag;

	memset(&res, 0, sizeof(un_bn));
	flag = 0;

	q = n/DIGIT_SIZE;
	r = n%DIGIT_SIZE;

	if(q > 0){
		flag = op.bn[SIZE-q]&0x1;
		for(i=SIZE-1; i>=q; i--)
			op.bn[i] = op.bn[i-q];

		for( ; i>=0; i--)
			op.bn[i] = 0;
	}

	if(r == 0){
		un_bn_cpy(&res, op);
		res.c = flag;
		return res;
	}

	flag = ( (op.bn[SIZE-1] & 0x0ffff) >> (DIGIT_SIZE-r) ) & 0x1;

	for(i=SIZE-1; i>=1; i--){
		t1 = op.bn[i];
		t2 = op.bn[i-1];

		t1 = (t1 << r ) & 0x0ffff;
		t2 = (t2 & 0x0ffff) >> (DIGIT_SIZE-r);
		op.bn[i] = (t1)|(t2);
	}
	op.bn[0] = ((op.bn[0])<<r)&0x0ffff;

	un_bn_cpy(&res, op);
	res.c = flag;

	return res;
}

//Right Shift on op with n bits.
un_bn un_bn_srl(un_bn op, int n){      //
	int i,q,r;
	unint t1,t2;
	un_bn res;

	q = n/DIGIT_SIZE;
	r = n%DIGIT_SIZE;

	memset(&res, 0, sizeof(un_bn));

	if(q >= SIZE){
		if(q==SIZE)
			res.bn[0] = op.c;
		return res;
	}

	if(q > 0){
		for(i=0; i<(SIZE-q); i++)
			res.bn[i] = op.bn[i+q];
		res.bn[SIZE-q] = op.c;
	}

	if(r == 0)
		return res;

	for(i=0; i<SIZE-1; i++){
		t1 = op.bn[i];
		t2 = op.bn[i+1];

		t1 = ( ( t1 & 0x0ffff ) >> r );
		t2 = ( t2 << (DIGIT_SIZE-r)) & 0x0ffff;
		res.bn[i] = (t1)|(t2);
	}
	if(op.c == 1)
		op.bn[SIZE-1] = op.bn[SIZE-1] + 0x10000;
	res.bn[SIZE-1] = ( (op.bn[SIZE-1])>>r ) & 0x0ffff;

	return res;
}

un_bn_64 add_64(un_bn_64 a, un_bn_64 b, unint cin){
	un_bn_64 res;
	undlong s0,a0,a1,b0,b1,s[2];
	unint c0,c[2];

	memset(&res, 0, sizeof(un_bn_64));

	a0 = a.bn[0]+a.bn[1]*65536;
	a1 = a.bn[2]+a.bn[3]*65536;
	b0 = b.bn[0]+b.bn[1]*65536;
	b1 = b.bn[2]+b.bn[3]*65536;

	s0 = a0+b0+cin;
	c0 = 0;
	if(s0 >= 4294967296){
		s0 -= 4294967296;
		c0 = 1;
	}
	//printf("%d \n", c0);
	s[0] = a1+b1;
	c[0] = 0;
	if(s[0] >= 4294967296){
		s[0] -= 4294967296;
		c[0] = 1;
	}

	s[1] = a1+b1+1;
	c[1] = 0;
	if(s[1] >= 4294967296){
		s[1] -= 4294967296;
		c[1] = 1;
	}

	//printf("%lld %lld %lld\n", s0, s[0],s[1]);
	res.bn[0] = s0&0x0ffff;
	res.bn[1] = (s0>>16)&0x0ffff;
	res.bn[2] = s[c0]&0x0ffff;
	res.bn[3] = (s[c0]>>16)&0x0ffff;
	res.c = c[c0];

	//printf("add_64: %d %d %d %d c=%d\n", res.bn[3],res.bn[2],res.bn[1],res.bn[0], res.c);

	return res;
}

un_bn_128 add_128(un_bn_128 a, un_bn_128 b, unint cin){
	un_bn_128 res;
	un_bn_64 s0,a0,a1,b0,b1,s[2];
	int i;

	memset(&res, 0, sizeof(un_bn_128));

	for(i=0; i<4; i++){
		a0.bn[i] = a.bn[i];
		a1.bn[i] = a.bn[i+4];
	}

	for(i=0; i<4; i++){
		b0.bn[i] = b.bn[i];
		b1.bn[i] = b.bn[i+4];
	}


	s0 = add_64(a0, b0, cin);

	s[0] = add_64(a1, b1, 0);
	s[1] = add_64(a1, b1, 1);

	for(i=0; i<4; i++){
		res.bn[i] = s0.bn[i];
		res.bn[i+4] = s[s0.c].bn[i];
	}
	res.c = s[s0.c].c;

	return res;
}

un_bn add_256(un_bn a, un_bn b, unint cin){
	un_bn res;
	un_bn_128 s0,a0,a1,b0,b1,s[2];
	int i;

	memset(&res, 0, sizeof(un_bn_128));

	for(i=0; i<8; i++){
		a0.bn[i] = a.bn[i];
		a1.bn[i] = a.bn[i+8];
	}

	for(i=0; i<8; i++){
		b0.bn[i] = b.bn[i];
		b1.bn[i] = b.bn[i+8];
	}

	s0 = add_128(a0, b0, cin);

	s[0] = add_128(a1, b1, 0);
	s[1] = add_128(a1, b1, 1);

	for(i=0; i<8; i++){
		res.bn[i] = s0.bn[i];
		res.bn[i+8] = s[s0.c].bn[i];
	}
	res.c = s[s0.c].c;

	return res;
}

/*****************************************************************************************
 * 								Basic operations in	prime field			    	    	 *
 *****************************************************************************************/
//Interleaved Montgomery ladder
//return a*b mod p
un_bn imml(un_bn a, un_bn b, un_bn p){

	//mul_counts++;

	un_bn s[2],t[2],tt[2],np;
	unint median[DIGIT_SIZE];
	int i,j,bi,nbi;

	memset(&median, 0, sizeof(unint)*DIGIT_SIZE);
	np = un_bn_negation(p);

	memset(&s[1], 0, sizeof(un_bn));			//s[1] = 0
	un_bn_cpy(&s[0], a);						//s[0] = a

	for(i=SIZE-1; i>=0; i--){
		un_bn_binary(median, b.bn[i]);
		for(j=DIGIT_SIZE-1; j>=0; j--){
			bi = median[j]&0x1;
			nbi = bi==0?1:0;					//~bi

			//printf("%d \n", bi);

			t[bi] = add_256(s[bi], s[nbi], 0);	//{c0, t[bi]} = add256(s[bi], s[~bi], 0)
			tt[bi] = add_256(t[bi], np, 1);		//{cc0, tt[bi]} = add256(t[bi], ~p, 1)
			t[nbi] = un_bn_sll(s[nbi], 1);		//{c1, t[~bi]} = s[~bi] << 1
			tt[nbi] = add_256(t[nbi], np, 1);	//{cc1, tt[~bi]} = add256(t[~bi], ~p, 1)

			if(t[bi].c || tt[bi].c)				//if(c0 or cc0) then s[bi] = tt[bi] else s[bi] = t[bi]
				un_bn_cpy(&s[bi], tt[bi]);
			else
				un_bn_cpy(&s[bi], t[bi]);

			if(t[nbi].c || tt[nbi].c)			//if(c1 or cc1) then s[~bi] = tt[~bi] else s[~bi] = t[~bi]
				un_bn_cpy(&s[nbi], tt[nbi]);
			else
				un_bn_cpy(&s[nbi], t[nbi]);
		}
	}

	//printf_res(a);
	//printf(" * ");
	//printf_res(b);
	//printf(" = ");
	//printf_res(s[1]);

	return s[1];
}

// a+b mod p
un_bn add_mod(un_bn a, un_bn b, un_bn p){
	un_bn res,s1,s2,v2,w2,np;
	int c2;

	memset(&res, 0, sizeof(un_bn));
	np = un_bn_negation(p);

	un_bn_cpy(&s1, a);
	un_bn_cpy(&s2, b);

	v2 = add_256(s1, s2, 0);
	w2 = add_256(v2, np, 1);

	c2 = v2.c|w2.c;

	if(c2 == 1)
		un_bn_cpy(&res, w2);
	else
		un_bn_cpy(&res, v2);

	return res;
}

//a-b mod p
un_bn sub_mod(un_bn a, un_bn b, un_bn p){
	un_bn res,s1,s2,v2,w2,ns2;
	int c2;

	memset(&res, 0, sizeof(un_bn));
	un_bn_cpy(&s1, a);
	un_bn_cpy(&s2, b);

	ns2 = un_bn_negation(s2);

	v2 = add_256(s1, ns2, 1);
	w2 = add_256(v2, p, 0);

	c2 = v2.c;
	//printf("%d\n", c2);
	if(c2 == 0)
		un_bn_cpy(&res, w2);
	else
		un_bn_cpy(&res, v2);

	return res;
}

// -a mod p
un_bn neg_mod(un_bn a, un_bn p){
	un_bn t,res;

	memset(&t, 0, sizeof(un_bn));
	res = sub_mod(t, a, p);

	return res;
}


//The interleaved multiplication based on Montgomery ladder
//return a*b mod p (This function is the same with the preceding function named imml)
un_bn mul_prime(un_bn a, un_bn b, un_bn p){

	mul_counts++;

	un_bn s1,s2, u, v1, v2, w1, w2, np, t1, t2;
	unint c1,c2, bi, mb[DIGIT_SIZE];
	int i,j;

	np = un_bn_negation(p);

	memset(&s1, 0, sizeof(un_bn));
	un_bn_cpy(&s2, a);

	for(i=SIZE-1; i>=0; i--){
		un_bn_binary(mb, b.bn[i]);

		for(j=DIGIT_SIZE-1; j>=0; j--){
			bi = mb[j]&0x1;

			//printf("bi = %d\n", bi);
			if(bi == 1)
				un_bn_cpy(&u, s2);
			else
				un_bn_cpy(&u, s1);

			v1 = un_bn_sll(u, 1);
			v2 = add_256(s1, s2, 0);

			w1 = add_256(v1, np, 1);
			w2 = add_256(v2, np, 1);

			c1 = (v1.c|w1.c)&0x1;
			c2 = (v2.c|w2.c)&0x1;

			if(c1 == 1)
				un_bn_cpy(&t1, w1);
			else
				un_bn_cpy(&t1, v1);

			if(c2 == 1)
				un_bn_cpy(&t2, w2);
			else
				un_bn_cpy(&t2, v2);

			if(bi == 0 )
				un_bn_cpy(&s1, t1);
			else
				un_bn_cpy(&s1, t2);

			if(bi == 0)
				un_bn_cpy(&s2, t2);
			else
				un_bn_cpy(&s2, t1);
		}
	}
	return s1;
}

/*
 * return (a^-1) mod p
 */
un_bn inverse_mod(un_bn a, un_bn p){
	un_bn u,v,x1,x2;

	inv_counts++;

	un_bn_cpy(&u, a);				//u = a
	un_bn_cpy(&v, p);				//v = p
	memset(&x1, 0, sizeof(un_bn));
	x1.bn[0] = 1;					//x1 = 1
	memset(&x2, 0, sizeof(un_bn));	//x2 = 0

	while( un_bn_equal1(u)==NO && un_bn_equal1(v)==NO ){
		/*printf("u = ");
		printf_res(u);
		printf("v = ");
		printf_res(v);*/
		while(un_bn_even(u)==YES){
			u = un_bn_srl(u,1);
			if(un_bn_even(x1) == YES){
				x1 = un_bn_srl(x1,1);
			}else{
				x1 = add_256(x1, p, 0);
				x1 = un_bn_srl(x1,1);
			}
		}
		while (un_bn_even(v) == YES) {
			v = un_bn_srl(v, 1);
			if (un_bn_even(x2) == YES) {
				x2 = un_bn_srl(x2, 1);
			} else {
				x2 = add_256(x2, p, 0);
				x2 = un_bn_srl(x2, 1);
			}
		}
		if(un_bn_cmp(u, v) != LESS){		// u >= v
			u = sub_mod(u, v, p);
			x1 = sub_mod(x1, x2, p);
		}else{
			v = sub_mod(v, u, p);
			x2 = sub_mod(x2, x1, p);
		}
	}
	if(un_bn_equal1(u)==YES)
		return x1;
	else
		return x2;
}


void printf_mul_counts(){
	printf("Prime field multiplication counts: %lld\n", mul_counts);
}

void printf_inv_counts(){
	printf("Prime field inversion counts: %lld\n", inv_counts);
}

/*****************************************************************************************
 * 								Basic operations in	GF(p^2)						 		 *
 *****************************************************************************************/
//copy the values of orign into assig
void ext2_cop(un_bn_ext2 *assig, un_bn_ext2 orign){
	un_bn_cpy(&(assig->bn0), orign.bn0);
	un_bn_cpy(&(assig->bn1), orign.bn1);
}

/*
 * Compute the negative of a.
 * return -a mod p
 */
un_bn_ext2 neg_ext2(un_bn_ext2 a, un_bn p){
	un_bn_ext2 res,t;

	memset(&t, 0, sizeof(un_bn_ext2));
	res = sub_ext2(t, a, p);

	return res;
}

/*The addition in GF(p^2)
 *RES = A+B mod (u^2 - β)
 *A = a[0] + a[1]u; B = b[0] + b[1]u;
 */
un_bn_ext2 add_ext2(un_bn_ext2 a, un_bn_ext2 b, un_bn p){
	un_bn_ext2 res;

	res.bn0 = add_mod(a.bn0, b.bn0, p);
	res.bn1 = add_mod(a.bn1, b.bn1, p);

	return res;
}


un_bn_ext2 ext2_add_prime(un_bn_ext2 a, un_bn b, un_bn p){
	un_bn_ext2 res;

	res.bn0 = add_mod(a.bn0, b, p);
	un_bn_cpy(&(res.bn1), a.bn1);

	return res;
}

/*The subtraction in GF(p^2)
 *RES = A-B mod (u^2 - β)
 *A = a[0]+a[1]u; B = b[0]+b[1]u;
 */
un_bn_ext2 sub_ext2(un_bn_ext2 a, un_bn_ext2 b, un_bn p){
	un_bn_ext2 res;

	res.bn0 = sub_mod(a.bn0, b.bn0, p);
	res.bn1 = sub_mod(a.bn1, b.bn1, p);

	return res;
}

/*
 * The multiplication in GF(p^2)
 *RES = A*B mod (u^2 - β) where β = -1
 *A = a[0] + a[1]u; B = b[0] + b[1]u;
 *return RES = res[0] + res[1]u
 *
 *Karatsuba method, 3M+5A+B
 */
un_bn_ext2 mul_ext2(un_bn_ext2 a, un_bn_ext2 b, un_bn p){
	un_bn t1,t2,t3,t4;
	un_bn_ext2 res;

	t3 = add_mod(a.bn0, a.bn1, p);	//t3 = a0+a1;
	t4 = add_mod(b.bn0, b.bn1, p);	//t4 = b0+b1;

	t1 = mul_prime(a.bn0, b.bn0, p);		//t1 = a0*b0;
	t2 = mul_prime(a.bn1, b.bn1, p);		//t2 = a1*b1;
	t3 = mul_prime(t3, t4, p);			//t3 = t3*t4; (t3 = (a0+a1)*(b0+b1))

	t4 = add_mod(t1, t2, p);		//t4 = t1+t2; (t4 = a0*b0+a1*b1)
	t1 = sub_mod(t1, t2, p);		//t1 = t1-t2; (t1 = a0*b0-a1*b1)

	//t1 = sub_mod(t1, t2, p);		//t1 = t1-t2; (t1 = a0*b0-2*a1*b1)
	t2 = sub_mod(t3, t4, p);		//t2 = t3-t4; (t2 = (a0+a1)*(b0+b1)-a0*b0-a1*b1)

	//printf_res(a.bn0);printf_res(a.bn1);
	//printf(" * \n");printf_res(b.bn0);printf_res(b.bn1);
	//printf(" = ");printf_res(t1);printf_res(t2);

	un_bn_cpy(&res.bn0, t1);
	un_bn_cpy(&res.bn1, t2);

	return res;
}

/*
 * a∈GF(p^2), b∈GF(p)
 *
 * 	return a*b mod (u^2 - β)
 */
un_bn_ext2 ext2_mul_prime(un_bn_ext2 a, un_bn b, un_bn p){
	un_bn_ext2 res;

	res.bn0 = mul_prime(a.bn0, b, p);
	res.bn1 = mul_prime(a.bn1, b, p);

	return res;
}

/*
 * return (a^2) mod f(x) where a∈GF(p^2)
 * Complex squaring method which only requires 2 multiplications.
 */
un_bn_ext2	square_ext2(un_bn_ext2 a, un_bn p){
	un_bn_ext2 res;
	un_bn v0,t0,t1;

	v0 = mul_prime(a.bn0, a.bn1, p);		//v0 = a0*a1

	//t0 = add_mod(a.bn1, a.bn1, p);		//t0 = 2*a1
	t1 = sub_mod(a.bn0, a.bn1, p);			//t1 = a0-a1 (a0+β*a1)
	t0 = add_mod(a.bn0, a.bn1, p);			//t0 = a0+a1
	res.bn0 = mul_prime(t0, t1, p);			//(a0+a1)*(a0+βa1)
	t1 = sub_mod(res.bn0, v0, p);			//(a0+a1)*(a0+β*a1)-v0

	res.bn1 = add_mod(v0, v0, p);			//2*v0
	//res.bn0 = add_mod(t1, v0, p);			//(a0+a1)*(a0+β*a1)-v0-β*v0

	return res;
}

/*
 * return (a^-1) mod f(x) where a∈GF(p^2)
 */
un_bn_ext2 inverse_ext2(un_bn_ext2 a, un_bn p){
	un_bn_ext2 res;
	un_bn t1,t2,t3;

	t1 = mul_prime(a.bn0, a.bn0, p);
	t2 = mul_prime(a.bn1, a.bn1, p);

	//t3 = add_mod(t2, t2, p);
	t3 = add_mod(t1, t2, p);			//t3 = t1+t2

	t1 = inverse_mod(t3, p);			//t1 = t3^-1
	res.bn0 = mul_prime(a.bn0, t1, p);	//res0 = a0/t3
	t2 = neg_mod(a.bn1, p);				//t2 = -a1
	res.bn1 = mul_prime(t2, t1, p);		//res1 = -a1/t3

	return res;
}

/*
 * f∈GF(p^2)
 *
 *     return f^a
 */
un_bn_ext2 ext2_power_ext(un_bn_ext2 f, un_bn a, un_bn p){
	un_bn_ext2 res;
	unint m[DIGIT_SIZE];
	int i,j,flag,ai;
	un_bn_ext2 q1,q2,t;

	flag = 0;

	for(i=SIZE-1; i>=0; i--){
		un_bn_binary(m, a.bn[i]);
		for(j=DIGIT_SIZE-1; j>=0; j--){
			ai = m[j]&0x1;
			//printf("%d ", ai);
			if(ai == 0 && flag == 0)
				continue;
			else{
				if(flag == 0){				//The first time to have '1'
					ext2_cop(&q1, f);
					t = square_ext2(f, p);
					//t = mul_ext2(f, f, p);
					ext2_cop(&q2, t);
					flag = 1;
					//printf("1--\n");
				}else{

					if(ai == 1){
						q1 = mul_ext2(q1, q2, p);
						q2 = square_ext2(q2, p);
						//q2 = mul_ext2(q2, q2, p);
						//printf("2--\n");
					}else{
						q2 = mul_ext2(q1, q2, p);
						q1 = square_ext2(q1, p);
						//q1 = mul_ext2(q1, q1, p);
						//printf("3--\n");
					}
				}
			}
		}
	}
	ext2_cop(&res, q1);

	return res;
}

/*
 * return ξ^p_sub1_div6
 */
un_bn_ext2 compute_xipow(un_bn p){
	un_bn_ext2 res,a;
	un_bn t;

	memset(&a, 0, sizeof(un_bn_ext2));
	memset(&t, 0, sizeof(un_bn));
	t.bn[0] = 1;
	a.bn0 = neg_mod(t, p);
	a.bn1 = neg_mod(t, p);

	//printf_res(a.bn0);
	//printf_res(a.bn1);
	//printf("p = ");printf_res(p);
	//printf("p_sub1_div6 = ");printf_res(p_sub1_div6);

	res = ext2_power_ext(a, p_sub1_div6, p);

	/*res = square_ext2(a, p);

	res = mul_ext2(a, a, p);
	printf_res(a.bn0);
	printf_res(a.bn1);*/

	return res;
}
/*****************************************************************************************
 * 								Basic operations in	GF(p^6)								 *
 *****************************************************************************************/
void ext6_cop(un_bn_ext6 *assig, un_bn_ext6 orign){
	ext2_cop(&(assig->V0), orign.V0);
	ext2_cop(&(assig->V1), orign.V1);
	ext2_cop(&(assig->V2), orign.V2);
}

/*
 * Compute the negative of a.
 * return -a mod p
 */
un_bn_ext6 neg_ext6(un_bn_ext6 a, un_bn p){
	un_bn_ext6 res,t;

	memset(&t, 0, sizeof(un_bn_ext6));
	res = sub_ext6(t, a, p);

	return res;
}

/*
 * The addition in GF(p^6)
 * return a+b
 */
un_bn_ext6 add_ext6(un_bn_ext6 a, un_bn_ext6 b, un_bn p){
	un_bn_ext6 res;

	res.V0 = add_ext2(a.V0, b.V0, p);
	res.V1 = add_ext2(a.V1, b.V1, p);
	res.V2 = add_ext2(a.V2, b.V2, p);

	return res;
}

/*
 * The subtraction in GF(p^6)
 * return a-b
 */
un_bn_ext6 sub_ext6(un_bn_ext6 a, un_bn_ext6 b, un_bn p){
	un_bn_ext6 res;

	res.V0 = sub_ext2(a.V0, b.V0, p);
	res.V1 = sub_ext2(a.V1, b.V1, p);
	res.V2 = sub_ext2(a.V2, b.V2, p);

	return res;
}

/*
 * The multiplication between the non-residue ξ (ξ=-1-u) and a which is in GF(p^2)
 *
 */
un_bn_ext2 nr_ext6(un_bn_ext2 a, un_bn p){
	un_bn_ext2 res;
	un_bn t0; //t1;

	/*//ξ=0+u
	memset(&t0, 0, sizeof(un_bn));
	t1 = mul_sn(a.bn0, 2, p);

	res.bn0 = sub_mod(t0, t1, p);
	un_bn_cpy(&res.bn1, a.bn0);*/

	//ξ=-1-u
	t0 = add_mod(a.bn0, a.bn1, p);
	res.bn1 = neg_mod(t0, p);		//res1 = -(a0+a1)
	//t1 = add_mod(a.bn1, a.bn1, p);
	res.bn0 = sub_mod(a.bn1, a.bn0, p);//res0 = -ξa1-a0

	return res;
}

/*
 * The multiplication in GF(p^6)
 * GF(P^2)[V]/(V^3-ξ)  where ξ∈GF(p^2)
 */
un_bn_ext6 mul_ext6(un_bn_ext6 a, un_bn_ext6 b, un_bn p){
	un_bn_ext2 v0, v1, v2;
	un_bn_ext2 t0,t1,t2;
	un_bn_ext6 res;

	v0 = mul_ext2(a.V0, b.V0, p);			//v0 = a0*b0
	v1 = mul_ext2(a.V1, b.V1, p);			//v1 = a1*b1
	v2 = mul_ext2(a.V2, b.V2, p);			//v2 = a2*b2

	t0 = add_ext2(a.V1, a.V2, p);			//t0 = a1+a2
	t1 = add_ext2(b.V1, b.V2, p);			//t1 = b1+b2
	t2 = mul_ext2(t0, t1, p);				//t2 = (a1+a2)*(b1+b2)
	t2 = sub_ext2(t2, v1, p);				//t2 = (a1+a2)*(b1+b2)-v1
	t2 = sub_ext2(t2, v2, p);				//t2 = (a1+a2)*(b1+b2)-v1-v2
	t2 = nr_ext6(t2, p);					//t2 = ξ*( (a1+a2)*(b1+b2)-v1-v2 )
	res.V0 = add_ext2(v0, t2, p);			//res0 = v0 + ξ*( (a1+a2)*(b1+b2)-v1-v2 )

	t0 = add_ext2(a.V0, a.V1, p);			//t0 = a0+a1
	t1 = add_ext2(b.V0, b.V1, p);			//t1 = b0+b1
	t2 = mul_ext2(t0, t1, p);				//t2 = (a0+a1)*(b0+b1)
	t0 = nr_ext6(v2, p);					//t0 = ξ*v2
	t0 = sub_ext2(t0, v1, p);				//t0 = ξ*v2-v1
	t1 = sub_ext2(t2, v0, p);				//t1 = (a0+a1)*(b0+b1)-v0
	res.V1 = add_ext2(t1, t0, p);			//res1 = (a0+a1)*(b0+b1)-v0+ξ*v2-v1

	t0 = add_ext2(a.V0, a.V2, p);			//t0 = a0+a2
	t1 = add_ext2(b.V0, b.V2, p);			//t1 = b0+b2
	t2 = mul_ext2(t0, t1, p);				//t2 = (a0+a2)*(b0+b2)
	t0 = sub_ext2(v1, v2, p);				//t0 = v1-v2
	t1 = sub_ext2(t2, v0, p);				//t1 = (a0+a2)*(b0+b2)-v0
	res.V2 = add_ext2(t0, t1, p);			//res2 = (a0+a2)*(b0+b2)-v0+v1-v2

	return res;
}

/*
 * a∈GF(p^6), l0∈GF(p), l1∈GF(p^2)
 *return
 *		a*(l0 + l1*V)
 */
un_bn_ext6 ext6_mul_sp1(un_bn_ext6 a, un_bn l0, un_bn_ext2 l1,un_bn p){
	un_bn_ext2 v0,v1;
	un_bn_ext2 t0,t1,t2;
	un_bn_ext6 res;

	v0 = ext2_mul_prime(a.V0, l0, p);
	v1 = mul_ext2(a.V1, l1, p);

	t0 = add_ext2(a.V1, a.V2, p);	//t0 = a1+a2
	t1 = mul_ext2(t0, l1, p);		//t1 = (a1+a2)*l1
	t0 = sub_ext2(t1, v1, p);		//t0 = (a1+a2)*l1-v1
	t1 = nr_ext6(t0, p);			//t1 = ξ*( (a1+a2)*l1-v1 )
	res.V0 = add_ext2(v0, t1, p);	//res0 = v0+ξ*( (a1+a2)*l1-v1 )

	t0 = add_ext2(a.V0, a.V1, p);	//t0 = a0+a1
	t1 = ext2_add_prime(l1, l0, p);	//t1 = l1+l0
	t2 = mul_ext2(t0, t1, p);		//t2 = (a0+a1)*(l1+l0)
	t0 = sub_ext2(t2, v0, p);		//t0 = ((a0+a1)*(l1+l0)) - v0
	res.V1 = sub_ext2(t0, v1, p);	//res1 = ((a0+a1)*(l1+l0)) - v0 - v1

	t1 = add_ext2(a.V0, a.V2, p);
	t0 = ext2_mul_prime(t1, l0, p);
	t2 = sub_ext2(t0, v0, p);
	res.V2 = add_ext2(t2, v1, p);	//res1 = (a0+a2)*l0-v0+v1

	return res;
}

/*
 * a∈GF(p^6), l2∈GF(p^2)
 *return
 *		a*(l2*V)
 */
un_bn_ext6 ext6_mul_sp2(un_bn_ext6 a, un_bn_ext2 l2,un_bn p){
	un_bn_ext2 v1;
	un_bn_ext2 t0,t1,t2;
	un_bn_ext6 res;

	v1 = mul_ext2(a.V1, l2, p);

	t0 = add_ext2(a.V1, a.V2, p);	//t0 = a1+a2
	t1 = mul_ext2(t0, l2, p);		//t1 = (a1+a2)*l2
	t2 = sub_ext2(t1, v1, p);		//t2 = (a1+a2)*l2-v1
	res.V0 = nr_ext6(t2, p);		//res0 = ξ*( (a1+a2)*l2-v1 )

	t0 = add_ext2(a.V0, a.V1, p);	//t0 = a0+a1
	t2 = mul_ext2(t0, l2, p);		//t2 = (a0+a1)*l2
	res.V1 = sub_ext2(t2, v1, p);	//res1 = (a0+a1)*l2 - v1

	ext2_cop(&(res.V2), v1);		//res2 = v1

	return res;
}

/*
 * a^2 a∈ GF(p^6)
 *
 * Chung-Hasan SQR2 method which requires 2M+3S+10A+2B
 */
un_bn_ext6 square_ext6(un_bn_ext6 a, un_bn p){
	un_bn_ext6 res;
	un_bn_ext2 t0,t1,t2,t3,t4,t5;

	t0 = square_ext2(a.V0, p);		//a0^2

	t2 = mul_ext2(a.V0, a.V1, p);	//a0*a1
	t1 = add_ext2(t2, t2, p);		//2*a0*a1

	t3 = sub_ext2(a.V0, a.V1, p);	//a0-a1
	t4 = add_ext2(t3, a.V2, p);		//a0-a1+a2
	t2 = square_ext2(t4, p);		//(a0-a1+a2)^2

	t4 = mul_ext2(a.V1, a.V2, p);	//a1*a2
	t3 = add_ext2(t4, t4, p);		//2*a1*a2

	t4 = square_ext2(a.V2, p);		//a2^2

	t5 = nr_ext6(t3, p);			//ξ*t3
	res.V0 = add_ext2(t0, t5, p);	//t0+ξ*t3
	t5 = nr_ext6(t4, p);			//ξ*t4
	res.V1 = add_ext2(t1, t5, p);	//t1+ξ*t4
	t1 = add_ext2(t1, t2, p);		//t1+t2
	t3 = sub_ext2(t3, t0, p);		//t3-t0
	t5 = add_ext2(t1, t3, p);		//t1+t2+t3-t0
	res.V2 = sub_ext2(t5, t4, p);	//t1+t2+t3-t0-t4

	return res;
}

/*
 * return (a^-1) mod f(x) where a∈GF(p^6)
 * 	A = a0^2-ξ*a1*a2; B = ξ*a2^2-a0*a1; C = a1^2-a0*a2
 * 	F = ξ*a1*C + a0*A + ξ*a2*B
 *  (a^-1) = (A+B*V+C*V^2)/F
 */
un_bn_ext6 inverse_ext6(un_bn_ext6 a, un_bn p){
	un_bn_ext6 res;
	un_bn_ext2 t1,t2,t3,t4;
	un_bn_ext2 A,B,C,F,FI;

	t1 = square_ext2(a.V0, p);		//t1 = a0^2
	t2 = mul_ext2(a.V1, a.V2, p);	//t2 = a1*a2
	t2 = nr_ext6(t2, p);			//t2 = ξ*a1*a2
	A = sub_ext2(t1, t2, p);		//A = a0^2 - ξ*a1*a2

	t1 = square_ext2(a.V2, p);		//t1 = a2^2
	t2 = nr_ext6(t1, p);			//t2 = ξ*(a2^2)
	t1 = mul_ext2(a.V0, a.V1, p);	//t1 = a0*a1
	B = sub_ext2(t2, t1, p);		//B = ξ*(a2^2) - a0*a1

	t1 = square_ext2(a.V1, p);		//t1 = a1^2
	t2 = mul_ext2(a.V0, a.V2, p);	//t2 = a0*a2
	C = sub_ext2(t1, t2, p);		//C = a1^2-a0*a2

	//ξ*a1*C + a0*A + ξ*a2*B
	t1 = mul_ext2(a.V1, C, p);
	t2 = nr_ext6(t1, p);			//t2 = ξ*a1*C
	t1 = mul_ext2(a.V2, B, p);
	t3 = nr_ext6(t1, p);			//t3 = ξ*a2*B
	t1 = mul_ext2(a.V0, A, p);		//t1 = a0*A
	t4 = add_ext2(t1, t2, p);		//t4 = a0*A+ξ*a1*C
	F  = add_ext2(t4, t3, p);		//F  = a0*A+ξ*a1*C+ξ*a2*B

	FI = inverse_ext2(F, p);		//FI = (a0*A+ξ*a1*C+ξ*a2*B)^-1

	res.V0 = mul_ext2(A, FI, p);
	res.V1 = mul_ext2(B, FI, p);
	res.V2 = mul_ext2(C, FI, p);

	return res;
}
/*****************************************************************************************
 * 							 	Basic operations in	GF(p^12)							 *
 *****************************************************************************************/
void ext12_cop(un_bn_ext12_t *assig, un_bn_ext12_t orign){
	ext6_cop(&(assig->W0), orign.W0);
	ext6_cop(&(assig->W1), orign.W1);
}

/*
 * Compute the negative of a.
 * return -a mod p
 */
un_bn_ext12_t neg_ext12(un_bn_ext12_t a, un_bn p){
	un_bn_ext12_t res;
	un_bn_ext6 t;

	memset(&t, 0, sizeof(un_bn_ext6));
	res.W0 = sub_ext6(t, a.W0, p);
	res.W1 = sub_ext6(t, a.W1, p);

	return res;
}

/*
 * The multiplication between the non-residue γ and a which is in GF(p^6)
 * γ = 0 + V + 0*V^2
 */
un_bn_ext6 nr_ext12(un_bn_ext6 a, un_bn p){
	un_bn_ext6 res;
	//un_bn_ext2 t0;

	res.V0 = nr_ext6(a.V2, p);			//res0 = γ*a2

	un_bn_cpy(&res.V1.bn0, a.V0.bn0);
	un_bn_cpy(&res.V1.bn1, a.V0.bn1);	//res1 = a0

	//t0 = add_ext2(a.V0, a.V1, p);
	ext2_cop(&res.V2, a.V1);			//res2 = a1

	return res;
}

/*
 * The square of a in GF(p^12)
 * GF(P^6)[W]/(W^2-γ)  where ∈GF(p^6), γ=0+V+0*V^2
 *
 * Complex squaring method.
 */
un_bn_ext12_t square_ext12(un_bn_ext12_t a, un_bn p){
	un_bn_ext6 f0,f1,v,t0,t1;
	un_bn_ext12_t res;

	ext6_cop(&f0, a.W0);
	ext6_cop(&f1, a.W1);

	v = mul_ext6(f0, f1, p);		//v = f0*f1

	t0 = add_ext6(f0, f1, p);		//t0 = f0+f1
	t1 = nr_ext12(f1, p);			//t1 = γ*f1
	t1 = add_ext6(f0, t1, p);		//t1 = f0+γ*f1
	t0 = mul_ext6(t0, t1, p);		//t0 = (f0+f1)*(f0+γ*f1)
	t0 = sub_ext6(t0, v, p);		//t0 = (f0+f1)*(f0+γ*f1)-v
	t1 = nr_ext12(v, p);			//t1 = γ*v
	res.W0 = sub_ext6(t0, t1, p);	//t0 = (f0+f1)*(f0+γ*f1)-v-γ*v

	res.W1 = add_ext6(v, v, p);		//t1 = 2*v

	return res;
}

/*
 * The multiplication of a in GF(p^12)
 * GF(P^6)[W]/(W^2-γ)  where γ∈GF(p^6), γ=0+V+0*V^2
 *
 * Karatsuba method
 */
un_bn_ext12_t mul_ext12(un_bn_ext12_t a, un_bn_ext12_t b, un_bn p){
	un_bn_ext6 a0,a1,b0,b1,v0,v1,t0,t1,t2;
	un_bn_ext12_t res;

	ext6_cop(&a0, a.W0);
	ext6_cop(&a1, a.W1);
	ext6_cop(&b0, b.W0);
	ext6_cop(&b1, b.W1);

	v0 = mul_ext6(a0, b0, p);		//v = a0*b0
	v1 = mul_ext6(a1, b1, p);		//v = a1*b1

	t0 = nr_ext12(v1, p);			//t0 = γ*v1
	res.W0 = add_ext6(v0, t0, p);	//t0 = v0+γ*v1

	t1 = add_ext6(a0, a1, p);		//t1 = a0+a1
	t2 = add_ext6(b0, b1, p);		//t2 = b0+b1
	t1 = mul_ext6(t1, t2, p);		//t1 = (a0+a1)*(b0+b1)
	t1 = sub_ext6(t1, v0, p);		//t1 = (a0+a1)*(b0+b1)-v0
	res.W1 = sub_ext6(t1, v1, p);	//t1 = (f0+f1)*(f0+γ*f1)-v0-v1

	return res;
}

/*
 * return (a^-1) mod f(x) where a∈GF(p^12)
 *
 * (a^-1) = (a0-a1*W)/(a0^2-γ*a1^2)
 */
un_bn_ext12_t inverse_ext12(un_bn_ext12_t a, un_bn p){
	un_bn_ext12_t res;
	un_bn_ext6 t1, t2, t3, t4;

	t1 = square_ext6(a.W1, p);
	t2 = nr_ext12(t1, p);			//t2 = γ*(a1^2)
	t1 = square_ext6(a.W0, p);		//t1 = a0^2

	t3 = sub_ext6(t1, t2, p);		//t3 = a0^2-γ*(a1^2)
	t4 = inverse_ext6(t3, p);		//t4 = (a0^2-γ*(a1^2))^-1



	//printf("(test---");print_ext6(t4);

	res.W0 = mul_ext6(a.W0, t4, p);	//res0 = a0/(a0^2-γ*(a1^2))
	t1 = neg_ext6(a.W1, p);
	res.W1 = mul_ext6(t1, t4, p);	//res1 = (-a1)/(a0^2-γ*(a1^2))

	return res;
}



/*
 * a∈GF(p^12), l0∈GF(p), l1∈GF(p^2), l2∈GF(p^2)
 *return
 *		a*((l0 + l1*V) + (l2*V)*W)
 */
un_bn_ext12_t ext12_mul_lineresult(un_bn_ext12_t a, line_result l, un_bn p){
	/*un_bn_ext12_t res, t;

	memset(&t, 0, sizeof(un_bn_ext12_t));

	ext2_cop(&t.W0.V0, l.l0);
	ext2_cop(&t.W0.V1, l.l1);
	ext2_cop(&t.W1.V1, l.l2);

	res = mul_ext12(t, a, p);

	return res;*/

	un_bn l0;
	un_bn_ext2 l1,l2;
	un_bn_ext6 a0, a1, v0, v1, t0, t1, t2;
	un_bn_ext12_t res;

	ext6_cop(&a0, a.W0);
	ext6_cop(&a1, a.W1);
	//ext6_cop(&b0, b.W0);
	//ext6_cop(&b1, b.W1);
	un_bn_cpy(&l0, l.l0);
	ext2_cop(&l1, l.l1);
	ext2_cop(&l2, l.l2);

	//v0 = mul_ext6(a0, b0, p);			//v = a0*b0
	//v1 = mul_ext6(a1, b1, p);			//v = a1*b1
	v0 = ext6_mul_sp1(a0, l0, l1, p);	//a0*(l0 + l1*V)
	v1 = ext6_mul_sp2(a1, l2, p);		//a1*(l2*V)

	t0 = nr_ext12(v1, p);				//t0 = γ*v1
	res.W0 = add_ext6(v0, t0, p);		//t0 = v0+γ*v1

	t1 = add_ext6(a0, a1, p);			//t1 = a0+a1
	l1 = add_ext2(l1, l2, p);			//l1 = l1+l2
	t2 = ext6_mul_sp1(t1, l0, l1, p);	//t2 = (a0+a1)*( (l0 + l1*V)+(l2*V) )
	t1 = sub_ext6(t2, v0, p);			//t1 = (a0+a1)*( (l0 + l1*V)+(l2*V) ) - v0
	res.W1 = sub_ext6(t1, v1, p);		//t1 = (f0+f1)*(f0+γ*f1)-v0-v1

	return res;
}

/*
 * a∈GF(p^12) -> f0+f1W where f0,f1∈GF(p^6)
 *
 *
 */
un_bn_ext12_t ext12_to_ext6(un_bn_ext12 a){
	un_bn_ext12_t f;

	un_bn_cpy(&(f.W0.V0.bn0), a.bn0);
	un_bn_cpy(&(f.W0.V0.bn1), a.bn6);
	un_bn_cpy(&(f.W0.V1.bn0), a.bn2);
	un_bn_cpy(&(f.W0.V1.bn1), a.bn8);
	un_bn_cpy(&(f.W0.V2.bn0), a.bn4);
	un_bn_cpy(&(f.W0.V2.bn1), a.bn10);

	un_bn_cpy(&(f.W1.V0.bn0), a.bn1);
	un_bn_cpy(&(f.W1.V0.bn1), a.bn7);
	un_bn_cpy(&(f.W1.V1.bn0), a.bn3);
	un_bn_cpy(&(f.W1.V1.bn1), a.bn9);
	un_bn_cpy(&(f.W1.V2.bn0), a.bn5);
	un_bn_cpy(&(f.W1.V2.bn1), a.bn11);

	return f;
}

/*
 * a+bW where a,b∈GF(p^6) -> res∈GF(p^12)
 */
un_bn_ext12 ext6_to_ext12(un_bn_ext12_t op){
	un_bn_ext12 res;

	un_bn_cpy(&res.bn0, op.W0.V0.bn0);
	un_bn_cpy(&res.bn6, op.W0.V0.bn1);
	un_bn_cpy(&res.bn2, op.W0.V1.bn0);
	un_bn_cpy(&res.bn8, op.W0.V1.bn1);
	un_bn_cpy(&res.bn4, op.W0.V2.bn0);
	un_bn_cpy(&res.bn10, op.W0.V2.bn1);

	un_bn_cpy(&res.bn1, op.W1.V0.bn0);
	un_bn_cpy(&res.bn7, op.W1.V0.bn1);
	un_bn_cpy(&res.bn3, op.W1.V1.bn0);
	un_bn_cpy(&res.bn9, op.W1.V1.bn1);
	un_bn_cpy(&res.bn5, op.W1.V2.bn0);
	un_bn_cpy(&res.bn11, op.W1.V2.bn1);

	return res;
}


void printf_res(un_bn a){
	int i;

	//for(i = 0; i<SIZE; i++)
	//	printf("0x%04X, ", a.bn[i]);

	printf("0x");
	for(i = SIZE-1; i>=0; i--)
		printf("%04X", a.bn[i]);
}

void print_ext6(un_bn_ext6 a){
	printf("( ");printf_res(a.V2.bn0);printf(" + ");printf_res(a.V2.bn1);printf("*u)*v^2\n");
	printf("( ");printf_res(a.V1.bn0);printf(" + ");printf_res(a.V1.bn1);printf("*u)*v\n");
	printf_res(a.V0.bn0);printf(" + ");printf_res(a.V0.bn1);printf("*u\n");
}


