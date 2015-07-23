/*
 * bn_curve_arithmetic.c
 *
 *  Created on: 2013-3-28
 *      Author: Suliangjian
 */

#include "bn_curve_arithmetic.h"


/*
 * Point double on BN curves in Jacobian coordinates.
 * T = (x, y, z)
 * 2T = (x3, y3, z3) where	x3 = 9*x^4-8*x*y^2,
 * 							y3 = (3x^2)(4x*y^2-x3)-8*y^4,
 * 							z3 = 2*y*z
 */
point_jc point_double(point_jc T, un_bn p){
	un_bn t0, t1, t2, t3, t4, t5, t6, t7, t8;
	point_jc res;

	//step 1
	t0 = mul_prime(T.x, T.x, p);	//t0 = x^2
	t1 = mul_prime(T.y, T.y, p);	//t1 = y^2
	t2 = mul_prime(T.y, T.z, p);	//t2 = y*z

	//step 2
	t3 = mul_prime(t0, t0, p);		//t3 = x^4
	t4 = mul_prime(T.x, t1, p);		//t4 = x*y^2
	t5 = mul_prime(t1, t1, p);		//t1 = y^4

	//step 3
	t4 = add_mod(t4, t4, p);		//t4 = 2*x*y^2
	t6 = add_mod(t3, t3, p);		//t6 = 2*x^4
	res.z = add_mod(t2, t2, p);		//z3 = 2*y*z

	//step 4
	t4 = add_mod(t4, t4, p);			//t4 = 4*x*y^2
	t6 = add_mod(t6, t6, p);			//t6 = 4*x^4
	t5 = add_mod(t5, t5, p);			//t5 = 2*y^4

	//step 5
	t8 = add_mod(t4, t4, p);			//t8 = 8*x*y^2
	t7 = add_mod(t0, t0, p);			//t7 = 2*x^2
	t6 = add_mod(t6, t6, p);			//t6 = 8*x^4

	//step 6
	t3 = add_mod(t3, t6, p);		//t3 = 9*x^4
	t5 = add_mod(t5, t5, p);			//t5 = 4*y^4

	//step 7
	res.x = sub_mod(t3, t8, p);		//x3 = 9*x^4-8*x*y^2
	t5 = add_mod(t5, t5, p);		//t5 = 8*y^4

	//step 8
	t3 = sub_mod(t4, res.x, p);		//t3 = 4*x*y^2 - x3
	t7 = add_mod(t7, t0, p);		//t7 = 3*x^2

	//step 9
	t7 = mul_prime(t7, t3, p);		//t7 = (3*x^2)*(4*x*y^2-x3)
	t4 = mul_prime(T.z, T.z, p);	//t4 = z^2
	t2 = mul_prime(T.x, t0, p);		//t2 = x^3

	//step 10
	res.y = sub_mod(t7, t5, p);		//y3 = (3*x^2)*(4*x*y^2-x3)-8*y^4

	return res;
}

/*
 *Point addition on BN curve in mixed Jacobian-affine coordinates.
 *T = (x1, y1, z1), P = (x2, y2), T+P = (x3, y3, z3)
 *return x3 = (y2*z1^3 - y1)^2 - ((x2*z1^2 - x1)^2)*(x1+x2*z1^2)
 *		 y3 = (y2*z1^3 - y1)*(x1*(x2*z1^2 - x1)^2 - x3)-y1*(x2*z1^2 - x1)^3
 *		 z3 = z1*(x2*z1^2 - x1)
 */
point_jc point_addition(point_jc T, point_ac P, un_bn p){
	point_jc res;
	un_bn t0, t1, t2, t3, t4, t5, t6;

	t0 = mul_prime(P.y, T.z, p);		//t0 = y2*z1
	t1 = mul_prime(T.z, T.z, p);		//t1 = z1^2

	t0 = mul_prime(t1, t0, p);			//t0 = y2^z1^3
	t1 = mul_prime(t1, P.x, p);			//t1 = x2*z1^2

	t4 = add_mod(t1, T.x, p);			//t4 = x2*z1^2 + x1
	t0 = sub_mod(t0, T.y, p);			//t0 = y2^z1^3 - y1
	t5 = sub_mod(t1, T.x, p);			//t5 = x2*z1^2-x1

	t3 = mul_prime(t0, t0, p);			//t3 = (y2^z1^3 - y1)^2
	res.z = mul_prime(t5, T.z, p);		//z3 = (x2*z1^2-x1)*z1
	t6 = mul_prime(t5, t5, p);			//t6 = (x2*z1^2-x1)^2

	t2 = mul_prime(t6, T.x, p);			//t2 = x1*(x2*z1^2-x1)^2
	t4 = mul_prime(t4, t6, p);			//t4 = (x2*z1^2 + x1)*(x2*z1^2-x1)^2
	t5 = mul_prime(t5, t6, p);			//t5 = (x2*z1^2-x1)^3

	res.x = sub_mod(t3, t4, p);			//x3 = (y2^z1^3 - y1)^2 - (x2*z1^2 + x1)*(x2*z1^2-x1)^2

	t2 = sub_mod(t2, res.x, p);			//t2 = x1*(x2*z1^2-x1)^2 - (y2^z1^3 - y1)^2 + (x2*z1^2 + x1)*(x2*z1^2-x1)^2

	t2 = mul_prime(t2, t0, p);			//t2 = (x1*(x2*z1^2-x1)^2 - x3)*(y2^z1^3 - y1)
	t5 = mul_prime(t5, T.y, p);			//t5 = y1*(x2*z1^2-x1)*(x2*z1^2-x1)^2

	res.y = sub_mod(t2, t5,  p);		//y3 = (x1*(x2*z1^2-x1)^2 - x3)*(y2^z1^3 - y1) - y1*(x2*z1^2-x1)^3

	return res;
}


/*
 * The computation of tangent line at point T (l(Q))
 * Q(x, y) where x,y∈GF(p^2), T(X, Y, Z) where X, Y, Z∈GF(p)
 * l(x, y) = 3*X^3 - 2*Y^2 - 3*X^2*Z^2*x + 2*YZ^3*y
 * Q(x, y) |-> Q`(x*W^2, y*W^3)
 */
line_result tangent_line(point_ac_ext2 Q, point_jc T, un_bn p){
	line_result res;
	un_bn t0,t1,t2,t3,t4,t5,t6;
	un_bn_ext2 l1,l2;
	un_bn l0;
	//un_bn_ext6 m0,m1;


	//step 1
	t0 = mul_prime(T.x, T.x, p);	//t0 = X^2
	t1 = mul_prime(T.y, T.y, p);	//t1 = Y^2
	t2 = mul_prime(T.y, T.z, p);	//t2 = Y*Z

	//step 2
	t3 = add_mod(t2, t2, p);		//t3 = 2*Y*Z
	t2 = mul_prime(T.x, t0, p);		//t2 = X^3
	t6 = mul_prime(T.z, T.z, p);	//t6 = Z^2

	//step 3
	t4 = mul_prime(t6, t0, p);		//t4 = X^2 * Z^2
	t5 = add_mod(t2, t2, p);		//t5 = 2*X^3
	t0 = mul_prime(t6, t3, p);		//t0 = 2*Y*Z^3

	//step 4
	t1 = add_mod(t1, t1, p);		//t1 = 2*Y^2
	t5 = add_mod(t2, t5, p);		//t5 = 3*X^3
	t2 = add_mod(t4, t4, p);		//t2 = 2*(X^2 * Z^2)

	//step 5
	t4 = add_mod(t4, t2, p);		//t4 = 3*(X^2 * Z^2)
	l0 = sub_mod(t5, t1, p);		//l0 = 3*X^3 - 2*Y^2

	//step 6
	l1.bn0 = mul_prime(t4, Q.x.bn0, p);		//l10 = 3*(X^2 * Z^2)*x0
	l1.bn1 = mul_prime(t4, Q.x.bn1, p);		//l11 = 3*(X^2 * Z^2)*x1
	l1.bn0 = neg_mod(l1.bn0, p);			//l10 = -3*(X^2 * Z^2)*x0
	l1.bn1 = neg_mod(l1.bn1, p);			//l11 = -3*(X^2 * Z^2)*x1

	//step 7
	l2.bn0 = mul_prime(t0, Q.y.bn0, p);		//l20 = 2*Y*Z^3 * y0
	l2.bn1 = mul_prime(t0, Q.y.bn1, p);		//l21 = 2*Y*Z^3 * y1

	//l1 = sub_etx2(l0, l1, p);
	//res = add_etx2(l1, l2, p);
	//memset(&m0, 0, sizeof(un_bn_ext6));
	//memset(&m1, 0, sizeof(un_bn_ext6));

	//un_bn_cpy(&m0.V0.bn0, l0.bn0);
	//un_bn_cpy(&m0.V1.bn0, l1.bn0);
	//un_bn_cpy(&m0.V1.bn1, l1.bn1);
	//un_bn_cpy(&m1.V1.bn0, l2.bn0);
	//un_bn_cpy(&m1.V1.bn1, l2.bn1);

	//ext6_cop(&res.W0, m0);
	//ext6_cop(&res.W1, m1);		//res = (l0+l1*V+0*V^2) + (0+l2*V+0*V^2)*W
	un_bn_cpy(&res.l0, l0);
	ext2_cop(&res.l1, l1);
	ext2_cop(&res.l2, l2);			//res = (l0 + l1*V + 0*V^2) + (0 + l2*V + 0*V^2)*W

	return res;
}

/*
 * The computation of secant line through  T and P. (l T,P(Q))
 * Q(x, y) where x,y∈GF(p^2), T(X1, Y1, Z1) is in Jacobian coordinates where X1, Y1, Z1∈GF(p)
 * and P(X2, Y2) is in affine coordinates where X2, Y2∈GF(p)
 * l T,P(x, y) = ( X2*(Y2*Z1^3 - Y1) - Y2Z3) - (Y2*Z1^3 - Y1)*x + Z3*y.  Z3 = Z1*(X2*X1^2 - X1)
 */
line_result secant_line(point_ac_ext2 Q, point_jc T, point_ac P, un_bn p){
	line_result res;
	un_bn t0,t1,t2,t3,t4;
	un_bn_ext2 l1,l2;
	un_bn l0;
	//un_bn_ext6 m0,m1;


	//step 1
	t0 = mul_prime(P.y, T.z, p);		//t0 = Y2*Z1
	t1 = mul_prime(T.z, T.z, p);		//t1 = Z1^2

	//step 2
	t0 = mul_prime(t1, t0, p);			//t0 = Y2*Z1^3
	t1 = mul_prime(t1, P.x, p);			//t1 = X2*Z1^2

	//step 3
	t2 = sub_mod(t1, T.x, p);			//t2 = X2*Z1^2 - X1
	t0 = sub_mod(t0, T.y, p);			//t0 = Y2*Z1^3 - Y1

	//step 4
	t3 = mul_prime(t2, T.z, p);			//t3 = (X2*Z1^2 - X1)*Z1 = Z3
	l1.bn0 = mul_prime(t0, Q.x.bn0, p);	//l10 = (Y2*Z1^3 - Y1)*x0
	l1.bn1 = mul_prime(t0, Q.x.bn1, p);	//l11 = (Y2*Z1^3 - Y1)*x1
	l1.bn0 = neg_mod(l1.bn0, p);		//l10 = -3*(X^2 * Z^2)*x0
	l1.bn1 = neg_mod(l1.bn1, p);		//l11 = -3*(X^2 * Z^2)*x1

	//step 5
	t4 = mul_prime(P.y, t3, p);			//t4 = Y2*Z3
	t1 = mul_prime(t0, P.x, p);			//t1 = (Y2*Z1^3 - Y1)*X2

	//step 6
	l0 = sub_mod(t1, t4, p);			//l0 = (Y2*Z1^3 - Y1)*X2 - Y2*Z3
	l2.bn0 = mul_prime(t3, Q.y.bn0, p);	//l20 = Z3*y0
	l2.bn1 = mul_prime(t3, Q.y.bn1, p);	//l21 = Z3*y1

	//l1 = sub_etx2(l0, l1, p);
	//res = add_etx2(l1, l2, p);
	//memset(&m0, 0, sizeof(un_bn_ext6));
	//memset(&m1, 0, sizeof(un_bn_ext6));

	//ext2_cop(&m0.V0, l0);
	//ext2_cop(&m0.V1, l1);
	//ext2_cop(&m1.V1, l2);

	//ext6_cop(&res.W0, m0);
	//ext6_cop(&res.W1, m1);	//res = (l0+l1*V+0*V^2) + (0+l2*V+0*V^2)*W

	un_bn_cpy(&res.l0, l0);
	ext2_cop(&res.l1, l1);
	ext2_cop(&res.l2, l2);		//res = (l0 + l1*V + 0*V^2) + (0 + l2*V + 0*V^2)*W

	return res;
}

/*
 * pd = 2*T
 * t1 = tangent_line(Q, T)
 * Total: 15M
 */
void point_double_and_tangentline(point_jc *pd,
								line_result *tl,
								point_ac_ext2 Q,
								point_jc T,
								un_bn p)
{
	un_bn t0,t1,t2,t3,t4,t5;

	t0 = mul_prime(T.x, T.x, p);			//t0 = X^2
	t1 = add_mod(t0, t0, p);				//t1 = 2*X^2
	t0 = add_mod(t1, t0, p);				//t0 = 3*X^2
	t1 = add_mod(T.y, T.y, p);				//t1 = 2*Y
	pd->z = mul_prime(t1, T.z, p);			//Z3 = 2*Y*Z

	t1 = mul_prime(t1, T.y, p);				//t1 = 2*Y^2
	t2 = mul_prime(t0, t0, p);				//t2 = 9*X^4
	t3 = add_mod(t1, t1, p);				//t3 = 4*Y^2
	t4 = mul_prime(T.x, t3, p);				//t4 = 4*X*Y^2
	t5 = add_mod(t4, t4, p);				//t5 = 8*X*Y^2
	pd->x = sub_mod(t2, t5, p);				//X3 = 9*X^4-8*X*Y^2

	t2 = mul_prime(t1, t3, p);				//t2 = 8*Y^4
	t4 = sub_mod(t4, pd->x, p);				//t4 = 4*X*Y^2-X3
	t5 = mul_prime(t0, t4, p);				//t5 = 3*X^2*(4*X*Y^2-X3)
	pd->y = sub_mod(t5, t2, p);				//Y3 = 3*X^2*(4*X*Y^2-X3)-8*Y^4


	t2 = mul_prime(t0, T.x, p);				//t2 = 3*X^3
	tl->l0 = sub_mod(t2, t1, p);			//l0 = 3*X^3-2*Y^2

	t3 = mul_prime(T.z, T.z, p);			//t3 = Z^2
	t4 = mul_prime(t0, t3, p);				//t4 = 3*X^2*Z^2
	t4 = neg_mod(t4, p);					//t4 = -3*X^2*Z^2
	tl->l1.bn0 = mul_prime(t4, Q.x.bn0, p);
	tl->l1.bn1 = mul_prime(t4, Q.x.bn1, p);

	t5 = mul_prime(pd->z, t3, p);			//Z3*Z^2
	tl->l2.bn0 = mul_prime(t5, Q.y.bn0, p);
	tl->l2.bn1 = mul_prime(t5, Q.y.bn1, p);
}

/*
 * pa = T+P
 * s1 = secant_line(Q, T, P)
 */
void point_addition_and_secantline(point_jc *pa,
								line_result *sl,
								point_ac_ext2 Q,
								point_jc T,
								point_ac P,
								un_bn p)
{
	un_bn t0,t1,t2,t3,t4,t5;

	t0 = mul_prime(T.z, T.z, p);				//t0 = z1^2
	t1 = mul_prime(P.x, t0, p);					//t1 = X2 * Z1^2
	t2 = sub_mod(t1, T.x, p);					//t2 = X2 * Z1^2 - X1
	t3 = add_mod(t1, T.x, p);					//t3 = X2 * Z1^2 + X1
	pa->z = mul_prime(T.z, t2, p);				//Z3 = Z1*(X2 * Z1^2 - X1)

	t4 = mul_prime(t2, t2, p);					//t4 = (X2 * Z1^2 - X1)^2
	t1 = mul_prime(t0, T.z, p);					//t1 = z1^3
	t1 = mul_prime(P.y, t1, p);					//t1 = Y2*z1^3
	t1 = sub_mod(t1, T.y, p);					//t1 = Y2*z1^3-Y1
	t5 = mul_prime(t1, t1, p);					//t5 = (Y2*z1^3-Y1)^2
	t3 = mul_prime(t3, t4, p);					//t3 = (X2 * Z1^2 - X1)^2 * (X2 * Z1^2 + X1)
	pa->x = sub_mod(t5, t3, p);					//X3 = (Y2*z1^3-Y1)^2 - (X2 * Z1^2 - X1)^2 * (X2 * Z1^2 + X1)

	t3 = mul_prime(t2, t4, p);					//t3 = (X2 * Z1^2 - X1)^3
	t4 = mul_prime(T.x, t4, p);					//t4 = x1 * (X2 * Z1^2 - X1)^2
	t4 = sub_mod(t4, pa->x, p);					//t4 = x1 * (X2 * Z1^2 - X1)^2 - X3
	t5 = mul_prime(t1, t4, p);					//t5 = (Y2*z1^3-Y1)*(x1 * (X2 * Z1^2 - X1)^2 - X3)
	t3 = mul_prime(T.y, t3, p);					//t3 = Y1*(X2 * Z1^2 - X1)^3
	pa->y = sub_mod(t5, t3, p);					//Y3 = (Y2*z1^3-Y1)*(x1 * (X2 * Z1^2 - X1)^2 - X3) - Y1*(X2 * Z1^2 - X1)^3

	t0 = mul_prime(P.x, t1, p);					//t0 = X2*(Y2*z1^3-Y1)
	t2 = mul_prime(P.y, pa->z, p);				//t2 = Y2*Z3
	sl->l0 = sub_mod(t0, t2, p);				//l0 = X2*(Y2*z1^3-Y1) - Y2*Z3

	t1 = neg_mod(t1, p);						//t1 = -(Y2*z1^3-Y1)
	sl->l1.bn0 = mul_prime(t1, Q.x.bn0, p);
	sl->l1.bn1 = mul_prime(t1, Q.x.bn1, p);

	sl->l2.bn0 = mul_prime(pa->z, Q.y.bn0, p);
	sl->l2.bn1 = mul_prime(pa->z, Q.y.bn1, p);
}

















