/*
 * bn_curve_arithmetic.h
 *
 *  Created on: 2013-3-28
 *      Author: Suliangjian
 */

#ifndef BN_CURVE_ARITHMETIC_H_
#define BN_CURVE_ARITHMETIC_H_


#include "bn_arithmetic.h"

typedef struct{
	un_bn x;
	un_bn y;
	un_bn z;
}point_jc;	//Jacobian coordinate

typedef struct{
	un_bn x;
	un_bn y;
}point_ac;	//affine coordinate

//The affine coordinate of the element in the field with extension degree 2
typedef struct{
	un_bn_ext2 x;
	un_bn_ext2 y;
}point_ac_ext2;



point_jc		point_double(point_jc, un_bn);
point_jc		point_addition(point_jc, point_ac, un_bn);
line_result		tangent_line(point_ac_ext2, point_jc, un_bn);
line_result		secant_line(point_ac_ext2, point_jc, point_ac, un_bn);

void point_double_and_tangentline(point_jc *,line_result *,point_ac_ext2,point_jc,un_bn);
void point_addition_and_secantline(point_jc *,line_result *,point_ac_ext2,point_jc,point_ac,un_bn);

#endif /* BN_CURVE_ARITHMETIC_H_ */
