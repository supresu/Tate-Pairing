/*
 * final_exp.h
 *
 *  Created on: 2013-4-3
 *      Author: Suliangjian
 */

#ifndef FINAL_EXP_H_
#define FINAL_EXP_H_

#include "bn_curve_arithmetic.h"


un_bn_ext12_t 	ext12_frobenius_p(un_bn_ext12_t, un_bn);
un_bn_ext12_t	ext12_frobenius_p2(un_bn_ext12_t, un_bn);
un_bn_ext12_t 	ext12_power_ext(un_bn_ext12_t, un_bn, un_bn);
un_bn_ext12_t 	ext12_power_e(un_bn_ext12_t, un_bn);
un_bn_ext12_t	first_step(un_bn_ext12_t, un_bn);
un_bn_ext12_t	second_step(un_bn_ext12_t, un_bn);
un_bn_ext12_t	hard_step(un_bn_ext12_t, un_bn);

#endif /* FINAL_EXP_H_ */
