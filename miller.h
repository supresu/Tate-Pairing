/*
 * miller.h
 *
 *  Created on: 2013-4-15
 *      Author: Suliangjian
 */

#ifndef MILLER_H_
#define MILLER_H_

#include "final_exp.h"

unsigned long long millercycles;

point_ac 		jc_to_ac(point_jc, un_bn);
un_bn_ext12_t 	miller_function(point_jc, point_ac_ext2, un_bn);
un_bn_ext12_t 	miller_tate(point_jc, point_ac_ext2, un_bn);
void 			print_cycles();

#endif /* MILLER_H_ */
