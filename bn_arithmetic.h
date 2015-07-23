/*
 * bn_arithmetic.h
 *
 *  Created on: 2013-3-25
 *      Author: Suliangjian
 */

#ifndef BN_ARITHMETIC_H_
#define BN_ARITHMETIC_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define SIZE 16
#define DIGIT_SIZE 16     //The big number is up to 16*16 bits
#define QUADRATIC_NR -2   //quadratic non-residue
#define LESS   -1
#define EQUAL  0
#define BIGGER 1
#define YES 0
#define NO  1


#define MAX(x,y) ((x)>(y)?(x):(y))
#define unint_swap(x, y)			\
	do {							\
		unint tmp = (x);			\
		(x) = (y);					\
		(y) = tmp;					\
	} while (0)

typedef unsigned int unint;
typedef unsigned long long undlong;
typedef int compareres;

typedef struct{
	unint bn[4];
	unint c;
}un_bn_64;

typedef struct{
	unint bn[8];
	unint c;
}un_bn_128;

typedef struct{
	unint bn[SIZE];
	unint c;
}un_bn;				//The element in prime field

typedef struct{
	un_bn bn0;
	un_bn bn1;
}un_bn_ext2;		//The element in the field with extension degree 2

typedef struct{
	un_bn_ext2 V0;
	un_bn_ext2 V1;
	un_bn_ext2 V2;
}un_bn_ext6;		//The element in the field with extension degree 6

typedef struct{
	un_bn_ext6 W0;
	un_bn_ext6 W1;
}un_bn_ext12_t;		//The element in the field with extension degree 12


typedef struct{
	un_bn bn0;
	un_bn bn1;
	un_bn bn2;
	un_bn bn3;
	un_bn bn4;
	un_bn bn5;
	un_bn bn6;
	un_bn bn7;
	un_bn bn8;
	un_bn bn9;
	un_bn bn10;
	un_bn bn11;
}un_bn_ext12;		//The element in the field with extension degree 12


typedef struct{
	un_bn l0;
	un_bn_ext2 l1;
	un_bn_ext2 l2;
}line_result;

unsigned long long mul_counts;
unsigned long long inv_counts;
//unsigned long long add_counts;

void		un_bn_cpy(un_bn *, un_bn);
compareres	un_bn_cmp(un_bn, un_bn);
compareres	un_bn_equal1(un_bn);
compareres	un_bn_even(un_bn);
un_bn		un_bn_negation(un_bn);
void		un_bn_binary(unint [], unint);
un_bn		un_bn_sll(un_bn, int);
un_bn 		un_bn_srl(un_bn, int);
un_bn_64	add_64(un_bn_64, un_bn_64, unint);
un_bn_128	add_128(un_bn_128, un_bn_128, unint);
un_bn		add_256(un_bn, un_bn, unint);

un_bn 		imml(un_bn , un_bn , un_bn );
un_bn		add_mod(un_bn, un_bn, un_bn);
un_bn		sub_mod(un_bn, un_bn, un_bn);
un_bn		neg_mod(un_bn, un_bn);
un_bn		mul_prime(un_bn, un_bn, un_bn);
un_bn		inverse_mod(un_bn, un_bn);
void 		printf_mul_counts();
void 		printf_inv_counts();

void 		ext2_cop(un_bn_ext2 *, un_bn_ext2);
un_bn_ext2	neg_ext2(un_bn_ext2, un_bn);
un_bn_ext2	add_ext2(un_bn_ext2, un_bn_ext2, un_bn);
un_bn_ext2 	ext2_add_prime(un_bn_ext2, un_bn, un_bn);		//1
un_bn_ext2	sub_ext2(un_bn_ext2, un_bn_ext2, un_bn);
un_bn_ext2	mul_ext2(un_bn_ext2, un_bn_ext2, un_bn);
un_bn_ext2 	ext2_mul_prime(un_bn_ext2, un_bn, un_bn);
un_bn_ext2	square_ext2(un_bn_ext2, un_bn);
un_bn_ext2	inverse_ext2(un_bn_ext2, un_bn);
un_bn_ext2	ext2_power_ext(un_bn_ext2, un_bn, un_bn);
un_bn_ext2 	compute_xipow(un_bn);

void 		ext6_cop(un_bn_ext6 *, un_bn_ext6);
un_bn_ext6	neg_ext6(un_bn_ext6, un_bn);
un_bn_ext6	add_ext6(un_bn_ext6, un_bn_ext6, un_bn);
un_bn_ext6	sub_ext6(un_bn_ext6, un_bn_ext6, un_bn);
un_bn_ext2	nr_ext6(un_bn_ext2, un_bn);
un_bn_ext6	mul_ext6(un_bn_ext6, un_bn_ext6, un_bn);
un_bn_ext6 	ext6_mul_sp1(un_bn_ext6, un_bn, un_bn_ext2,un_bn);	//1
un_bn_ext6 	ext6_mul_sp2(un_bn_ext6, un_bn_ext2,un_bn);			//1
un_bn_ext6	square_ext6(un_bn_ext6, un_bn);
un_bn_ext6	inverse_ext6(un_bn_ext6, un_bn);

void 			ext12_cop(un_bn_ext12_t *, un_bn_ext12_t);
un_bn_ext12_t	neg_ext12(un_bn_ext12_t, un_bn);
un_bn_ext6		nr_ext12(un_bn_ext6, un_bn);
un_bn_ext12_t	square_ext12(un_bn_ext12_t, un_bn);
un_bn_ext12_t	mul_ext12(un_bn_ext12_t, un_bn_ext12_t, un_bn);
un_bn_ext12_t 	inverse_ext12(un_bn_ext12_t, un_bn);
un_bn_ext12_t 	ext12_mul_lineresult(un_bn_ext12_t, line_result, un_bn);	//1
un_bn_ext12_t	ext12_to_ext6(un_bn_ext12);
un_bn_ext12		ext6_to_ext12(un_bn_ext12_t);

void			printf_res(un_bn);
void 			print_ext6(un_bn_ext6);


#endif /* BN_ARITHMETIC_H_ */
