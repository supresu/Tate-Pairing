#include "miller.h"
#include <stdio.h>


un_bn scanf_prime(FILE *fp){
	un_bn res;
	int i;
	unint tmp;
	//FILE *fp;

	memset(&res, 0, sizeof(un_bn));

	//fp = fopen("./data/data.txt", "r");

	for (i = SIZE - 1; i >= 0; i--) {
		fscanf(fp, "%4X", &tmp);
		res.bn[i] = tmp;
	}

	//fclose(fp);

	return res;
}

int main(){
	//p=0x2523648240000001BA344D80000000086121000000000013A700000000000013
	un_bn p = {{0x13,0,0,0xA700,0x13,0,0,0x6121,0x8,0,0x4D80,0xBA34,0x1,0x4000,0x6482,0x2523}, 0};

	point_jc P,tP;
	point_ac_ext2 Q, tQ;
	un_bn_ext12_t res,ts;

	FILE *fp;

	fp = fopen("./data/data.txt", "r");

	if(fp==NULL)
		printf("Cann't find the data file in the direction \"./data\" !\n");
	else{
		P.x = scanf_prime(fp);
		P.y = scanf_prime(fp);
		memset(&P.z, 0, sizeof(un_bn));
		P.z.bn[0] = 1;
		printf("P = (");printf_res(P.x);printf(" :\n");
		printf("     ");printf_res(P.y);printf(" :\n");
		printf("     ");printf("1 )\n");

		Q.x.bn1 = scanf_prime(fp);
		Q.x.bn0 = scanf_prime(fp);
		Q.y.bn1 = scanf_prime(fp);
		Q.y.bn0 = scanf_prime(fp);
		printf("Q = (");printf_res(Q.x.bn0);printf(" + ");printf_res(Q.x.bn1);printf("*u :\n");
		printf("     ");printf_res(Q.y.bn0);printf(" + ");printf_res(Q.y.bn1);printf("*u )\n");

		printf("\nComputing the tate pairing...\n\n");

		res = miller_tate(P, Q, p);

		printf("T(P, Q) = \n(\n");print_ext6(res.W1);printf(")*w\n");
		print_ext6(res.W0);

		printf("\n");
		printf_mul_counts();
		printf_inv_counts();



		printf("\nTesting the property...\n\n");

		//T(P, Q)^6
		ts  = square_ext12(res, p);
		res = square_ext12(ts,  p);
		res = mul_ext12(res, ts, p);
		printf("T(P, Q)^6 = \n(\n");print_ext6(res.W1);printf(")*w\n");
		print_ext6(res.W0);
		printf("\n");

		//T(2*P, 3*Q)
		tP.x = scanf_prime(fp);
		tP.y = scanf_prime(fp);
		memset(&tP.z, 0, sizeof(un_bn));
		tP.z.bn[0] = 1;
		printf("2*P = (");printf_res(tP.x);printf(" :\n");
		printf("       ");printf_res(tP.y);printf(" :\n");
		printf("       ");printf("1 )\n");

		tQ.x.bn1 = scanf_prime(fp);
		tQ.x.bn0 = scanf_prime(fp);
		tQ.y.bn1 = scanf_prime(fp);
		tQ.y.bn0 = scanf_prime(fp);
		printf("3*Q = (");printf_res(tQ.x.bn0);printf(" + ");printf_res(tQ.x.bn1);printf("*u :\n");
		printf("       ");printf_res(tQ.y.bn0);printf(" + ");printf_res(tQ.y.bn1);printf("*u )\n");

		printf("\n");

		res = miller_tate(tP, tQ, p);

		printf("T(2*P, 3*Q) = \n(\n");print_ext6(res.W1);printf(")*w\n");
		print_ext6(res.W0);

		fclose(fp);
	}

	//system("pause");

	return 0;
}
