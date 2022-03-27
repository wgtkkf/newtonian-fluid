#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#define size0 3
#define size1 2
#define nmax 191 /* 190 + 1 */
#define emax 325 /* 324 + 1 */

#define bcmax1 50
#define bcmax2 50
#define bcmax3 50
#define bcmax4 50
#define bcmax5 50
#define bcmax6 50
#define bcmax  50

void start();
void end();
void check();
void output();

/* main program */
int main(void)
{
	/* counters */
	int i, j, k, n;
	int counter=0;
	int stepnum=0;
	double time, tmax=0.5;
	
	/* matrixs */
	double area, area2, a2, area3, area12;
	int n1, n2, n3;
	double nop[emax][size0]={};
	double cord[nmax][size1]={};
	double b[size0]={}, c[size0]={}, dd[size0][size0]={};
	double eb[emax][size0]={}, ec[emax][size0]={};
	double x1, x2, x3, y1, y2, y3;
	double emm[emax][size0]={}, ehx[emax][size0][size0]={}, ehy[emax][size0][size0]={}, ess[emax][size0][size0]={};
	double essxbb[emax][size0][size0]={}, essxbc[emax][size0][size0]={}, essycc[emax][size0][size0]={}, essycb[emax][size0][size0]={};
	double lmm[nmax]={}, ilmm[nmax]={};
	double ekx[emax][size0][size0][size0]={}, eky[emax][size0][size0][size0]={};
	int elem;
	
	/* step1, 2, 3 */
	double nbc1, nnbc1[bcmax1];
	double nbc2, nnbc2[bcmax2];
	double nbc3, nnbc3[bcmax3];
	double nbc4, nnbc4[bcmax4];
	double nbc5, nnbc5[bcmax5];
	double nbc6, nnbc6[bcmax6];
	double vbc1[bcmax1], vbc2[bcmax2], vbc3[bcmax3], vbc4[bcmax4], vbc5[bcmax5], vbc6[bcmax6]; 
		
	/* the parameters for newton flow */
	FILE *fp0;
	double rho;
	double length;
	double velocity;
	double mu;
	double re;
	double rei;
	
	/* element informations */
	FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8;
	char s[20];
	int np;
	int ne;
	double dt;
	
	/* velocitys, pressures */
	double vx0[nmax]={}, vy0[nmax]={}, pr0[nmax]={};
	double vx1[nmax]={}, vy1[nmax]={}, pr1[nmax]={};
	double xx[nmax]={};
	double vvx[nmax]={}, vvy[nmax]={};
	double dudx, dvdy, duyvx;
	double ii, kk;
	/* stress tensor */
	double txx1[nmax]={}, txy1[nmax]={}, tyy1[nmax]={};
	
	/* cg method */
	double r[nmax];
	double p[nmax];
	double vectorb[nmax];
	double b2;
	double eps, eps1, eps2;
	double delta;
	int kend, icount;
	double ap[nmax];
	double rur0, rur1;
	double res2;
	double pap;
	double alpha;
	double beta;
	double v0, v1;
	double tmp;
	double vmean;
	double res;
	
	double nbcn;
	double nnbcn[bcmax][size1]={};/* check the size of table */
	double vbcn[bcmax];/* check the size of table */
	
	int N, M, PG;
	double LL;
	
	double std = 10E-8;/* convergence */
	
	/* output file */
	char filename[256];
	FILE *fp;
	div_t d;
	
	/* cpu time */
	clock_t stime, ftime;
	
	/* start of the main program */
	start();
	stime=clock();

	/* 1. read the parameters for newton flow */
	fp0=fopen("parameters.txt","r");
	if(fp0==NULL)
	{
		printf("failure\n");
        return -1;
	}
		
	fscanf(fp0,"%lf",&rho);
	fscanf(fp0,"%lf",&length);
	fscanf(fp0,"%lf",&velocity);
	fscanf(fp0,"%lf",&mu);
	fclose(fp0);
	re = (rho*length*velocity)/mu;
		  
	printf("parameters for newton flow\n");
	printf("density:%lf\n", rho);
	printf("characteristic length:%lf\n", length);
	printf("characteristic velocity:%lf\n", velocity);
	printf("viscosity:%lf\n", mu);

	printf("reynolds number:%lf\n", re);

	printf("\n");

	/* 2. read the element information */
	/* 2.1 read the mesh information */
	fp1=fopen("rectangle1.msh","r");
	if(fp1==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	fp2=fopen("rectangle2.msh","r");
	if(fp2==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	fp3=fopen("rectangle3.msh","r");
	if(fp3==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	
	/* *fp1: rectangle1.msh */
	while(fgets(s, 20, fp1)!=NULL){
		if(counter==1)
		{
			np = atoi(s);
		}else if(counter==2){
			ne = atoi(s);
		}else if(counter==3){
			dt = atof(s);
		}
		counter++;
	};
	counter=0;
	
	printf("number of the nodes:%d\n", np);
	printf("number of the elements:%d\n", ne);
	printf("dt:%lf\n", dt);
	
	/* show the cord information */
	/*printf("No. x y\n");*/
	/* read the coordinates of the mesh */
	for(i=1;i<=np;i++)
	{
		fscanf(fp2, "%d %lf %lf", &n, &cord[i][0], &cord[i][1]);
		/*printf("%d %lf %lf\n", n, cord[i][0], cord[i][1]);*/
	}
		
	/* show the node information */
	/*printf("nodes information\n");*/
	for(i=1;i<=ne;i++)
	{
		fscanf(fp3, "%d %lf %lf %lf", &n, &nop[i][0], &nop[i][1], &nop[i][2]);
		/*printf("%d %lf %lf %lf\n", n, nop[i][0], nop[i][1], nop[i][2]);*/
	}
	
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);

	/* 3. read the initial conditions */
	/* initialize */
	for(i=1;i<=np;i++)
	{
		vx0[i] = 0.0;
		vy0[i] = 0.0;
		pr0[i] = 0.0;
		vx1[i] = 0.0;
		vy1[i] = 0.0;
		pr1[i] = 0.0;
	}
	
	/* input initial data */
		
	/* 4. read the boundary conditions */
	fp4=fopen("rectangle-num.bc","r");
	if(fp4==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	fp5=fopen("rectangle-u.bc","r");
	if(fp4==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	fp6=fopen("rectangle-v.bc","r");
	if(fp4==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	fp7=fopen("rectangle-p.bc","r");
	if(fp6==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	fp8=fopen("rectangle-n.bc","r");
	if(fp6==NULL){
		printf("cannot open the file\n");
		return -1;
	}
	
	/* read the number of each condition U, V, P, N.B.*/
	fscanf(fp4, "%lf", &nbc1);
	fscanf(fp4, "%lf", &nbc2);
	fscanf(fp4, "%lf", &nbc3);
	fscanf(fp4, "%lf", &nbcn);
	fscanf(fp4, "%lf", &vmean);

	/* boundary condition for U */
	for(i=1;i<=(int)nbc1;i++)
	{
		fscanf(fp5, "%lf %lf", &nnbc1[i], &vbc1[i]);
	}
	/* boundary condition for V */
	for(i=1;i<=(int)nbc2;i++)
	{
		fscanf(fp6, "%lf %lf", &nnbc2[i], &vbc2[i]);
	}
	/* boundary condition for P */
	for(i=1;i<=(int)nbc3;i++)
	{
		fscanf(fp7, "%lf %lf", &nnbc3[i], &vbc3[i]);
	}
	/* boundary condition for natural B.C */
	for(i=1;i<=(int)nbcn;i++)
	{
		fscanf(fp8, "%d %d %d", &N, &M, &PG);
		if(M==1)
		{
			n1 = nop[n][1];
			n2 = nop[n][2];
		}else if(M==2)
		{
			n1 = nop[n][0];
			n2 = nop[n][2];
		}else if(M==3)
		{
			n1 = nop[n][0];
			n2 = nop[n][1];
		}
		nnbcn[i][0]=n1;
		nnbcn[i][1]=n2;
		
		x1 = cord[n1][0];
		y1 = cord[n1][1];
		x2 = cord[n2][0];
		y2 = cord[n2][1];
		LL = sqrt(pow((x1-x2),2) + pow((y1-y2),2));
		vbcn[i] = PG*LL;
	}
	fclose(fp4);
	fclose(fp5);
	fclose(fp6);
	fclose(fp7);
	fclose(fp8);

	/* 5. calculation of the element matrix */
	for(elem=1;elem<=ne;elem++)
	{
		/* initialization */
		for(i=0;i<size0;i++)
		{
			emm[elem][i] = 0.0;
			for(j=0;j<size0;j++)
			{
				ehx[elem][i][j] = 0.0;
				ehy[elem][i][j] = 0.0;
				ess[elem][i][j] = 0.0;
				if(i==j)
				{
					dd[i][j] = 2.0;
				}else{
					dd[i][j] = 1.0;
				}
				for(k=0;k<size0;k++)
				{
					ekx[elem][i][j][k] = 0.0;
					eky[elem][i][j][k] = 0.0;
				}
			}
		}
		
		/* calculation area of element */
		/* n1 -- number of 1st node */
		/* n2 -- number of 1st node */
		/* n3 -- number of 1st node */
		/* area -- area of element */
		n1 = nop[elem][0];
		n2 = nop[elem][1];
		n3 = nop[elem][2];
		x1 = cord[n1][0];
		y1 = cord[n1][1];
		x2 = cord[n2][0];
		y2 = cord[n2][1];
		x3 = cord[n3][0];
		y3 = cord[n3][1];
		
		area2 = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
		area = 0.5*area2;		
		
		if(area<0)
		{
			printf("error: area is less than zero\n");
			exit(1);
		}
		
		a2 = 1/area2;
		area3 = area/3;
		area12 = area/12;
				
		/*
		printf("area2:%lf\n", area2);
		printf("a2:%lf\n", a2);
		printf("area3:%lf\n", area3);
		printf("area12:%lf\n", area12);		
		printf("x1:%lf x2:%lf x3:%lf\n", x1, x2, x3);
		printf("y1:%lf y2:%lf y3:%lf\n", y1, y2, y3);
		 */
		
		/* differential of shape function */
		b[0] = (y2-y3)*a2;
		b[1] = (y3-y1)*a2;
		b[2] = (y1-y2)*a2;
		c[0] = (x3-x2)*a2;
		c[1] = (x1-x3)*a2;
		c[2] = (x2-x1)*a2;
		
		/*		
		printf("b[0]:%lf b[1]:%lf b[2]:%lf\n", b[0], b[1], b[2]);
		printf("c[0]:%lf c[1]:%lf c[2]:%lf\n", c[0], c[1], c[2]);
		 */
	
		eb[elem][0] = b[0];
		eb[elem][1] = b[1];
		eb[elem][2] = b[2];
		ec[elem][0] = c[0];
		ec[elem][1] = c[1];
		ec[elem][2] = c[2];
		
		/* calculation of matrix */
		/* Hx, Hy */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				ehx[elem][i][j] = area3*b[j];
				ehy[elem][i][j] = area3*c[j];
			}
		}
		
		/* Kxx, Kyy */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				for(k=0;k<size0;k++)
				{
					ekx[elem][i][j][k] = area12*dd[i][j]*b[k];
					eky[elem][i][j][k] = area12*dd[i][j]*c[k];
				}
			}
		}
		
		/* S */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				ess[elem][i][j] = area*(b[i]*b[j]+c[i]*c[j]);
			}
		}
		
		/* Sxbb */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				essxbb[elem][i][j] = 2*area*(b[i]*b[j]);
			}
		}

		/* Sxbc */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				essxbc[elem][i][j] = area*(b[i]*c[j]);
			}
		}
				
		/* Sycc */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				essycc[elem][i][j] = 2*area*(c[i]*c[j]);
			}
		}
		
		/* Sycb */
		for(i=0;i<size0;i++)
		{
			for(j=0;j<size0;j++)
			{
				essycb[elem][i][j] = area*(c[i]*b[j]);
			}
		}
		
		/*  -                    */
		/* [M]lumped mass matrix */
		for(i=0;i<size0;i++)
		{
			emm[elem][i] = area3;
			n = nop[elem][i];
			lmm[n] = lmm[n] + area3;
			/*
			 printf("emm[%d][%d]:%lf\n", elem, i, emm[elem][i]);
			 printf("n:%d\n", n);
			 printf("lmm[%d]:%lf\n", n, lmm[n]);
			 */
		}
	}
	
	/* inverse matrix of the element matrix */
	for(i=1;i<=np;i++)
	{
		ilmm[i] = 1.0/lmm[i];
	}	
		
	/* start calculation */
	start();
	
	/* 4. start time loop */
	while(time<=tmax){
		
		if(tmax<time)
		{
			break;
		}
		
		/* 5. step 1, 2, 3 */
		/* 5.1 step1 */
		/* calculate velocity */
		
		/* calculate r.h.s vector (u) */
		/* initialization */
		for(i=1;i<=np;i++)
		{
			vvx[i] = 0.0;
		}
		rei = 1.0/re;
		
		for(n=1;n<=ne;n++)
		{
			/*
			n1 = nop[n][0];
			n2 = nop[n][1];
			n3 = nop[n][2];
			 */
			
			for(i=0;i<size0;i++)
			{
				n1 = nop[n][i];
				for(j=0;j<size0;j++)
				{
					n2 = nop[n][j];
					vvx[n1] = vvx[n1] - rei*essxbb[n][i][j]*vx1[n2]- rei*essxbc[n][i][j]*vy1[n2];
					for(k=0;k<size0;k++)
					{
						n3 = nop[n][k];
						vvx[n1] = vvx[n1] - ekx[n][i][j][k]*vx1[n2]*vx1[n3] - eky[n][i][j][k]*vy1[n2]*vx1[n3];
					}
				}
			}
		}
		
		/* boundary condition */
		for(i=1;i<=nbc1;i++)
		{
			vvx[(int)nnbc1[i]] = 0.0;
		}
		
		/* calculate r.h.s vector (v) */
		/* initialization */
		for(i=1;i<=np;i++)
		{
			vvy[i] = 0.0;
		}
		rei = 1.0/re;
	
		for(n=1;n<=ne;n++)
		{
			/*
			n1 = nop[n][0];
			n2 = nop[n][1];
			n3 = nop[n][2];
			*/
			
			for(i=0;i<size0;i++)
			{
				n1 = nop[n][i];
				for(j=0;j<size0;j++)
				{
					n2 = nop[n][j];
					vvy[n1] = vvy[n1]  - rei*essycb[n][i][j]*vx1[n2] - rei*essycc[n][i][j]*vy1[n2];
					for(k=0;k<size0;k++)
					{
						n3 = nop[n][k];
						vvy[n1] = vvy[n1] - ekx[n][i][j][k]*vx1[n2]*vx1[n3] - eky[n][i][j][k]*vy1[n2]*vy1[n3];
					}
				}
			}
		}
		
		/* boundary condition */
		for(i=1;i<=nbc1;i++)
		{
			vvy[(int)nnbc2[i]]=0.0;
		}
		
		/* boundary condition (u) */
		for(i=1;i<=nbc1;i++)
		{
			n = nnbc1[i];
			vx0[n] = vbc1[i];
		}
		
		/* explicit method */
		for(i=1;i<=np;i++)
		{
			vx1[i] = vx0[i] + dt*vvx[i]*ilmm[i];
		}
		
		/* boundary condition (v) */
		for(i=1;i<=nbc2;i++)
		{
			n = nnbc2[i];
			vy0[(int)n] = vbc2[i];
		}
		
		/* explicit method */
		for(i=1;i<=np;i++)
		{
			vy1[i] = vy0[i] + dt*vvy[i]*ilmm[i];
		}
		
		/* boundary condition (tij) */
		for(i=1;i<=nbc4;i++)
		{
			n = nnbc4[i];
			txx1[n] = vbc4[i];
		}
		
		for(i=1;i<=nbc5;i++)
		{
			n = nnbc5[i];
			txy1[n] = vbc5[i];
		}
		
		for(i=1;i<=nbc6;i++)
		{
			n = nnbc6[i];
			tyy1[n] = vbc6[i];
		}
		
		/* 5.2 step2 */
		/* calculate the poison equation by cg method */
		/* cg method */
		/* initialization */
		
		eps1 = 0.0001;
		eps2 =eps1*eps1;
		kend=1000;
		delta = -1.0/dt;
		for(i=1;i<=np;i++)
		{
			r[i] = 0.0;
			vectorb[i] = 0.0;
			xx[i] = 0.0;
		}
		
		/* r0=b-Ax0 */
		for(n=1;n<=ne;n++)
		{
			for(i=0;i<size0;i++)
			{
				n1=nop[n][i];
				for(j=0;j<size0;j++)
				{
					n2 = nop[n][j];
					vectorb[n1]=vectorb[n1]+delta*(ehx[n][i][j]*vx1[n2]+ehy[n][i][j]*vy1[n2]);
				}
			}
		}
		
		/* natural boundary condition */
		for(i=1;i<=nbcn;i++)
		{
			n1 = nnbcn[i][0];
			n2 = nnbcn[i][1];
			vectorb[n1] = vectorb[n1] + 0.5*vbcn[i];
			vectorb[n2] = vectorb[n2] + 0.5*vbcn[i];
		}
		
		/* boundary condition(dp=0) */
		for(i=1;i<=nbc3;i++)
		{
			vectorb[(int)nnbc3[i]] = 0.0;
		}
		
		/* p0 = r0 */
		for(i=1;i<=np;i++)
		{
			r[i] = vectorb[i];
			p[i] = vectorb[i];
		}
	
		for(i=1;i<=np;i++)
		{
			b2 = b2 + vectorb[i]*vectorb[i];
		}
		
		if(b2<eps2)
		{
			eps = 0.0;
			icount = 1.0;
			kend = icount;
		}
		
		/* interative loop */
		icount=0;
		for(k=0;k<kend;k++)
		{
			icount=icount+1;
			for(i=1;i<=np;i++)
			{
				ap[i] = 0.0;
			}
			
			/* AP: ApK, P: pk */
			/* PAP: (pk, Apk) */
			for(n=1;n<=ne;n++)
			{
				for(i=0;i<size0;i++)
				{
					n1 = nop[n][i];
					for(j=0;j<size0;j++)
					{
						n2 = nop[n][j];
						ap[n1] = ap[n1] + ess[n][i][j]*p[n2];
					}
				}
			}
			
			/* Boundary Condition */
			for(i=1;i<=nbc3;i++)
			{
				ap[(int)nnbc3[i]] = 0.0;
			}	
			
			rur0 = 0.0;
			for(i=1;i<=np;i++)
			{
				rur0 = rur0 + p[i]*r[i];
			}
			
			for(i=1;i<=np;i++)
			{
				pap = pap + p[i]*ap[i];
			}
			
			if(pap==0)
			{
				alpha = 0.0;
			}else{
				alpha = rur0/pap;
			}
			
			for(i=1;i<=np;i++)
			{
				xx[i] = xx[i] + alpha*p[i];
				r[i] = r[i] - alpha*ap[i];
			}
			
			rur1=0.0;
			res2=0.0;
			for(i=1;i<=np;i++)
			{
				rur1 = rur1 + r[i]*ap[i];
				res2 = res2 + r[i]*r[i];
			}
			
			eps = res2/b2;
			if(eps<eps2)
			{
				kend =  icount;
			}
			
			/* not converged */
			beta = -rur1/pap;
			/* p(k+1) = r(k+1) + beta*p(k) */
			for(i=1;i<=np;i++)
			{
				p[i] = r[i] + beta*p[i];
			}
		}
		
		/* update of pressure */
		for(i=1;i<=np;i++)
		{
			pr1[i] = xx[i];
		}
		
		/* 5.3 step3 */
		/* initialization */
		for(i=1;i<=np;i++)
		{
			vvx[i] = 0.0;
			vvy[i] = 0.0;
		}
		
		/* calculate vector */
		for(n=1;n<=ne;n++)
		{
			for(i=0;i<size0;i++)
			{
				n1 = nop[n][i];
				for(j=0;j<size0;j++)
				{
					n2 = nop[n][j];
					vvx[n1] = vvx[n1] + ehx[n][i][j]*pr1[n2];
					vvy[n1] = vvy[n1] + ehy[n][i][j]*pr1[n2];
				}
			}
		}
		
		for(i=1;i<=np;i++)
		{
			vx1[i] = vx1[i] - dt*vvx[i]*ilmm[i];
			vy1[i] = vy1[i] - dt*vvy[i]*ilmm[i];
		}
		
		/* boundary condition */
		for(i=1;i<=nbc1;i++)
		{
			vx1[(int)nnbc1[i]] = vbc1[i];
		}
	
		for(i=1;i<=nbc2;i++)
		{
			vy1[(int)nnbc2[i]] = vbc2[i];
		}
		
		/* boundary condition */
		/*
		for(i=1;i<=nbc1;i++)
		{
			vvx[(int)nnbc1[i]] = 0.0;
		}*/
		
		/* 6. check of convergence */
		for(i=1;i<=np;i++)
		{
			v0 = sqrt(vx0[i]*vx0[i] + vy0[i]*vy0[i]);
			v1 = sqrt(vx1[i]*vx1[i] + vy1[i]*vy1[i]);
			tmp = fabs(v1-v0)/vmean;
			res = fmax(res, tmp);
		}
		
		if(time<tmax && res > std)
		{
			/* 7. Update Calculated Data */
			for(i=1;i<=np;i++)
			{
				vx0[i] = vx1[i];
				vy0[i] = vy1[i];
				pr0[i] = pr1[i];
			}
			
			/* format of AVese */		
			d = div(stepnum, 500);
			
			/* write to files */
			if(stepnum==0)
			{
				/* define the file name */
				sprintf(filename,"./files/step%d.dat",stepnum);
				/* creation of files */
				if((fp=fopen(filename, "w"))==NULL){
					fprintf(stderr, "%s\n", filename);
				}
				for(i=1;i<=np;i++)
				{
					fprintf(fp, "%lf %lf 0.0 0.0\n", cord[i][0], cord[i][1]);
				}
			}else if(stepnum>0 && d.rem==0)
			{
				/* define the file name */
				sprintf(filename,"./files/step%d.dat",stepnum);
				/* creation of files */
				if((fp=fopen(filename, "w"))==NULL){
					fprintf(stderr, "%s\n", filename);
				}
				for(i=1;i<=np;i++)
				{
					fprintf(fp, "%lf %lf %lf %lf\n", cord[i][0], cord[i][1], vx1[i], vy1[i]);
				}
			}
			
			
			fclose(fp);
			
			/* show to display */		
			for(i=1;i<=np;i++)
			{
				printf("node:%d vx1:%lf vy1:%lf 0.0\n",i, vx1[i], vy1[i]);
			}
			
			/* prpceed time */				
			time = time + dt;
			
			/* proceed the step */
			stepnum++;
		}
		
		/* 8. output of the result */
		/*
		sprintf(filename,"./files/step%d.vtk",stepnum);
		 */
		
		/* 8.1 creation of files */
		/*
		if((fp=fopen(filename, "w"))==NULL){
			fprintf(stderr, "%s\n", filename);
		}
		 */
		
		/* 8.2 write words to files */
		/* format of PARAVIEW */
		/*
		fprintf(fp, "# vtk DataFile Version 2.0\n");
		fprintf(fp, "vector\n");
		fprintf(fp, "ASCII\n");
		fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
		fprintf(fp, "POINTS %d float\n", np);
		 */
		/* cord information */
		/*
		for(i=1;i<=np;i++)
		{
			fprintf(fp, "%lf %lf 0.0\n", cord[i][0], cord[i][1]);
		}
		 */

		/* nord information */
		/*
		fprintf(fp, "CELLS %d 1296\n", ne);
		for(i=1;i<=ne;i++)
		{
			fprintf(fp, "3 %d %d %d\n", (int)nop[i][0], (int)nop[i][1], (int)nop[i][2]);
		}
		 */

		/*
		fprintf(fp, "CELL_TYPES %d\n",ne);
		for(i=1;i<=ne;i++)
		{
			fprintf(fp, "5\n");
		}
		
		fprintf(fp, "POINT_DATA %d\n", np);
		fprintf(fp, "VECTORS point_vectors float\n");
		*/
		
		/* velocity data */
		/*
		for(i=1;i<=np;i++)
		{
			fprintf(fp, "%lf %lf 0.0\n", vx1[i], vy1[i]);
		}
		 */
		 
		/*
		fprintf(fp, "CELL_DATA 1\n");
		fprintf(fp, "VECTORS cell_vectors float\n");
		 */
	}

	/* 9. end of the maim program */
	ftime = clock();
	printf("%lf second\n", (double)(ftime-stime)/(double)CLOCKS_PER_SEC);
	end();
	return 0;
}

/* functions */
void start()
{
	printf("start of the calculation\n\n");	
}

void end()
{
	printf("end of the calculation\n");	
}

void check()
{
	printf("proceeding..\n");	
}

void output(FILE *fp)
{
}