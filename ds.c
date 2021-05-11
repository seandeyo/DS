#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

	
int nodes,D,*deg,(**inx)[2],totdeg;
double epsilon,eta,beta,**x,**xA,**xB,**xR,xerr,*y,*yA,*yB,*yR,yerr,(*yrank)[2];
char errfile[50],statsfile[50],solfile[50];


static inline double sq(double diff)
  {
  return diff*diff;
  }

	
int compare_error( const void *pa, const void *pb )
	{
	const double (*a)[2] = pa;
	const double (*b)[2] = pb;
	if ( (*a)[1] < (*b)[1] ) return +1;
	if ( (*a)[1] > (*b)[1] ) return -1;
	return 0;
	}


void projA(double **xo,double *yo)
	{
	int n,i,imax[nodes];
	double xmax[nodes];
	for(n=0;n<nodes;++n)
		{
		for(i=0;i<deg[n];++i)
			{
			xA[inx[n][i][0]][inx[n][i][1]]=(xo[inx[n][i][0]][inx[n][i][1]]<.5?0.:1.);
			if(xo[inx[n][i][0]][inx[n][i][1]]>xmax[n] || i==0)
				{
				imax[n]=i;
				xmax[n]=xo[inx[n][i][0]][inx[n][i][1]];
				}
			}
		xA[inx[n][imax[n]][0]][inx[n][imax[n]][1]]=1.;
		yrank[n][0]=n;
		yrank[n][1]=2*yo[n]*eta+(xmax[n]<.5?1-2*xmax[n]:0.);	
		yA[n]=0.;
		}

	qsort(yrank, nodes, 2*sizeof(double), compare_error);

	for(i=0;i<D;++i)
		if(yrank[i][1]>eta*eta)
			{
			n=(int)(yrank[i][0]+.5);
			yA[n]=eta;
			if(xmax[n]<.5)
				xA[inx[n][imax[n]][0]][inx[n][imax[n]][1]]=0.;
			}
	}

	
void projB(double **xo,double *yo)
	{
	double c;
	int n,i;

	for(n=0;n<nodes;++n)
		{
		c=eta*yo[n];
		for(i=0;i<deg[n];++i)
			c+=xo[n][i];
		c/=eta*eta+deg[n];
		yB[n]=c*eta;
		for(i=0;i<deg[n];++i)
			xB[n][i]=c;
		}
	}
	
	
void ref(double **xo,double *yo)
	{
	int n,i;
	
	for(n=0;n<nodes;++n)
		{
		yR[n]=2.*yo[n]-y[n];
		for(i=0;i<deg[n];++i)
			xR[n][i]=2.*xo[n][i]-x[n][i];
		}
	}
	
	
void RRR()
	{
	int n,i;
	double diff;
	if(beta>0.)
		{
		projA(x,y);
		ref(xA,yA);
		projB(xR,yR);
		}
	else
		{
		projB(x,y);
		ref(xB,yB);
		projA(xR,yR);
		}
		
	xerr=0.;
	yerr=0.;
	
	for(n=0;n<nodes;++n)
		{
		diff=yB[n]-yA[n];	
		y[n]+=beta*diff;	
		yerr+=sq(diff/eta);
		y[n]/=eta;
		
		for(i=0;i<deg[n];++i)
			{
			diff=xB[n][i]-xA[n][i];
			x[n][i]+=beta*diff;	
			xerr+=sq(diff);
			}
		}
	
	xerr/=totdeg;
	yerr/=nodes;
	
	eta*=(1.+epsilon*(sqrt(2*yerr/(xerr+yerr))-1.));

	for(n=0;n<nodes;++n)
		y[n]*=eta;

	xerr=sqrt(xerr);
	yerr=sqrt(yerr);
	}
	
	
int getnet(char *netfile)
	{
	int n,i,in;
	FILE *fp;

	fp=fopen(netfile,"r");
	if(!fp)
		{
		printf("netfile not found\n");
		return 0;
		}
	
	fscanf(fp,"%d",&nodes);

	deg=malloc(nodes*sizeof(int));
	inx=malloc(nodes*sizeof(int*[2]));

	for(n=0;n<nodes;++n)
		deg[n]=0;
	totdeg=0;
	for(n=0;n<nodes;++n)
		{
		fscanf(fp,"%*d%d",&deg[n]);
		totdeg+=deg[n];
		
		inx[n]=malloc(deg[n]*sizeof(int[2]));
		for(i=0;i<deg[n];++i)
			fscanf(fp,"%d",&inx[n][i][0]);
		}
	
	fclose(fp);

	int index[nodes];
	for(n=0;n<nodes;++n)
		index[n]=0;
	for(n=0;n<nodes;++n)
		for(i=0;i<deg[n];++i)
			inx[n][i][1]=index[inx[n][i][0]]++;

	return 1;
	}
	
	
void makevars()
	{
	int n;
	
	x=malloc(nodes*sizeof(double*));
	xA=malloc(nodes*sizeof(double*));
	xB=malloc(nodes*sizeof(double*));
	xR=malloc(nodes*sizeof(double*));

	y=malloc(nodes*sizeof(double));
	yA=malloc(nodes*sizeof(double));
	yB=malloc(nodes*sizeof(double));
	yR=malloc(nodes*sizeof(double));
	yrank=malloc(nodes*sizeof(double[2]));
	
	for(n=0;n<nodes;++n)
		{
		x[n]=malloc(deg[n]*sizeof(double));
		xA[n]=malloc(deg[n]*sizeof(double));
		xB[n]=malloc(deg[n]*sizeof(double));
		xR[n]=malloc(deg[n]*sizeof(double));
		}
	}
	
	
double urand()
  {
  return ((double)rand())/RAND_MAX;
  }


void randstart()
	{
	int n,i;
	
	for(n=0;n<nodes;++n)
		{
		y[n]=eta*urand();
		for(i=0;i<deg[n];++i)
			x[n][i]=urand();
		}
	}

	
int solve(int maxiter,int iterstride,double stoperr)
	{
	int n,i,iter;
	FILE *fp;
	
	fp=fopen(errfile,"w");
	
	for(iter=1;iter<=maxiter;++iter)
		{
		RRR();
		
		if(iter%iterstride==0)
			fprintf(fp,"%.6f\n",sqrt(sq(xerr)+sq(yerr)));
			//{
			//for(n=0;n<nodes;++n)
			//	fprintf(fp,"%.6f,",yB[n]/eta);
			//fprintf(fp,"%.6f,%.6f,%.6f\n",xerr,yerr,eta);
			//}
		
		if(xerr+yerr<stoperr)
			{
			fprintf(fp,"%.6f\n",sqrt(sq(xerr)+sq(yerr)));
			//for(n=0;n<nodes;++n)
			//	fprintf(fp,"%.6f,",yB[n]/eta);
			//fprintf(fp,"%.6f,%.6f,%.6f\n",xerr,yerr,eta);
			fclose(fp);
			return iter;
			}
		}
			
	fclose(fp);
			
	return 0;
	}

	
void printsol(char *solfile)
	{
	FILE *fp;
	int n,i,wsol;
	
	fp=fopen(solfile,"a");
	
	for(n=0;n<nodes;++n)
		if(yA[n]>eta/2.)
			fprintf(fp,"%d ",n);
	fprintf(fp,"\n----------\n");
	
	fclose(fp);
	}
	
	
int main(int argc,char* argv[])
	{
	char *netfile,*id;
	int iter,maxiter,iterstride,trials,t,c,solcount;
	double stoperr,aveiter,elapsed,iterpersec;
	FILE *fp;
	clock_t start;
	
	if(argc==10)
		{
		netfile=argv[1];
		D=atoi(argv[2]);
		epsilon=atof(argv[3]);
		beta=atof(argv[4]);
		maxiter=atoi(argv[5]);
		iterstride=atoi(argv[6]);
		stoperr=atof(argv[7]);
		trials=atoi(argv[8]);
		id=argv[9];
		}
	else
		{
		printf("expected nine arguments: netfile, D, epsilon, beta, maxiter, iterstride, stoperr, trials, id\n");
		return 1;
		}
		
	sprintf(errfile,"%s.err",id);
	sprintf(statsfile,"%s.stats",id);
	sprintf(solfile,"%s.sol",id);
	
	if(!getnet(netfile))
		return 1;
	makevars();

	fp=fopen(solfile,"w");
	fclose(fp);
	
	fp=fopen(statsfile,"w");
	for(c=0;c<argc;++c)
    		fprintf(fp,"%s ",argv[c]);
	fprintf(fp,"\n\n");
	fclose(fp);
	
	srand(time(NULL));
	start=clock();
	
	solcount=0;
	aveiter=0.;
	for(t=0;t<trials;++t)
		{
		eta=1.;
		randstart();
		//sprintf(errfile,"%s%d.err",id,t);
		iter=solve(maxiter,iterstride,stoperr);
		if(iter)
			{
			++solcount;
			aveiter+=iter;
			}
		else
			aveiter+=maxiter;
		
		printsol(solfile);
		fp=fopen(statsfile,"a");
		fprintf(fp,"%3d%12d%12.2f\n",t,iter,eta);
		fclose(fp);
		}
	
	elapsed=((double)(clock()-start))/CLOCKS_PER_SEC;
	iterpersec=aveiter/elapsed;
	aveiter/=solcount;
	
	fp=fopen(statsfile,"a");
	fprintf(fp,"\n %d/%d solutions%10.2e iterations/solution\n",solcount,trials,aveiter);
	fprintf(fp,"%10.2lf iterations/sec\n",iterpersec);
	fclose(fp);
	
	return 0;
	}
