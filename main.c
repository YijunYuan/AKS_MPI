#include <mpi.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/arith.h>
#include <mpfr.h>
#include <mpir.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#define MPFR_DEFAULT_PREC 512
#define MAX_STRING_LENGTH 7000
//#define EXPERIMENTAL

void write_log(int step,fmpz_t n,double time,ulong r,ulong a,ulong b,int numprocs){
    char name[MAX_STRING_LENGTH];fmpz_get_str(name,16,n);
    strcat(name,".log");
    FILE* log=fopen(name,"wr");
    fprintf(log,"This is the test log\n\n");
    fprintf(log,"n=\t%s\n",fmpz_get_str(NULL,10,n));
    fprintf(log,"Time used:\t%lf\n",time);
    if(step==1){
        fprintf(log,"Test result:\t0\n");
        fprintf(log,"The test failed in Step\t1\tn is a power.\n");
        fclose(log);return NULL;
    }
    fprintf(log,"Step 1 pass!\n");
    if(step==2){
        fprintf(log,"Test abort. r can't be calculated\n");
        return NULL;
    }
    fprintf(log,"Step 2 pass!\n");
    fprintf(log,"r=\t%llu\n",r);
    if(step==3){
        fprintf(log,"Test result:\t0\n");
        fprintf(log,"The test failed in Step\t3\tGCD of n and %llu is greater than 1.\n",a);
        fclose(log);return NULL;
    }
    fprintf(log,"Step 3 pass!\n");
    if(step==4){
        fprintf(log,"Test result:\t1\n");
        fprintf(log,"The test finished in Step\t4, n<=r\n");
        fclose(log);return NULL;
    }
    fprintf(log,"Step 4 pass!\n");
    fprintf(log,"Parallel part is activated.\n%d threads are involved.\n",numprocs);
    fprintf(log,"%llu tasks in total\n",a);
    if(step==5){
        fprintf(log,"Test result:\t0\n");
        fprintf(log,"The test failded when a=\t%llu\n",b);
        fclose(log);return NULL;
    }
    fprintf(log,"Step 5 pass!\n");
    fprintf(log,"All test pass! n is a prime.\n");
    fclose(log);return NULL;
}

ulong ui_log_2_n_sqr(mpz_t n_z){
#ifndef EXPERIMENTAL
    mpfr_prec_t prec=MPFR_DEFAULT_PREC;

    mpz_t res_u,res_d;
    mpz_init(res_d);
    mpz_init(res_u);

    mpfr_t n_f;
    mpfr_init2(n_f,prec);
    mpfr_set_z(n_f,n_z,MPFR_RNDN);

    mpfr_t log2_n_f_u,log2_n_f_d,temp_f;
    mpfr_init2(log2_n_f_d,prec);
    mpfr_init2(log2_n_f_u,prec);
    mpfr_init2(temp_f    ,prec);

    mpfr_log2 (log2_n_f_u,n_f       ,MPFR_RNDU);
    mpfr_sqr  (temp_f    ,log2_n_f_u,MPFR_RNDU);
    mpfr_get_z(res_u     ,temp_f    ,MPFR_RNDU);

    mpfr_log2 (log2_n_f_d,n_f       ,MPFR_RNDD);
    mpfr_sqr  (temp_f    ,log2_n_f_d,MPFR_RNDD);
    mpfr_get_z(res_d     ,temp_f    ,MPFR_RNDU);

    while(mpz_cmp(res_d,res_u)){
        if(prec>MPFR_PREC_MAX-257){
            printf("Precision is not enough when calculating. Abort.\n");
            MPI_Abort(MPI_COMM_WORLD,-1);
        }
        prec+=256;
        mpfr_prec_round(log2_n_f_u,prec,MPFR_RNDU);
        mpfr_prec_round(log2_n_f_d,prec,MPFR_RNDD);
        mpfr_prec_round(temp_f    ,prec,MPFR_RNDN);

        mpfr_log2 (log2_n_f_u,n_f       ,MPFR_RNDU);
        mpfr_sqr  (temp_f    ,log2_n_f_u,MPFR_RNDU);
        mpfr_get_z(res_u     ,temp_f    ,MPFR_RNDU);

        mpfr_log2 (log2_n_f_d,n_f       ,MPFR_RNDD);
        mpfr_sqr  (temp_f    ,log2_n_f_d,MPFR_RNDD);
        mpfr_get_z(res_d     ,temp_f    ,MPFR_RNDU);
    }

    mpfr_clear(n_f       );
    mpfr_clear(log2_n_f_u);
    mpfr_clear(log2_n_f_d);
    mpfr_clear(temp_f    );
    mpz_clear (res_d     );

    return mpz_get_ui(res_u);
#else
    mpfr_t q;
    mpfr_init2(q,1024);
    mpfr_set_z(q,n_z,MPFR_RNDN);
    mpfr_log2(q,q,MPFR_RNDN);
    mpfr_sqr(q,q,MPFR_RNDU);
    return ceill(mpfr_get_ld(q,MPFR_RNDU));
#endif
}

ulong calc_pama(slong r_si,mpz_t n_z){
    r_si=n_euler_phi(r_si);
#ifndef EXPERIMENTAL
    mpfr_prec_t prec=MPFR_DEFAULT_PREC;

    mpfr_t r_elr,n_f;

    mpfr_init2(r_elr,prec);
    mpfr_set_ui(r_elr,r_si,MPFR_RNDN);

    mpfr_init2(n_f,prec);
    mpfr_set_z(n_f,n_z,MPFR_RNDN);

    mpfr_t log2_n_u,log2_n_d,temp_f;
    mpfr_init2(log2_n_u,prec);
    mpfr_init2(log2_n_d,prec);
    mpfr_init2(temp_f,prec);

    mpz_t res_u,res_d;
    mpz_init(res_u);mpz_init(res_d);

    mpfr_sqrt  (temp_f  ,r_elr ,MPFR_RNDU);
    mpfr_log2 (log2_n_u,n_f   ,MPFR_RNDU);
    mpfr_mul  (temp_f  ,temp_f,log2_n_u ,MPFR_RNDU);
    mpfr_get_z(res_u   ,temp_f,MPFR_RNDD);

    mpfr_sqrt  (temp_f  ,r_elr ,MPFR_RNDD);
    mpfr_log2 (log2_n_d,n_f   ,MPFR_RNDD);
    mpfr_mul  (temp_f  ,temp_f,log2_n_d ,MPFR_RNDD);
    mpfr_get_z(res_d   ,temp_f,MPFR_RNDD);

    while(mpz_cmp(res_d,res_u)){
        if(prec>MPFR_PREC_MAX-257){
            printf("Precision is not enough when calculating. Abort.\n");
            MPI_Abort(MPI_COMM_WORLD,-1);
        }
        prec+=256;
        mpfr_prec_round(log2_n_u,prec,MPFR_RNDU);
        mpfr_prec_round(log2_n_d,prec,MPFR_RNDD);
        mpfr_prec_round(temp_f,prec,MPFR_RNDN);

        mpfr_sqrt  (temp_f  ,r_elr ,MPFR_RNDU);
        mpfr_log2 (log2_n_u,n_f   ,MPFR_RNDU);
        mpfr_mul  (temp_f  ,temp_f,log2_n_u ,MPFR_RNDU);
        mpfr_get_z(res_u   ,temp_f,MPFR_RNDD);

        mpfr_sqrt  (temp_f  ,r_elr ,MPFR_RNDD);
        mpfr_log2 (log2_n_d,n_f   ,MPFR_RNDD);
        mpfr_mul  (temp_f  ,temp_f,log2_n_d ,MPFR_RNDD);
        mpfr_get_z(res_d   ,temp_f,MPFR_RNDD);
    }

    mpfr_clear(r_elr   );
    mpfr_clear(n_f     );
    mpfr_clear(log2_n_u);
    mpfr_clear(log2_n_d);
    mpfr_clear(temp_f  );
    mpz_clear (res_d   );

    return mpz_get_ui(res_u);

#else
    long double res;
    res=sqrtl(r_si);
    mpfr_t t_f;
    mpfr_init2(t_f,1024);mpfr_init_set_z(t_f,n_z,MPFR_RNDN);
    mpfr_log2(t_f,t_f,MPFR_RNDN);res=mpfr_get_ld(t_f,MPFR_RNDN);
    res*=sqrtl(r_si);
    mpfr_clear(t_f);
    return floorl(res);
#endif // EXPERIMENTAL
}

int main(int argc,char* argv[]){
    char* n_str;
    ulong pama;
    slong r_si;
    slong n_mod_r;
    int myid,numprocs;
    double start_time=0,end_time=0;
    _Bool tag=true;
    fmpz_t n;fmpz_init(n);

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    printf("This is thread %d of %d in total\n",myid,numprocs);
    if(myid==0){
        ulong temp_ui=0;
        /*Read the test data from file*/
        printf("Reading the test data from file...\n");
        FILE* data=fopen("aks.data","rw");
        fmpz_fread (data, n);
        fclose(data);
        start_time=MPI_Wtime();
        /*Step 1-- Check the prefect power*/
        mpz_t n_z;mpz_init(n_z);fmpz_get_mpz(n_z,n);
        if(mpz_perfect_power_p(n_z)!=0){
            end_time=MPI_Wtime();
            //printf("\n\nResult:\t0\nTotal time:\t%.12lf\n\n",end_time-start_time);
            write_log(1,n,end_time-start_time,0,0,0,numprocs);
            printf("\nnot prime. Details are in the log.\n");
            mpz_clear(n_z);
            fmpz_clear(n);
            MPI_Abort(MPI_COMM_WORLD,1);
        }

        /*Step 2-- Find the r*/
        r_si=0;
        ulong res_ui=ui_log_2_n_sqr(n_z);
        if(fmpz_bits(n)<6205){//r can be contained in a 63-bit slong type safely
            __uint128_t res;///This marco is only supported in the x86_64 GCC, to ensure the safety
            ulong count;
            for(r_si=res_ui;;r_si++){
                count=0;res=1;
                temp_ui=fmpz_fdiv_ui(n,r_si);
                if(n_gcd(r_si,temp_ui)>1)continue;
                while(1){
                    res=(res*temp_ui)%r_si;
                    count++;
                    if(res==1&&count<=res_ui)
                        break;
                    if(count==res_ui)
                        goto ext;
                }
            }
            ext:;

        }
        else{
            ///This blocks is still under typing!!!
            printf("Not completed!\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        /*Step 3*/
        slong a=2;
        for(a=2;a<=r_si;a++){
            temp_ui=fmpz_fdiv_ui(n,a);
            if(n_gcd(a,temp_ui)>1){
                end_time=MPI_Wtime();
                write_log(3,n,end_time-start_time,r_si,a,0,numprocs);
                printf("\nnot prime. Details are in the log.\n");
                mpz_clear(n_z);
                fmpz_clear(n);
                MPI_Abort(MPI_COMM_WORLD,1);
            }
        }

        /*Step 4*/
        if(fmpz_cmp_ui(n,5690034)<=0&&fmpz_cmp_si(n,r_si)<=0){
            end_time=MPI_Wtime();
            write_log(4,n,end_time-start_time,r_si,0,0,numprocs);
            printf("\nprime. Details are in the log.\n");
            mpz_clear(n_z);
            fmpz_clear(n);
            MPI_Abort(MPI_COMM_WORLD,0);
        }

        /*Step 5*/
        printf("r=%lld\n",r_si);
        pama=calc_pama(r_si,n_z);
        printf("pama=%llu\n",pama);
        getchar();
        n_str=fmpz_get_str(NULL,10,n);
        n_mod_r=fmpz_fdiv_ui(n,r_si);
        //fmpz_clear(n);
        mpz_clear(n_z);

    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&pama   ,1              ,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
    MPI_Bcast(&r_si   ,1              ,MPI_LONG_LONG         ,0,MPI_COMM_WORLD);
    MPI_Bcast(&n_mod_r,1              ,MPI_LONG_LONG         ,0,MPI_COMM_WORLD);

    int sl=0;
    if(myid==0)
        sl=strlen(n_str);
    MPI_Bcast(&sl,1,MPI_INT,0,MPI_COMM_WORLD);
    if(myid!=0)n_str=(char*)malloc(sizeof(char)*sl);
    MPI_Bcast(n_str   ,sl+1,MPI_CHAR              ,0,MPI_COMM_WORLD);

    if(myid==0){
        ulong test_bad_point_a;
        MPI_Status sta;
        long i;
        for(i=1;i<numprocs;i++){
            MPI_Recv(&tag     ,1,MPI_C_BOOL,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&sta);
            MPI_Recv(&end_time,1,MPI_DOUBLE,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&sta);
            if(tag==false){
                MPI_Recv(&test_bad_point_a,1,MPI_UNSIGNED_LONG_LONG,MPI_ANY_SOURCE,3,MPI_COMM_WORLD,&sta);
                break;
            }
        }
        if(tag==false){
            write_log(5,n,end_time-start_time,r_si,pama,test_bad_point_a,numprocs);
            printf("\nnot prime. Details are in the log.\n");
            MPI_Abort(MPI_COMM_WORLD,0);
        }
        printf("\nprime. Details are in the log.\n");
        write_log(6,n,end_time-start_time,r_si,pama,0,numprocs);
        MPI_Finalize();
        return 0;

    }
    else{
        fmpz_set_str(n,n_str,10);
        fmpz_t n_1;fmpz_init(n_1);fmpz_set_si(n_1,-1);
        fmpz_mod_poly_t poly1,poly2,poly_tmp,mod_poly;
        fmpz_mod_poly_init(poly1   ,n);
        fmpz_mod_poly_init(poly2   ,n);
        fmpz_mod_poly_init(poly_tmp,n);

        fmpz_mod_poly_set_coeff_ui(poly1,1, 1);

        fmpz_mod_poly_set_coeff_ui(poly2,n_mod_r,1);
        fmpz_mod_poly_init(mod_poly,n);
        fmpz_mod_poly_set_coeff_ui(mod_poly,r_si, 1);
        fmpz_mod_poly_set_coeff_fmpz(mod_poly,0   ,n_1);

        ulong i;
        for(i=myid-1;i<=pama;i+=numprocs){
            printf("%llu\n",i);
            fmpz_mod_poly_set_coeff_ui(poly1,0,i);
            fmpz_mod_poly_powmod_fmpz_binexp(poly_tmp,poly1,n,mod_poly);
            fmpz_mod_poly_set_coeff_ui(poly2,0,i);
            if(!fmpz_mod_poly_equal(poly_tmp,poly2)){
                tag=false;
                break;
            }
        }
        end_time=MPI_Wtime();
        MPI_Send(&tag     ,1,MPI_C_BOOL            ,0,1,MPI_COMM_WORLD);
        MPI_Send(&end_time,1,MPI_DOUBLE            ,0,2,MPI_COMM_WORLD);
        if(tag==false)
            MPI_Send(&i       ,1,MPI_UNSIGNED_LONG_LONG,0,3,MPI_COMM_WORLD);
        fmpz_clear(n);
        fmpz_mod_poly_clear(poly1   );
        fmpz_mod_poly_clear(poly2   );
        fmpz_mod_poly_clear(poly_tmp);
        fmpz_mod_poly_clear(mod_poly);
    }
    MPI_Finalize();
    return 0;
}

