#include <mpi.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/arith.h>
#include <mpfr.h>
#include <mpir.h>
#include <math.h>
#include <stdio.h>

#define MPFR_DEFAULT_PREC 512
//#define EXPERIMENTAL

ulong ui_log_2_n_sqr(mpz_t n_z){///The function is safe in the considerate range
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

    mpfr_sqr  (temp_f  ,r_elr ,MPFR_RNDU);
    mpfr_log2 (log2_n_u,n_f   ,MPFR_RNDU);
    mpfr_mul  (temp_f  ,temp_f,log2_n_u ,MPFR_RNDU);
    mpfr_get_z(res_u   ,temp_f,MPFR_RNDD);

    mpfr_sqr  (temp_f  ,r_elr ,MPFR_RNDD);
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

        mpfr_sqr  (temp_f  ,r_elr ,MPFR_RNDU);
        mpfr_log2 (log2_n_u,n_f   ,MPFR_RNDU);
        mpfr_mul  (temp_f  ,temp_f,log2_n_u ,MPFR_RNDU);
        mpfr_get_z(res_u   ,temp_f,MPFR_RNDD);

        mpfr_sqr  (temp_f  ,r_elr ,MPFR_RNDD);
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

    int myid,numprocs;
    double start_time,end_time;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    printf("This is thread %d of %d in total\n",myid,numprocs);
    if(myid==0){
        ulong temp_ui;
        /*Read the test data from file*/
        printf("Reading the test data from file...\n");
        fmpz_t n;fmpz_init(n);
        FILE* data=fopen("aks.data","r");
        fmpz_fread (data, n);
        fclose(data);

        start_time=MPI_Wtime();

        /*Step 1-- Check the prefect power*/
        mpz_t n_z;mpz_init(n_z);fmpz_get_mpz(n_z,n);
        if(mpz_perfect_power_p(n_z)!=0){
            end_time=MPI_Wtime();
            printf("\n\nResult:\t0\nTotal time:\t%.12lf\n\n",end_time-start_time);
            mpz_clear(n_z);
            fmpz_clear(n);
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        /*Step 2-- Find the r*/
        slong r_si=0;
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
            if(n_gcd(a,temp_ui)>1){
                end_time=MPI_Wtime();
                printf("\n\nResult:\t0\nTotal time:\t%.12lf\n\n",end_time-start_time);
                mpz_clear(n_z);
                fmpz_clear(n);
                MPI_Abort(MPI_COMM_WORLD,1);
            }
        }

        /*Step 4*/
        if(fmpz_cmp_ui(n,5690034)<=0&&fmpz_cmp_si(n,r_si)<=0){
            end_time=MPI_Wtime();
            printf("\n\nResult:\t1\nTotal time:\t%.12lf\n\n",end_time-start_time);
            mpz_clear(n_z);
            fmpz_clear(n);
            MPI_Abort(MPI_COMM_WORLD,0);
        }


    }else{
;
    }
    MPI_Finalize();
}

