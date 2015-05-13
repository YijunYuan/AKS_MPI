#include <mpi.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/arith.h>
#include <mpfr.h>
#include <mpir.h>
#include <math.h>
#include <stdio.h>
//#define DEBUG
/*
__thread static fmpz_t temp1;
static fmpz_t temp2;
static fmpz_t temp3;
static ulong  temp_ui;*/
static fmpz_t tmp_td1_1;
static __inline__ void MULTIPLICATIVE_ORDER(fmpz_t out,const fmpz_t n,const fmpz_t mod){
    fmpz_gcd(out,n,mod);
    if(!fmpz_equal_ui(out,1)){
        fmpz_set_si(out,-1);return ;
    }
    fmpz_t i;fmpz_init(i);
    for(fmpz_one(i);;fmpz_add_ui(i,i,1)){
        fmpz_powm(out,mod,i,n);
        if(fmpz_equal_ui(out,1)){
                fmpz_set(out,i);
                fmpz_clear(i);
                return ;
        }
    }
}*/
/*int ui_mul_order_larger(ulong r,fmpz_t n_t,ulong bound){
    ulong n_ui;fmpz_mod_ui(tmp_td1_1,n_t,r);n_ui=fmpz_get_ui(tmp_td1_1);
    ulong count,res;
    res=1;count=0;
    while(1){
        res=(res*n_ui)%r;
        count++;
        if(count==bound)
            return 1;
        if(res==1)
    }
}

int fmpz_mul_order_larger(fmpz_t r,fmpz_t n_t,ulong bound){
    fmpz_set(temp2,n_t);fmpz_mod(temp2,n_t,r);
    fmpz_gcd(temp1,temp2,r);if(fmpz_cmp_ui(temp1,1)>0)return 0;
    fmpz_zero(temp1);fmpz_set_ui(temp3,1);
    while(1){
        fmpz_mul(temp3,temp3,temp2);fmpz_mod(temp3,temp3,r);fmpz_add_ui(temp1,temp1,1);
        if(fmpz_equal_ui(temp3,1)){
            if(fmpz_cmp_ui(temp1,bound)<=0)
                return 0;
            else
                return 1;
        }
    }

}
*/
ulong ui_log_2_n_sqr(fmpz_t& n,mpz_t& n_z){///The function is safe in the considerate range
    mpfr_prec_t prec=512;
    mpz_t log2_n_sqr_u,log2_n_sqr_d;mpz_init(log2_n_sqr_d);mpz_init(log2_n_sqr_u);
    mpfr_t n_f;mpfr_init2(n_f,prec);mpfr_set_z(n_f,n_z,MPFR_RNDN);
    mpfr_t log2_n_f_u,log2_n_f_d,temp_f;mpfr_init2(log2_n_f_d,prec);mpfr_init2(log2_n_f_u,prec);mpfr_init2(temp_f,prec);
    mpfr_log2(log2_n_f_u,n_f,MPFR_RNDU);mpfr_sqr(temp_f,log2_n_f_u,MPFR_RNDU);mpfr_get_z(log2_n_sqr_u,temp_f,MPFR_RNDU);
    mpfr_log2(log2_n_f_d,n_f,MPFR_RNDD);mpfr_sqr(temp_f,log2_n_f_d,MPFR_RNDD);mpfr_get_z(log2_n_sqr_d,temp_f,MPFR_RNDU);
    while(mpz_cmp(log2_n_sqr_d,log2_n_sqr_u)){//Ensure the precision is enough, add 64-bits each time
        if(prec>MPFR_PREC_MAX-257){
            printf("Precision is not enough when calculating. Abort.\n");
            MPI_Abort(MPI_COMM_WORLD,-1);
        }
        prec+=256;
        mpfr_prec_round(log2_n_f_u,prec,MPFR_RNDU);
        mpfr_prec_round(log2_n_f_d,prec,MPFR_RNDU);
        mpfr_prec_round(temp_f    ,prec,MPFR_RNDU);
        mpfr_log2(log2_n_f_u,n_f,MPFR_RNDU);mpfr_sqr(temp_f,log2_n_f_u,MPFR_RNDU);mpfr_get_z(log2_n_sqr_u,temp_f,MPFR_RNDU);
        mpfr_log2(log2_n_f_d,n_f,MPFR_RNDD);mpfr_sqr(temp_f,log2_n_f_d,MPFR_RNDD);mpfr_get_z(log2_n_sqr_d,temp_f,MPFR_RNDU);
    }
    return mpz_get_ui(log2_n_sqr_u);
}

int main(int argc,char* argv[]){

    int myid,numprocs;
    double start_time,end_time;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    printf("This is thread %d of %d in total\n",myid,numprocs);
    if(myid==0){
        //fmpz_init(tmp_td1_1);
        fmpz_init(temp1);
        ulong temp_ui;
        /*
        fmpz_init(temp2);
        fmpz_init(temp3);*/
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
            //fmpz_init(tmp_td1_1);
            fmpz_clear(temp1);/*
            fmpz_clear(temp2);
            fmpz_clear(temp3);*/
            fmpz_clear(n);
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        /*Step 2-- Find the r*/
        ulong r_ui=0,res,count;
        ulong log2_n_sqr_ui=ui_log_2_n_sqr(n,n_z);
        if(fmpz_bits(n)<7001){//r can be contained in a 64-bit ulong type safely
            for(r_ui=log2_n_sqr_ui;;r_ui++){
                count=0;res=1;
                temp_ui=fmpz_fdiv_ui(n,r_ui);
                if(n_gcd(r_ui,temp_ui)>1)continue;
                while(1){
                    res=(tes*temp_ui)%r_ui;
                    count++;
                    if(res==1&&count<=log2_n_sqr_ui)
                        break;
                    if(count==log2_n_sqr_ui)
                        goto ext;

                }
            }
            ext:;

        }
        else{
            for(fmpz_set_ui(r,log2_n_sqr_ui+1);;fmpz_add_ui(r,r,1)){
                if(fmpz_mul_order_larger(r,n,log2_n_sqr_ui)){
                    if(fmpz_cmp_ui(r,18446744073709551615)>=0){
                        printf("r overflow, abort.\n");
                        MPI_Abort(MPI_COMM_WORLD,-2);
                    }
                    r_ui=fmpz_get_ui(r);
                    break;
                }
            }
        }
#ifdef DEBUG
        ulong r_ui=0;
        if(fmpz_bits(n)<7001){//r can be contained in a 64-bit ulong type safely
            for(r_ui=log2_n_sqr_ui+1;;r_ui++){
                if(ui_mul_order_larger(r_ui,n,log2_n_sqr_ui))
                    break;
            }
        }
        else{
            for(fmpz_set_ui(r,log2_n_sqr_ui+1);;fmpz_add_ui(r,r,1)){
                if(fmpz_mul_order_larger(r,n,log2_n_sqr_ui)){
                    if(fmpz_cmp_ui(r,18446744073709551615)>=0){
                        printf("r overflow, abort.\n");
                        MPI_Abort(MPI_COMM_WORLD,-2);
                    }
                    r_ui=fmpz_get_ui(r);
                    break;
                }
            }
        }

        /*Step 3*/
        #pragma omp parallel for private(temp_ui),private(temp1)
        for(temp_ui=2;temp_ui<=r;temp_ui++){
            temp_ui=n_gcd(fmpz_mod_ui(n,temp_ui));
            if(temp_ui>1){
                #pragma omp critical
                {
                    end_time=MPI_Wtime();
                    printf("\n\nResult:\t0\nTotal time:\t%.12lf\n\n",end_time-start_time);
                    MPI_Abort(MPI_COMM_WORLD,1);
                }
            }
        }

        /*Step 4*/
        if(fmpz_cmp_ui(n,r_ui)<0){
            end_time=MPI_Wtime();
            printf("\n\nResult:\t1\nTotal time:\t%.12lf\n\n",end_time-start_time);
            MPI_Abort(MPI_COMM_WORLD,2);
        }
#endif // DEBUG
    }else{

    }
    MPI_Finalize();
}

