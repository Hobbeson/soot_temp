#include "mass_routines.h"
#include <cstring>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <algorithm>

// =============================================================================
// mass 模块 C++ 完整实现
// =============================================================================

namespace {


static double pbar_at(const double* PBAR_RAW, int PBAR_NK, int K, int PZ) {
    return PBAR_RAW[K + PBAR_NK * PZ];
}

struct MassB {
    int ibar,jbar,kbar,ibp1,jbp1,kbp1,ibm1,jbm1,kbm1;
    int nts,n_total,n_zone,n_passive,zeta_idx;
    int nw_ext,nw_int,nw;
    bool predictor,first_pass,cc_ibm,store_species_flux;
    bool flux_limiter_mw_correction,solid_phase_only;
    bool grav_settling,ther_settling,have_m_dot_ppp;
    int i_flux_limiter,period_test,icyc;
    double rhomin,rhomax;
    IX3 i3,i3p1;
    IX4 i4_zz,i4_fx,i4_adv,i4_drdz,i4_mdot;
    double *RHO,*RHOS,*TMP,*RSUM,*ZZ,*ZZS;
    double *DEL_RHO_D_DEL_Z,*FX,*FY,*FZ,*M_DOT_PPP;
    double *ADV_FX,*ADV_FY,*ADV_FZ;
    const double *U,*V,*W,*US,*VS,*WS;
    const double *R,*RRN,*DX,*RDX,*DY,*RDY,*DZ,*RDZ,*XC,*ZC;
    double *WORK1,*WORK2,*WORK3,*WORK4,*WORK5,*SWORK4;
    double *D_SOURCE;
    double *PBAR,*PBAR_S;
    const int *PRESSURE_ZONE;
    const double *D_PBAR_DT,*D_PBAR_DT_S;
    const int *CELL_INDEX,*CELL_WALL_INDEX;
    const double *CELL_SOLID,*CELL_EXTERIOR,*UVW_SAVE;
    const int *WALL_BC_TYPE,*WALL_B1_INDEX;
    const int *BC_IIG,*BC_JJG,*BC_KKG,*BC_IOR;
    const int *BC_II,*BC_JJ,*BC_KK;
    const double *BP1_RHO_F,*BP1_ZZ_F;
    const int *EWC_BOUNDARY_TYPE_PREVIOUS;
    const double *SM_MW,*SM_DEPOSITING;
    inline int ci(int I,int J,int K)const{return CELL_INDEX[i3p1(I,J,K)]-1;}
    inline bool iss(int I,int J,int K)const{int idx=ci(I,J,K);return idx>=0&&CELL_SOLID&&CELL_SOLID[idx]!=0.0;}
    inline bool depo()const{for(int i=0;i<nts;++i)if(SM_DEPOSITING[i]!=0.0)return true;return false;}

    MassB(const MassMeshPointersC* c){
        ibar=c->IBAR;jbar=c->JBAR;kbar=c->KBAR;
        ibp1=ibar+1;jbp1=jbar+1;kbp1=kbar+1;
        ibm1=ibar-1;jbm1=jbar-1;kbm1=kbar-1;
        nts=c->N_TRACKED_SPECIES;n_total=c->N_TOTAL_SCALARS;
        n_zone=c->N_ZONE;n_passive=c->N_PASSIVE_SCALARS;zeta_idx=c->ZETA_INDEX;
        nw_ext=c->N_EXTERNAL_WALL_CELLS;nw_int=c->N_INTERNAL_WALL_CELLS;nw=nw_ext+nw_int;
        predictor=c->PREDICTOR_FLAG!=0;first_pass=c->FIRST_PASS_FLAG!=0;
        cc_ibm=c->CC_IBM_FLAG!=0;store_species_flux=c->STORE_SPECIES_FLUX_FLAG!=0;
        flux_limiter_mw_correction=c->FLUX_LIMITER_MW_CORRECTION_FLAG!=0;
        solid_phase_only=c->SOLID_PHASE_ONLY_FLAG!=0;
        grav_settling=c->GRAVITATIONAL_SETTLING_FLAG!=0;
        ther_settling=c->THERMOPHORETIC_SETTLING_FLAG!=0;
        have_m_dot_ppp=c->M_DOT_PPP!=nullptr;
        i_flux_limiter=c->I_FLUX_LIMITER;period_test=c->PERIODIC_TEST_VAL;icyc=c->ICYC_VAL;
        rhomin=c->RHOMIN_VAL;rhomax=c->RHOMAX_VAL;
        i3=IX3(ibar+1,jbar+1,kbar+1);i3p1=IX3(ibp1+1,jbp1+1,kbp1+1);

        i4_zz=IX4(ibp1+1,jbp1+1,kbp1+1,n_total,1);
        i4_fx=IX4(ibp1+1,jbp1+1,kbp1+1,n_total+(flux_limiter_mw_correction?1:0),flux_limiter_mw_correction?0:1);
        i4_adv=IX4(ibp1+1,jbp1+1,kbp1+1,n_total,1);
        i4_drdz=IX4(ibp1+1,jbp1+1,kbp1+1,n_total,1);
        i4_mdot=IX4(ibp1+1,jbp1+1,kbp1+1,nts,1);
        auto d=[](void*p)->double*{return static_cast<double*>(p);};
        auto cd=[](const void*p)->const double*{return static_cast<const double*>(p);};
        auto ci2=[](const void*p)->const int*{return static_cast<const int*>(p);};
        RHO=d(c->RHO);RHOS=d(c->RHOS);TMP=d(c->TMP);RSUM=d(c->RSUM);
        ZZ=d(c->ZZ);ZZS=d(c->ZZS);DEL_RHO_D_DEL_Z=d(c->DEL_RHO_D_DEL_Z);
        FX=d(c->FX);FY=d(c->FY);FZ=d(c->FZ);
        M_DOT_PPP=d(c->M_DOT_PPP);
        ADV_FX=d(c->ADV_FX);ADV_FY=d(c->ADV_FY);ADV_FZ=d(c->ADV_FZ);
        U=cd(c->U);V=cd(c->V);W=cd(c->W);
        US=cd(c->US);VS=cd(c->VS);WS=cd(c->WS);
        R=cd(c->R);RRN=cd(c->RRN);
        DX=cd(c->DX);RDX=cd(c->RDX);
        DY=cd(c->DY);RDY=cd(c->RDY);
        DZ=cd(c->DZ);RDZ=cd(c->RDZ);
        XC=cd(c->XC);ZC=cd(c->ZC);
        WORK1=d(c->WORK1);WORK2=d(c->WORK2);WORK3=d(c->WORK3);
        WORK4=d(c->WORK4);WORK5=d(c->WORK5);SWORK4=d(c->SWORK4);
        D_SOURCE=d(c->D_SOURCE);
        PBAR=d(c->PBAR);PBAR_S=d(c->PBAR_S);
        PRESSURE_ZONE=ci2(c->PRESSURE_ZONE);
        D_PBAR_DT=cd(c->D_PBAR_DT);D_PBAR_DT_S=cd(c->D_PBAR_DT_S);
        CELL_INDEX=ci2(c->CELL_INDEX);CELL_WALL_INDEX=ci2(c->CELL_WALL_INDEX);
        CELL_SOLID=cd(c->CELL_SOLID);CELL_EXTERIOR=cd(c->CELL_EXTERIOR);
        UVW_SAVE=cd(c->UVW_SAVE);
        WALL_BC_TYPE=ci2(c->WALL_BC_TYPE);WALL_B1_INDEX=ci2(c->WALL_B1_INDEX);
        BC_IIG=ci2(c->BC_IIG);BC_JJG=ci2(c->BC_JJG);
        BC_KKG=ci2(c->BC_KKG);BC_IOR=ci2(c->BC_IOR);
        BC_II=ci2(c->BC_II);BC_JJ=ci2(c->BC_JJ);BC_KK=ci2(c->BC_KK);
        BP1_RHO_F=cd(c->BP1_RHO_F);BP1_ZZ_F=cd(c->BP1_ZZ_F);
        EWC_BOUNDARY_TYPE_PREVIOUS=ci2(c->EWC_BOUNDARY_TYPE_PREVIOUS);
        SM_MW=cd(c->SM_MW);SM_DEPOSITING=cd(c->SM_DEPOSITING);
    }
};

// =============================================================================
// MASS_FINITE_DIFFERENCES (行 20-333)
// =============================================================================
void mass_fd_impl(const MassMeshPointersC* c){
    {
        static FILE* fp=std::fopen("/tmp/mass_c_fd_called.txt","w");
        if(fp){std::fprintf(fp,"mass_fd_impl called\n");std::fflush(fp);}
    }
    MassB m(c); if(m.solid_phase_only) return;
    int ibar=m.ibar,jbar=m.jbar,kbar=m.kbar;
    int ibp1=m.ibp1,jbp1=m.jbp1,kbp1=m.kbp1;
    int ibm1=m.ibm1,jbm1=m.jbm1,kbm1=m.kbm1;
    int nts=m.nts,n_total=m.n_total;
    const double* UU=m.predictor?m.U:m.US;
    const double* VV=m.predictor?m.V:m.VS;
    const double* WW=m.predictor?m.W:m.WS;
    const double* RHOP=m.predictor?m.RHO:m.RHOS;
    const double* ZZP=m.predictor?m.ZZ:m.ZZS;

    for(int N=1;N<=n_total;++N){
        double* RHO_Z_P=m.WORK1;
        // RHO_Z_P=RHOP*ZZP (行 67-75)
        #pragma omp parallel for
        for(int K=0;K<=kbp1;++K) for(int J=0;J<=jbp1;++J) for(int I=0;I<=ibp1;++I){
            int idx=m.i3p1(I,J,K);
            RHO_Z_P[idx]=RHOP[idx]*ZZP[m.i4_zz(I,J,K,N)];
        }
        // GET_SCALAR_FACE_VALUE (行 79-81)
        fds_get_scalar_face_value(UU,RHO_Z_P,&m.FX[m.i4_fx(0,0,0,N)],1,ibm1,1,jbar,1,kbar,1,m.i_flux_limiter,ibp1+1,jbp1+1,kbp1+1);
        fds_get_scalar_face_value(VV,RHO_Z_P,&m.FY[m.i4_fx(0,0,0,N)],1,ibar,1,jbm1,1,kbar,2,m.i_flux_limiter,ibp1+1,jbp1+1,kbp1+1);
        fds_get_scalar_face_value(WW,RHO_Z_P,&m.FZ[m.i4_fx(0,0,0,N)],1,ibar,1,jbar,1,kbm1,3,m.i_flux_limiter,ibp1+1,jbp1+1,kbp1+1);
        // WALL_LOOP_2 (行 83-177)
        #pragma omp parallel for
        for(int IW=0;IW<m.nw;++IW){
            int bt=m.WALL_BC_TYPE[IW]; if(bt==NULL_BOUNDARY) continue;
            int II=m.BC_II?m.BC_II[IW]:0, JJ=m.BC_JJ?m.BC_JJ[IW]:0, KK=m.BC_KK?m.BC_KK[IW]:0;
            int IIG=m.BC_IIG[IW],JJG=m.BC_JJG[IW],KKG=m.BC_KKG[IW],IOR=m.BC_IOR[IW];
            // IC = CELL_INDEX(II,JJ,KK)-1 (ghost cell index, 0-based)
            int IC=-1;
            if (m.BC_II && m.BC_JJ && m.BC_KK) IC=m.ci(II,JJ,KK);

            if(bt==SOLID_BOUNDARY && IC>=0 && m.CELL_SOLID && m.CELL_EXTERIOR &&
               m.CELL_SOLID[IC]==0.0 && m.CELL_EXTERIOR[IC]==0.0){
                // 内部固壁
                switch(IOR){
                    case 1: m.FX[m.i4_fx(IIG-1,JJG,KKG,N)]=0.0;break;
                    case-1: m.FX[m.i4_fx(IIG,JJG,KKG,N)]=0.0;break;
                    case 2: m.FY[m.i4_fx(IIG,JJG-1,KKG,N)]=0.0;break;
                    case-2: m.FY[m.i4_fx(IIG,JJG,KKG,N)]=0.0;break;
                    case 3: m.FZ[m.i4_fx(IIG,JJG,KKG-1,N)]=0.0;break;
                    case-3: m.FZ[m.i4_fx(IIG,JJG,KKG,N)]=0.0;break;
                }
            } else {
                // 其他边界
                double fv=0.0;
                if(m.BP1_RHO_F) fv=m.BP1_RHO_F[IW];
                if(m.BP1_ZZ_F)  fv*=m.BP1_ZZ_F[IW*nts+(N-1)];
                switch(IOR){
                    case 1: m.FX[m.i4_fx(IIG-1,JJG,KKG,N)]=fv;break;
                    case-1: m.FX[m.i4_fx(IIG,JJG,KKG,N)]=fv;break;
                    case 2: m.FY[m.i4_fx(IIG,JJG-1,KKG,N)]=fv;break;
                    case-2: m.FY[m.i4_fx(IIG,JJG,KKG,N)]=fv;break;
                    case 3: m.FZ[m.i4_fx(IIG,JJG,KKG-1,N)]=fv;break;
                    case-3: m.FZ[m.i4_fx(IIG,JJG,KKG,N)]=fv;break;
                }
            }

            if(bt!=INTERPOLATED_BOUNDARY && bt!=OPEN_BOUNDARY &&
               m.BC_II && m.BC_JJ && m.BC_KK && m.CELL_WALL_INDEX){
                
                double zt[64]={},ut[64]={},ft[64]={};
                switch(IOR){
                    case 1: // FX/UU(II+1)
                        if(UU[m.i3p1(II+1,JJ,KK)]>0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II+1,JJ,KK)*7)+4]<=0){
                            zt[20]=RHO_Z_P[m.i3p1(II+1,JJ,KK)];
                            zt[21]=RHO_Z_P[m.i3p1(II+1,JJ,KK)];
                            zt[22]=RHO_Z_P[m.i3p1(II+2,JJ,KK)];
                            zt[23]=0.0;
                            ut[21]=UU[m.i3p1(II+1,JJ,KK)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,1,m.i_flux_limiter,4,4,4);
                            m.FX[m.i4_fx(II+1,JJ,KK,N)]=ft[21];
                        }
                        break;
                    case -1:
                        if(UU[m.i3p1(II-2,JJ,KK)]<0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II-1,JJ,KK)*7)+2]<=0){
                            zt[20]=0.0;
                            zt[21]=RHO_Z_P[m.i3p1(II-2,JJ,KK)];
                            zt[22]=RHO_Z_P[m.i3p1(II-1,JJ,KK)];
                            zt[23]=RHO_Z_P[m.i3p1(II-1,JJ,KK)];
                            ut[21]=UU[m.i3p1(II-2,JJ,KK)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,1,m.i_flux_limiter,4,4,4);
                            m.FX[m.i4_fx(II-2,JJ,KK,N)]=ft[21];
                        }
                        break;
                    case 2:
                        if(VV[m.i3p1(II,JJ+1,KK)]>0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II,JJ+1,KK)*7)+5]<=0){
                            zt[17]=RHO_Z_P[m.i3p1(II,JJ+1,KK)];
                            zt[21]=RHO_Z_P[m.i3p1(II,JJ+1,KK)];
                            zt[25]=RHO_Z_P[m.i3p1(II,JJ+2,KK)];
                            zt[29]=0.0;
                            ut[21]=VV[m.i3p1(II,JJ+1,KK)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,2,m.i_flux_limiter,4,4,4);
                            m.FY[m.i4_fx(II,JJ+1,KK,N)]=ft[21];
                        }
                        break;
                    case -2:
                        if(VV[m.i3p1(II,JJ-2,KK)]<0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II,JJ-1,KK)*7)+1]<=0){
                            zt[17]=0.0;
                            zt[21]=RHO_Z_P[m.i3p1(II,JJ-2,KK)];
                            zt[25]=RHO_Z_P[m.i3p1(II,JJ-1,KK)];
                            zt[29]=RHO_Z_P[m.i3p1(II,JJ-1,KK)];
                            ut[21]=VV[m.i3p1(II,JJ-2,KK)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,2,m.i_flux_limiter,4,4,4);
                            m.FY[m.i4_fx(II,JJ-2,KK,N)]=ft[21];
                        }
                        break;
                    case 3:
                        if(WW[m.i3p1(II,JJ,KK+1)]>0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II,JJ,KK+1)*7)+6]<=0){
                            zt[5]=RHO_Z_P[m.i3p1(II,JJ,KK+1)];
                            zt[21]=RHO_Z_P[m.i3p1(II,JJ,KK+1)];
                            zt[37]=RHO_Z_P[m.i3p1(II,JJ,KK+2)];
                            zt[53]=0.0;
                            ut[21]=WW[m.i3p1(II,JJ,KK+1)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,3,m.i_flux_limiter,4,4,4);
                            m.FZ[m.i4_fx(II,JJ,KK+1,N)]=ft[21];
                        }
                        break;
                    case -3:
                        if(WW[m.i3p1(II,JJ,KK-2)]<0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II,JJ,KK-1)*7)+0]<=0){
                            zt[5]=0.0;
                            zt[21]=RHO_Z_P[m.i3p1(II,JJ,KK-2)];
                            zt[37]=RHO_Z_P[m.i3p1(II,JJ,KK-1)];
                            zt[53]=RHO_Z_P[m.i3p1(II,JJ,KK-1)];
                            ut[21]=WW[m.i3p1(II,JJ,KK-2)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,3,m.i_flux_limiter,4,4,4);
                            m.FZ[m.i4_fx(II,JJ,KK-2,N)]=ft[21];
                        }
                        break;
                }
            }
        }
    }
    // FACE_CORRECTION_IF
    if(m.flux_limiter_mw_correction){
        double* RHO_RMW=m.WORK1; double zz_get[512],mw_g;
        #pragma omp parallel for private(zz_get,mw_g)
        for(int K=0;K<=kbp1;++K) for(int J=0;J<=jbp1;++J) for(int I=0;I<=ibp1;++I){
            for(int ns=0;ns<nts;++ns) zz_get[ns]=ZZP[m.i4_zz(I,J,K,ns+1)];
            fds_get_molecular_weight(zz_get,&mw_g);
            RHO_RMW[m.i3p1(I,J,K)]=RHOP[m.i3p1(I,J,K)]/mw_g;
        }
        fds_get_scalar_face_value(UU,RHO_RMW,&m.FX[m.i4_fx(0,0,0,0)],1,ibm1,1,jbar,1,kbar,1,m.i_flux_limiter,ibp1+1,jbp1+1,kbp1+1);
        fds_get_scalar_face_value(VV,RHO_RMW,&m.FY[m.i4_fx(0,0,0,0)],1,ibar,1,jbm1,1,kbar,2,m.i_flux_limiter,ibp1+1,jbp1+1,kbp1+1);
        fds_get_scalar_face_value(WW,RHO_RMW,&m.FZ[m.i4_fx(0,0,0,0)],1,ibar,1,jbar,1,kbm1,3,m.i_flux_limiter,ibp1+1,jbp1+1,kbp1+1);
        // WALL_LOOP_3 
        #pragma omp parallel for
        for(int IW=0;IW<m.nw;++IW){
            int bt=m.WALL_BC_TYPE[IW]; if(bt==NULL_BOUNDARY) continue;
            int II=m.BC_II?m.BC_II[IW]:0, JJ=m.BC_JJ?m.BC_JJ[IW]:0, KK=m.BC_KK?m.BC_KK[IW]:0;
            int IIG=m.BC_IIG[IW],JJG=m.BC_JJG[IW],KKG=m.BC_KKG[IW],IOR=m.BC_IOR[IW];
            int IC=-1;
            if (m.BC_II && m.BC_JJ && m.BC_KK) IC=m.ci(II,JJ,KK);

            if(bt==SOLID_BOUNDARY && IC>=0 && m.CELL_SOLID && m.CELL_EXTERIOR &&
               m.CELL_SOLID[IC]==0.0 && m.CELL_EXTERIOR[IC]==0.0){
                // 内部固壁
                switch(IOR){
                    case 1: m.FX[m.i4_fx(IIG-1,JJG,KKG,0)]=0.0;break;
                    case-1: m.FX[m.i4_fx(IIG,JJG,KKG,0)]=0.0;break;
                    case 2: m.FY[m.i4_fx(IIG,JJG-1,KKG,0)]=0.0;break;
                    case-2: m.FY[m.i4_fx(IIG,JJG,KKG,0)]=0.0;break;
                    case 3: m.FZ[m.i4_fx(IIG,JJG,KKG-1,0)]=0.0;break;
                    case-3: m.FZ[m.i4_fx(IIG,JJG,KKG,0)]=0.0;break;
                }
            } else {
                // 其他边界
                double fv=0.0;
                if(m.BP1_RHO_F && m.BP1_ZZ_F){
                    double zz_get[512],mw_f;
                    for(int ns=0;ns<nts;++ns) zz_get[ns]=m.BP1_ZZ_F[IW*nts+ns];
                    fds_get_molecular_weight(zz_get,&mw_f);
                    fv=m.BP1_RHO_F[IW]/mw_f;
                }
                switch(IOR){
                    case 1: m.FX[m.i4_fx(IIG-1,JJG,KKG,0)]=fv;break;
                    case-1: m.FX[m.i4_fx(IIG,JJG,KKG,0)]=fv;break;
                    case 2: m.FY[m.i4_fx(IIG,JJG-1,KKG,0)]=fv;break;
                    case-2: m.FY[m.i4_fx(IIG,JJG,KKG,0)]=fv;break;
                    case 3: m.FZ[m.i4_fx(IIG,JJG,KKG-1,0)]=fv;break;
                    case-3: m.FZ[m.i4_fx(IIG,JJG,KKG,0)]=fv;break;
                }
            }

            // OFF_WALL_IF_3

            if(bt!=INTERPOLATED_BOUNDARY && bt!=OPEN_BOUNDARY &&
               m.BC_II && m.BC_JJ && m.BC_KK && m.CELL_WALL_INDEX){
                double zt[64]={},ut[64]={},ft[64]={};
                switch(IOR){
                    case 1:
                        if(UU[m.i3p1(II+1,JJ,KK)]>0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II+1,JJ,KK)*7)+4]<=0){
                            zt[20]=RHO_RMW[m.i3p1(II+1,JJ,KK)];
                            zt[21]=RHO_RMW[m.i3p1(II+1,JJ,KK)];
                            zt[22]=RHO_RMW[m.i3p1(II+2,JJ,KK)];
                            zt[23]=0.0;
                            ut[21]=UU[m.i3p1(II+1,JJ,KK)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,1,m.i_flux_limiter,4,4,4);
                            m.FX[m.i4_fx(II+1,JJ,KK,0)]=ft[21];
                        }
                        break;
                    case -1:
                        if(UU[m.i3p1(II-2,JJ,KK)]<0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II-1,JJ,KK)*7)+2]<=0){
                            zt[20]=0.0;
                            zt[21]=RHO_RMW[m.i3p1(II-2,JJ,KK)];
                            zt[22]=RHO_RMW[m.i3p1(II-1,JJ,KK)];
                            zt[23]=RHO_RMW[m.i3p1(II-1,JJ,KK)];
                            ut[21]=UU[m.i3p1(II-2,JJ,KK)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,1,m.i_flux_limiter,4,4,4);
                            m.FX[m.i4_fx(II-2,JJ,KK,0)]=ft[21];
                        }
                        break;
                    case 2:
                        if(VV[m.i3p1(II,JJ+1,KK)]>0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II,JJ+1,KK)*7)+5]<=0){
                            zt[17]=RHO_RMW[m.i3p1(II,JJ+1,KK)];
                            zt[21]=RHO_RMW[m.i3p1(II,JJ+1,KK)];
                            zt[25]=RHO_RMW[m.i3p1(II,JJ+2,KK)];
                            zt[29]=0.0;
                            ut[21]=VV[m.i3p1(II,JJ+1,KK)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,2,m.i_flux_limiter,4,4,4);
                            m.FY[m.i4_fx(II,JJ+1,KK,0)]=ft[21];
                        }
                        break;
                    case -2:
                        if(VV[m.i3p1(II,JJ-2,KK)]<0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II,JJ-1,KK)*7)+1]<=0){
                            zt[17]=0.0;
                            zt[21]=RHO_RMW[m.i3p1(II,JJ-2,KK)];
                            zt[25]=RHO_RMW[m.i3p1(II,JJ-1,KK)];
                            zt[29]=RHO_RMW[m.i3p1(II,JJ-1,KK)];
                            ut[21]=VV[m.i3p1(II,JJ-2,KK)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,2,m.i_flux_limiter,4,4,4);
                            m.FY[m.i4_fx(II,JJ-2,KK,0)]=ft[21];
                        }
                        break;
                    case 3:
                        if(WW[m.i3p1(II,JJ,KK+1)]>0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II,JJ,KK+1)*7)+6]<=0){
                            zt[5]=RHO_RMW[m.i3p1(II,JJ,KK+1)];
                            zt[21]=RHO_RMW[m.i3p1(II,JJ,KK+1)];
                            zt[37]=RHO_RMW[m.i3p1(II,JJ,KK+2)];
                            zt[53]=0.0;
                            ut[21]=WW[m.i3p1(II,JJ,KK+1)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,3,m.i_flux_limiter,4,4,4);
                            m.FZ[m.i4_fx(II,JJ,KK+1,0)]=ft[21];
                        }
                        break;
                    case -3:
                        if(WW[m.i3p1(II,JJ,KK-2)]<0.0 &&
                           m.CELL_WALL_INDEX[(m.ci(II,JJ,KK-1)*7)+0]<=0){
                            zt[5]=0.0;
                            zt[21]=RHO_RMW[m.i3p1(II,JJ,KK-2)];
                            zt[37]=RHO_RMW[m.i3p1(II,JJ,KK-1)];
                            zt[53]=RHO_RMW[m.i3p1(II,JJ,KK-1)];
                            ut[21]=WW[m.i3p1(II,JJ,KK-2)];
                            fds_get_scalar_face_value(ut,zt,ft,1,1,1,1,1,1,3,m.i_flux_limiter,4,4,4);
                            m.FZ[m.i4_fx(II,JJ,KK-2,0)]=ft[21];
                        }
                        break;
                }
            }
        }
        // 最大面值修正
        #pragma omp parallel for private(zz_get,mw_g)
        for(int K=0;K<=kbar;++K) for(int J=0;J<=jbar;++J) for(int I=0;I<=ibar;++I){
            for(int dim=0;dim<3;++dim){
                double* F=(dim==0?m.FX:(dim==1?m.FY:m.FZ));
                int nmax=1;double fmax=F[m.i4_fx(I,J,K,1)];
                for(int ns=2;ns<=nts;++ns){double fv=F[m.i4_fx(I,J,K,ns)];if(fv>fmax){fmax=fv;nmax=ns;}}
                mw_g=m.SM_MW[nmax-1];double sm=0;
                for(int ns=1;ns<=nts;++ns) if(ns!=nmax) sm+=F[m.i4_fx(I,J,K,ns)]/m.SM_MW[ns-1];
                F[m.i4_fx(I,J,K,nmax)]=mw_g*std::max(0.0,F[m.i4_fx(I,J,K,0)]-sm);
            }
        }
    }
}

// =============================================================================
// DENSITY (行 341-952)
// =============================================================================
void density_impl(double T,double DT,const MassMeshPointersC* c){
    {
        static FILE* fp=std::fopen("/tmp/mass_c_density_called.txt","w");
        if(fp){std::fprintf(fp,"density_impl called T=%.6e DT=%.6e\n",T,DT);std::fflush(fp);}
    }
    MassB m(c); if(m.solid_phase_only) return;
    if(m.period_test==0){if(m.icyc<=1) return;}
    else if(m.period_test==5||m.period_test==8) return;
    int ibar=m.ibar,jbar=m.jbar,kbar=m.kbar;
    int ibp1=m.ibp1,jbp1=m.jbp1,kbp1=m.kbp1;
    int nts=m.nts,n_total=m.n_total,n_zone=m.n_zone;
    double* UU_W=m.WORK1; double* VV_W=m.WORK2; double* WW_W=m.WORK3;

    if(m.predictor){
        // FIRST_PASS
        if(m.first_pass){
            if(m.depo()&&(m.grav_settling||m.ther_settling)){
                // SETTLING_VELOCITY 在 Fortran 侧调用
            }
            // DEL_RHO_D_DEL_Z 分配为 (0:IBP1,0:JBP1,0:KBP1,N_TOTAL_SCALARS)
            std::memcpy(m.SWORK4,m.DEL_RHO_D_DEL_Z,(ibp1+1)*(jbp1+1)*(kbp1+1)*n_total*sizeof(double));
        }
        // UU=U,VV=V,WW=W (行 393-395) — U/V/W 分配为 (0:IBP1,0:JBP1,0:KBP1), 步长 IBP1+1
        std::memcpy(UU_W,m.U,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));
        std::memcpy(VV_W,m.V,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));
        std::memcpy(WW_W,m.W,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));

        // WALL_LOOP 
        #pragma omp parallel
        {
            #pragma omp for
            for(int IW=0;IW<m.nw_ext;++IW){
                if(m.WALL_BC_TYPE[IW]!=INTERPOLATED_BOUNDARY) continue;
                int IIG=m.BC_IIG[IW],JJG=m.BC_JJG[IW],KKG=m.BC_KKG[IW],IOR=m.BC_IOR[IW];
                double uvs=m.UVW_SAVE[IW];
                switch(IOR){
                    case 1: UU_W[m.i3p1(IIG-1,JJG,KKG)]=uvs;break;
                    case-1: UU_W[m.i3p1(IIG,JJG,KKG)]=uvs;break;
                    case 2: VV_W[m.i3p1(IIG,JJG-1,KKG)]=uvs;break;
                    case-2: VV_W[m.i3p1(IIG,JJG,KKG)]=uvs;break;
                    case 3: WW_W[m.i3p1(IIG,JJG,KKG-1)]=uvs;break;
                    case-3: WW_W[m.i3p1(IIG,JJG,KKG)]=uvs;break;
                }
            }
            // RHS+ZZS
            #pragma omp for
            for(int N=1;N<=n_total;++N) for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
                if(m.iss(I,J,K)) continue;
                double RHS=-m.SWORK4[m.i4_drdz(I,J,K,N)]
                    +((m.FX[m.i4_fx(I,J,K,N)]*UU_W[m.i3p1(I,J,K)]*m.R[I]-m.FX[m.i4_fx(I-1,J,K,N)]*UU_W[m.i3p1(I-1,J,K)]*m.R[I-1])*m.RDX[I]*m.RRN[I])
                    +((m.FY[m.i4_fx(I,J,K,N)]*VV_W[m.i3p1(I,J,K)]-m.FY[m.i4_fx(I,J-1,K,N)]*VV_W[m.i3p1(I,J-1,K)])*m.RDY[J])
                    +((m.FZ[m.i4_fx(I,J,K,N)]*WW_W[m.i3p1(I,J,K)]-m.FZ[m.i4_fx(I,J,K-1,N)]*WW_W[m.i3p1(I,J,K-1)])*m.RDZ[K]);
                m.ZZS[m.i4_zz(I,J,K,N)]=m.RHO[m.i3p1(I,J,K)]*m.ZZ[m.i4_zz(I,J,K,N)]-DT*RHS;
            }
        }
        // STORE_SPECIES_FLUX 
        if(m.store_species_flux) for(int N=1;N<=n_total;++N) for(int K=0;K<=kbp1;++K)
            for(int J=0;J<=jbp1;++J) for(int I=0;I<=ibp1;++I){
                m.ADV_FX[m.i4_adv(I,J,K,N)]=m.FX[m.i4_fx(I,J,K,N)]*UU_W[m.i3p1(I,J,K)];
                m.ADV_FY[m.i4_adv(I,J,K,N)]=m.FY[m.i4_fx(I,J,K,N)]*VV_W[m.i3p1(I,J,K)];
                m.ADV_FZ[m.i4_adv(I,J,K,N)]=m.FZ[m.i4_fx(I,J,K,N)]*WW_W[m.i3p1(I,J,K)];
            }
        // 源项
        if(m.have_m_dot_ppp) for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I)
            for(int N=1;N<=nts;++N) m.ZZS[m.i4_zz(I,J,K,N)]+=DT*m.M_DOT_PPP[m.i4_mdot(I,J,K,N)];

        // RHOS=SUM(ZZS) (行 472-481)
        #pragma omp parallel for
        for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
            if(m.iss(I,J,K)) continue;
            double s=0; for(int N=1;N<=nts;++N) s+=m.ZZS[m.i4_zz(I,J,K,N)];
            m.RHOS[m.i3p1(I,J,K)]=s;
        }

        // CHECK_MASS_DENSITY 
        {
            bool clip_rhomin=false,clip_rhomax=false;
            double* DELTA_RHO=m.WORK4;
            std::memset(DELTA_RHO,0,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));
            double* RHO_ZZ_ptr=m.ZZS; double* RHOP_ptr=m.RHOS;
            double RHO_MIN=m.rhomin,RHO_MAX=m.rhomax;
            // 密度裁剪
            for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J){
                double vc1[7]={m.DY[J]*m.DZ[K],0,0,m.DY[J-1]*m.DZ[K],m.DY[J+1]*m.DZ[K],m.DY[J]*m.DZ[K-1],m.DY[J]*m.DZ[K+1]};
                vc1[1]=vc1[0];vc1[2]=vc1[0];
                for(int I=1;I<=ibar;++I){
                    double rh=RHOP_ptr[m.i3p1(I,J,K)];
                    if(rh>=RHO_MIN&&rh<=RHO_MAX) continue;
                    int IC=m.ci(I,J,K); if(IC>=0&&m.CELL_SOLID[IC]!=0.0) continue;
                    double RHO_CUT,SIGN_FACTOR;
                    if(rh<RHO_MIN){RHO_CUT=RHO_MIN;SIGN_FACTOR=1.0;clip_rhomin=true;}
                    else{RHO_CUT=RHO_MAX;SIGN_FACTOR=-1.0;clip_rhomax=true;}
                    double vc[7]={m.DX[I]*vc1[0],m.DX[I-1]*vc1[1],m.DX[I+1]*vc1[2],m.DX[I]*vc1[3],m.DX[I]*vc1[4],m.DX[I]*vc1[5],m.DX[I]*vc1[6]};
                    double MASS_C=std::abs(RHO_CUT-rh)*vc[0];
                    double mn[7]={0};
                    if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+2]==0) mn[1]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I-1,J,K)]))-RHO_CUT)*vc[1];
                    if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+4]==0) mn[2]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I+1,J,K)]))-RHO_CUT)*vc[2];
                    if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+1]==0) mn[3]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I,J-1,K)]))-RHO_CUT)*vc[3];
                    if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+5]==0) mn[4]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I,J+1,K)]))-RHO_CUT)*vc[4];
                    if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+0]==0) mn[5]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I,J,K-1)]))-RHO_CUT)*vc[5];
                    if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+6]==0) mn[6]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I,J,K+1)]))-RHO_CUT)*vc[6];
                    double sum_mn=mn[1]+mn[2]+mn[3]+mn[4]+mn[5]+mn[6];
                    if(sum_mn<=2.0*2.22e-16) continue;
                    double CONST=SIGN_FACTOR*std::min(1.0,MASS_C/sum_mn);
                    int off[7]={0};
                    DELTA_RHO[m.i3p1(I,J,K)]  +=CONST*sum_mn/vc[0];
                    if(mn[1]>0) DELTA_RHO[m.i3p1(I-1,J,K)]-=CONST*mn[1]/vc[1];
                    if(mn[2]>0) DELTA_RHO[m.i3p1(I+1,J,K)]-=CONST*mn[2]/vc[2];
                    if(mn[3]>0) DELTA_RHO[m.i3p1(I,J-1,K)]-=CONST*mn[3]/vc[3];
                    if(mn[4]>0) DELTA_RHO[m.i3p1(I,J+1,K)]-=CONST*mn[4]/vc[4];
                    if(mn[5]>0) DELTA_RHO[m.i3p1(I,J,K-1)]-=CONST*mn[5]/vc[5];
                    if(mn[6]>0) DELTA_RHO[m.i3p1(I,J,K+1)]-=CONST*mn[6]/vc[6];
                }
            }
            if(clip_rhomin||clip_rhomax){
                for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
                    double val=RHOP_ptr[m.i3p1(I,J,K)]+DELTA_RHO[m.i3p1(I,J,K)];
                    RHOP_ptr[m.i3p1(I,J,K)]=std::min(RHO_MAX,std::max(RHO_MIN,val));
                }
            }
            if(nts==1){
                if(clip_rhomin||clip_rhomax)
                    for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I)
                        RHO_ZZ_ptr[m.i4_zz(I,J,K,1)]=RHOP_ptr[m.i3p1(I,J,K)];
                goto end_check; // N_TRACKED_SPECIES==1 时直接跳到结束
            }
          
            bool clip_rho_zz=false;
            double RHO_ZZ_MIN=0.0;
            for(int N=1;N<=nts;++N){
                double* DELTA_RHO_ZZ=m.WORK5;
                std::memset(DELTA_RHO_ZZ,0,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));
                bool clip_rz=false;
                for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J){
                    double vc1[7]={m.DY[J]*m.DZ[K],0,0,m.DY[J-1]*m.DZ[K],m.DY[J+1]*m.DZ[K],m.DY[J]*m.DZ[K-1],m.DY[J]*m.DZ[K+1]};
                    vc1[1]=vc1[0];vc1[2]=vc1[0];
                    for(int I=1;I<=ibar;++I){
                        int IC=m.ci(I,J,K); if(IC>=0&&m.CELL_SOLID[IC]!=0.0) continue;
                        double rhomax=RHOP_ptr[m.i3p1(I,J,K)];
                        double rzv=RHO_ZZ_ptr[m.i4_zz(I,J,K,N)];
                        if(rzv>=RHO_ZZ_MIN&&rzv<=rhomax) continue;
                        clip_rz=true; clip_rho_zz=true;
                        double RHO_ZZ_CUT,SIGN_FACTOR2;
                        if(rzv<RHO_ZZ_MIN){RHO_ZZ_CUT=RHO_ZZ_MIN;SIGN_FACTOR2=1.0;}
                        else{RHO_ZZ_CUT=rhomax;SIGN_FACTOR2=-1.0;}
                        double vc[7]={m.DX[I]*vc1[0],m.DX[I-1]*vc1[1],m.DX[I+1]*vc1[2],m.DX[I]*vc1[3],m.DX[I]*vc1[4],m.DX[I]*vc1[5],m.DX[I]*vc1[6]};
                        double MASS_C=std::abs(RHO_ZZ_CUT-rzv)*vc[0];
                        double mn[7]={0};
                        if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+2]==0) mn[1]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I-1,J,K,N)]))-RHO_ZZ_CUT)*vc[1];
                        if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+4]==0) mn[2]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I+1,J,K,N)]))-RHO_ZZ_CUT)*vc[2];
                        if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+1]==0) mn[3]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I,J-1,K,N)]))-RHO_ZZ_CUT)*vc[3];
                        if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+5]==0) mn[4]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I,J+1,K,N)]))-RHO_ZZ_CUT)*vc[4];
                        if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+0]==0) mn[5]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I,J,K-1,N)]))-RHO_ZZ_CUT)*vc[5];
                        if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+6]==0) mn[6]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I,J,K+1,N)]))-RHO_ZZ_CUT)*vc[6];
                        double sum_mn=mn[1]+mn[2]+mn[3]+mn[4]+mn[5]+mn[6];
                        if(sum_mn<=2.0*2.22e-16) continue;
                        double CONST=SIGN_FACTOR2*std::min(1.0,MASS_C/sum_mn);
                        DELTA_RHO_ZZ[m.i3p1(I,J,K)]  +=CONST*sum_mn/vc[0];
                        if(mn[1]>0) DELTA_RHO_ZZ[m.i3p1(I-1,J,K)]-=CONST*mn[1]/vc[1];
                        if(mn[2]>0) DELTA_RHO_ZZ[m.i3p1(I+1,J,K)]-=CONST*mn[2]/vc[2];
                        if(mn[3]>0) DELTA_RHO_ZZ[m.i3p1(I,J-1,K)]-=CONST*mn[3]/vc[3];
                        if(mn[4]>0) DELTA_RHO_ZZ[m.i3p1(I,J+1,K)]-=CONST*mn[4]/vc[4];
                        if(mn[5]>0) DELTA_RHO_ZZ[m.i3p1(I,J,K-1)]-=CONST*mn[5]/vc[5];
                        if(mn[6]>0) DELTA_RHO_ZZ[m.i3p1(I,J,K+1)]-=CONST*mn[6]/vc[6];
                    }
                }
                if(clip_rz){
                    for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
                        double val=RHO_ZZ_ptr[m.i4_zz(I,J,K,N)]+DELTA_RHO_ZZ[m.i3p1(I,J,K)];
                        RHO_ZZ_ptr[m.i4_zz(I,J,K,N)]=std::min(RHOP_ptr[m.i3p1(I,J,K)],std::max(RHO_ZZ_MIN,val));
                    }
                }
            }
            // 最终检查和归一化
            if(clip_rhomin||clip_rhomax||clip_rho_zz){
                for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
                    if(m.iss(I,J,K)) continue;
                    double sum_rz=0; int nmax=1; double fmax=RHO_ZZ_ptr[m.i4_zz(I,J,K,1)];
                    for(int N=1;N<=nts;++N){
                        double v=RHO_ZZ_ptr[m.i4_zz(I,J,K,N)];
                        sum_rz+=v; if(v>fmax){fmax=v;nmax=N;}
                    }
                    double test=RHO_ZZ_ptr[m.i4_zz(I,J,K,nmax)]+RHOP_ptr[m.i3p1(I,J,K)]-sum_rz;
                    if(test<0.0||test>RHOP_ptr[m.i3p1(I,J,K)]){
                        for(int N=1;N<=nts;++N)
                            RHO_ZZ_ptr[m.i4_zz(I,J,K,N)]=RHOP_ptr[m.i3p1(I,J,K)]*RHO_ZZ_ptr[m.i4_zz(I,J,K,N)]/sum_rz;
                    } else {
                        RHO_ZZ_ptr[m.i4_zz(I,J,K,nmax)]=test;
                    }
                }
            }
       }
       end_check:;

       // ZZS = ZZS / RHOS
       #pragma omp parallel for
       for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
           if(m.iss(I,J,K)) continue;
           for(int N=1;N<=n_total;++N)
               m.ZZS[m.i4_zz(I,J,K,N)]=m.ZZS[m.i4_zz(I,J,K,N)]/m.RHOS[m.i3p1(I,J,K)];
       }

       // CLIP_PASSIVE_SCALARS
       if(m.n_passive>0){
           int zi=m.zeta_idx;
           for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
               if(m.iss(I,J,K)) continue;
               double val=m.ZZS[m.i4_zz(I,J,K,zi)];
               m.ZZS[m.i4_zz(I,J,K,zi)]=std::max(0.0,std::min(1.0,val));
           }
       }

   
       for(int I=1;I<=n_zone;++I){
           for(int K=0;K<=kbp1;++K){
               m.PBAR_S[K+(kbp1+1)*I]=m.PBAR[K+(kbp1+1)*I]+m.D_PBAR_DT[I-1]*DT;
           }
       }

       // RSUM 计算 (行 523-533)
       #pragma omp parallel for
       for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
           if(m.iss(I,J,K)) continue;
           double zz_get_loc[512];
           for(int N=1;N<=nts;++N) zz_get_loc[N-1]=m.ZZS[m.i4_zz(I,J,K,N)];
           fds_get_specific_gas_constant(zz_get_loc,&m.RSUM[m.i3p1(I,J,K)]);
       }

       // TMP = PBAR_S/(RSUM*RHOS) (行 537-546)
       #pragma omp parallel for
       for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
           if(m.iss(I,J,K)) continue;
           int PZ=m.PRESSURE_ZONE[m.i3p1(I,J,K)];
           double pbs=pbar_at(m.PBAR_S,kbp1+1,K,PZ);
           m.TMP[m.i3p1(I,J,K)]=pbs/(m.RSUM[m.i3p1(I,J,K)]*m.RHOS[m.i3p1(I,J,K)]);
       }

       {
           static int dbg_cnt=0;
           if(dbg_cnt<20){
               int II=ibar/2, JJ=jbar/2, KK=kbar/2;
               int idx=m.i3p1(II,JJ,KK);
               int PZ=m.PRESSURE_ZONE[idx];
               FILE* fp=std::fopen("/tmp/mass_c_tmp_dbg.txt","a");
               if(fp){
                   std::fprintf(fp,"PRED T=%.4e PZ=%d TMP=%.6e RHOS=%.6e RSUM=%.6e PBAR_S=%.6e\n",
                       T,PZ,m.TMP[idx],m.RHOS[idx],m.RSUM[idx],pbar_at(m.PBAR_S,kbp1+1,KK,PZ));
                   std::fclose(fp);
               }
               ++dbg_cnt;
           }
       }

   // ====================================================================
   // CORRECTOR 步
   // ====================================================================
   } else {
       // UU=US, VV=VS, WW=WS 
       std::memcpy(UU_W,m.US,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));
       std::memcpy(VV_W,m.VS,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));
       std::memcpy(WW_W,m.WS,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));

       // WALL_LOOP_2 
       #pragma omp parallel for
       for(int IW=0;IW<m.nw_ext;++IW){
           int bt_prev=m.EWC_BOUNDARY_TYPE_PREVIOUS?m.EWC_BOUNDARY_TYPE_PREVIOUS[IW]:0;
           if(bt_prev!=INTERPOLATED_BOUNDARY) continue;
           int IIG=m.BC_IIG[IW],JJG=m.BC_JJG[IW],KKG=m.BC_KKG[IW],IOR=m.BC_IOR[IW];
           double uvs=m.UVW_SAVE[IW];
           switch(IOR){
               case 1: UU_W[m.i3p1(IIG-1,JJG,KKG)]=uvs;break;
               case-1: UU_W[m.i3p1(IIG,JJG,KKG)]=uvs;break;
               case 2: VV_W[m.i3p1(IIG,JJG-1,KKG)]=uvs;break;
               case-2: VV_W[m.i3p1(IIG,JJG,KKG)]=uvs;break;
               case 3: WW_W[m.i3p1(IIG,JJG,KKG-1)]=uvs;break;
               case-3: WW_W[m.i3p1(IIG,JJG,KKG)]=uvs;break;
           }
       }
       // SETTLING_VELOCITY 
       // RHS+ZZ
       #pragma omp parallel for
       for(int N=1;N<=n_total;++N) for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
           if(m.iss(I,J,K)) continue;
           double RHS=-m.DEL_RHO_D_DEL_Z[m.i4_drdz(I,J,K,N)]
               +((m.FX[m.i4_fx(I,J,K,N)]*UU_W[m.i3p1(I,J,K)]*m.R[I]-m.FX[m.i4_fx(I-1,J,K,N)]*UU_W[m.i3p1(I-1,J,K)]*m.R[I-1])*m.RDX[I]*m.RRN[I])
               +((m.FY[m.i4_fx(I,J,K,N)]*VV_W[m.i3p1(I,J,K)]-m.FY[m.i4_fx(I,J-1,K,N)]*VV_W[m.i3p1(I,J-1,K)])*m.RDY[J])
               +((m.FZ[m.i4_fx(I,J,K,N)]*WW_W[m.i3p1(I,J,K)]-m.FZ[m.i4_fx(I,J,K-1,N)]*WW_W[m.i3p1(I,J,K-1)])*m.RDZ[K]);
           m.ZZ[m.i4_zz(I,J,K,N)]=0.5*(m.RHO[m.i3p1(I,J,K)]*m.ZZ[m.i4_zz(I,J,K,N)]+m.RHOS[m.i3p1(I,J,K)]*m.ZZS[m.i4_zz(I,J,K,N)]-DT*RHS);
       }
       // 源项
       if(m.have_m_dot_ppp){
           for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I)
               for(int N=1;N<=nts;++N) m.ZZ[m.i4_zz(I,J,K,N)]+=0.5*DT*m.M_DOT_PPP[m.i4_mdot(I,J,K,N)];
           if(!m.cc_ibm){
               // M_DOT_PPP 分配为 (0:IBP1,0:JBP1,0:KBP1,1:N_TRACKED_SPECIES)
               std::memset(m.M_DOT_PPP,0,(ibp1+1)*(jbp1+1)*(kbp1+1)*nts*sizeof(double));
               std::memset(m.D_SOURCE,0,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));
           }
       }
       // STORE_SPECIES_FLUX (行 616-622) — Fortran ADV_FX(:,:,:,N)=0.5*(ADV_FX(:,:,:,N)+FX(:,:,:,N)*UU(:,:,:)) 全数组
       if(m.store_species_flux) for(int N=1;N<=n_total;++N) for(int K=0;K<=kbp1;++K)
           for(int J=0;J<=jbp1;++J) for(int I=0;I<=ibp1;++I){
               m.ADV_FX[m.i4_adv(I,J,K,N)]=0.5*(m.ADV_FX[m.i4_adv(I,J,K,N)]+m.FX[m.i4_fx(I,J,K,N)]*UU_W[m.i3p1(I,J,K)]);
               m.ADV_FY[m.i4_adv(I,J,K,N)]=0.5*(m.ADV_FY[m.i4_adv(I,J,K,N)]+m.FY[m.i4_fx(I,J,K,N)]*VV_W[m.i3p1(I,J,K)]);
               m.ADV_FZ[m.i4_adv(I,J,K,N)]=0.5*(m.ADV_FZ[m.i4_adv(I,J,K,N)]+m.FZ[m.i4_fx(I,J,K,N)]*WW_W[m.i3p1(I,J,K)]);
           }

       // RHO = SUM(ZZ) 
       #pragma omp parallel for
       for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
           if(m.iss(I,J,K)) continue;
           double s=0; for(int N=1;N<=nts;++N) s+=m.ZZ[m.i4_zz(I,J,K,N)];
           m.RHO[m.i3p1(I,J,K)]=s;
       }

       // CHECK_MASS_DENSITY — 用 ZZ(此时为 RHO*ZZ)/RHO 
           // 重复 CHECK_MASS_DENSITY 逻辑, 但使用 m.ZZ 和 m.RHO
           bool clip_rhomin=false,clip_rhomax=false;
           double* DELTA_RHO=m.WORK4;
           std::memset(DELTA_RHO,0,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));
           double* RHO_ZZ_ptr=m.ZZ; double* RHOP_ptr=m.RHO;
           double RHO_MIN=m.rhomin,RHO_MAX=m.rhomax;
          
           for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J){
               double vc1[7]={m.DY[J]*m.DZ[K],0,0,m.DY[J-1]*m.DZ[K],m.DY[J+1]*m.DZ[K],m.DY[J]*m.DZ[K-1],m.DY[J]*m.DZ[K+1]};
               vc1[1]=vc1[0];vc1[2]=vc1[0];
               for(int I=1;I<=ibar;++I){
                   double rh=RHOP_ptr[m.i3p1(I,J,K)];
                   if(rh>=RHO_MIN&&rh<=RHO_MAX) continue;
                   int IC=m.ci(I,J,K); if(IC>=0&&m.CELL_SOLID[IC]!=0.0) continue;
                   double RHO_CUT,SIGN_FACTOR;
                   if(rh<RHO_MIN){RHO_CUT=RHO_MIN;SIGN_FACTOR=1.0;clip_rhomin=true;}
                   else{RHO_CUT=RHO_MAX;SIGN_FACTOR=-1.0;clip_rhomax=true;}
                   double vc[7]={m.DX[I]*vc1[0],m.DX[I-1]*vc1[1],m.DX[I+1]*vc1[2],m.DX[I]*vc1[3],m.DX[I]*vc1[4],m.DX[I]*vc1[5],m.DX[I]*vc1[6]};
                   double MASS_C=std::abs(RHO_CUT-rh)*vc[0];
                   double mn[7]={0};
                   if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+2]==0) mn[1]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I-1,J,K)]))-RHO_CUT)*vc[1];
                   if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+4]==0) mn[2]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I+1,J,K)]))-RHO_CUT)*vc[2];
                   if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+1]==0) mn[3]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I,J-1,K)]))-RHO_CUT)*vc[3];
                   if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+5]==0) mn[4]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I,J+1,K)]))-RHO_CUT)*vc[4];
                   if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+0]==0) mn[5]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I,J,K-1)]))-RHO_CUT)*vc[5];
                   if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+6]==0) mn[6]=std::abs(std::min(RHO_MAX,std::max(RHO_MIN,RHOP_ptr[m.i3p1(I,J,K+1)]))-RHO_CUT)*vc[6];
                   double sum_mn=mn[1]+mn[2]+mn[3]+mn[4]+mn[5]+mn[6];
                   if(sum_mn<=2.0*2.22e-16) continue;
                   double CONST=SIGN_FACTOR*std::min(1.0,MASS_C/sum_mn);
                   DELTA_RHO[m.i3p1(I,J,K)]  +=CONST*sum_mn/vc[0];
                   if(mn[1]>0) DELTA_RHO[m.i3p1(I-1,J,K)]-=CONST*mn[1]/vc[1];
                   if(mn[2]>0) DELTA_RHO[m.i3p1(I+1,J,K)]-=CONST*mn[2]/vc[2];
                   if(mn[3]>0) DELTA_RHO[m.i3p1(I,J-1,K)]-=CONST*mn[3]/vc[3];
                   if(mn[4]>0) DELTA_RHO[m.i3p1(I,J+1,K)]-=CONST*mn[4]/vc[4];
                   if(mn[5]>0) DELTA_RHO[m.i3p1(I,J,K-1)]-=CONST*mn[5]/vc[5];
                   if(mn[6]>0) DELTA_RHO[m.i3p1(I,J,K+1)]-=CONST*mn[6]/vc[6];
               }
           }
           if(clip_rhomin||clip_rhomax){
               for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
                   double val=RHOP_ptr[m.i3p1(I,J,K)]+DELTA_RHO[m.i3p1(I,J,K)];
                   RHOP_ptr[m.i3p1(I,J,K)]=std::min(RHO_MAX,std::max(RHO_MIN,val));
               }
           }
           if(nts==1){
               if(clip_rhomin||clip_rhomax)
                   for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I)
                       RHO_ZZ_ptr[m.i4_zz(I,J,K,1)]=RHOP_ptr[m.i3p1(I,J,K)];
               goto end_check2;
           }
           
           bool clip_rho_zz=false;
           double RHO_ZZ_MIN=0.0;
           for(int N=1;N<=nts;++N){
               double* DELTA_RHO_ZZ=m.WORK5;
               std::memset(DELTA_RHO_ZZ,0,(ibp1+1)*(jbp1+1)*(kbp1+1)*sizeof(double));
               bool clip_rz=false;
               for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J){
                   double vc1[7]={m.DY[J]*m.DZ[K],0,0,m.DY[J-1]*m.DZ[K],m.DY[J+1]*m.DZ[K],m.DY[J]*m.DZ[K-1],m.DY[J]*m.DZ[K+1]};
                   vc1[1]=vc1[0];vc1[2]=vc1[0];
                   for(int I=1;I<=ibar;++I){
                       int IC=m.ci(I,J,K); if(IC>=0&&m.CELL_SOLID[IC]!=0.0) continue;
                       double rhomax=RHOP_ptr[m.i3p1(I,J,K)];
                       double rzv=RHO_ZZ_ptr[m.i4_zz(I,J,K,N)];
                       if(rzv>=RHO_ZZ_MIN&&rzv<=rhomax) continue;
                       clip_rz=true; clip_rho_zz=true;
                       double RHO_ZZ_CUT,SIGN_FACTOR2;
                       if(rzv<RHO_ZZ_MIN){RHO_ZZ_CUT=RHO_ZZ_MIN;SIGN_FACTOR2=1.0;}
                       else{RHO_ZZ_CUT=rhomax;SIGN_FACTOR2=-1.0;}
                       double vc[7]={m.DX[I]*vc1[0],m.DX[I-1]*vc1[1],m.DX[I+1]*vc1[2],m.DX[I]*vc1[3],m.DX[I]*vc1[4],m.DX[I]*vc1[5],m.DX[I]*vc1[6]};
                       double MASS_C=std::abs(RHO_ZZ_CUT-rzv)*vc[0];
                       double mn[7]={0};
                       if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+2]==0) mn[1]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I-1,J,K,N)]))-RHO_ZZ_CUT)*vc[1];
                       if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+4]==0) mn[2]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I+1,J,K,N)]))-RHO_ZZ_CUT)*vc[2];
                       if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+1]==0) mn[3]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I,J-1,K,N)]))-RHO_ZZ_CUT)*vc[3];
                       if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+5]==0) mn[4]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I,J+1,K,N)]))-RHO_ZZ_CUT)*vc[4];
                       if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+0]==0) mn[5]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I,J,K-1,N)]))-RHO_ZZ_CUT)*vc[5];
                       if(!m.CELL_WALL_INDEX||m.CELL_WALL_INDEX[IC*7+6]==0) mn[6]=std::abs(std::min(rhomax,std::max(RHO_ZZ_MIN,RHO_ZZ_ptr[m.i4_zz(I,J,K+1,N)]))-RHO_ZZ_CUT)*vc[6];
                       double sum_mn=mn[1]+mn[2]+mn[3]+mn[4]+mn[5]+mn[6];
                       if(sum_mn<=2.0*2.22e-16) continue;
                       double CONST=SIGN_FACTOR2*std::min(1.0,MASS_C/sum_mn);
                       DELTA_RHO_ZZ[m.i3p1(I,J,K)]  +=CONST*sum_mn/vc[0];
                       if(mn[1]>0) DELTA_RHO_ZZ[m.i3p1(I-1,J,K)]-=CONST*mn[1]/vc[1];
                       if(mn[2]>0) DELTA_RHO_ZZ[m.i3p1(I+1,J,K)]-=CONST*mn[2]/vc[2];
                       if(mn[3]>0) DELTA_RHO_ZZ[m.i3p1(I,J-1,K)]-=CONST*mn[3]/vc[3];
                       if(mn[4]>0) DELTA_RHO_ZZ[m.i3p1(I,J+1,K)]-=CONST*mn[4]/vc[4];
                       if(mn[5]>0) DELTA_RHO_ZZ[m.i3p1(I,J,K-1)]-=CONST*mn[5]/vc[5];
                       if(mn[6]>0) DELTA_RHO_ZZ[m.i3p1(I,J,K+1)]-=CONST*mn[6]/vc[6];
                   }
               }
               if(clip_rz){
                   for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
                       double val=RHO_ZZ_ptr[m.i4_zz(I,J,K,N)]+DELTA_RHO_ZZ[m.i3p1(I,J,K)];
                       RHO_ZZ_ptr[m.i4_zz(I,J,K,N)]=std::min(RHOP_ptr[m.i3p1(I,J,K)],std::max(RHO_ZZ_MIN,val));
                   }
               }
           }
           if(clip_rhomin||clip_rhomax||clip_rho_zz){
               for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
                   if(m.iss(I,J,K)) continue;
                   double sum_rz=0; int nmax=1; double fmax=RHO_ZZ_ptr[m.i4_zz(I,J,K,1)];
                   for(int N=1;N<=nts;++N){
                       double v=RHO_ZZ_ptr[m.i4_zz(I,J,K,N)];
                       sum_rz+=v; if(v>fmax){fmax=v;nmax=N;}
                   }
                   double test=RHO_ZZ_ptr[m.i4_zz(I,J,K,nmax)]+RHOP_ptr[m.i3p1(I,J,K)]-sum_rz;
                   if(test<0.0||test>RHOP_ptr[m.i3p1(I,J,K)]){
                       for(int N=1;N<=nts;++N)
                           RHO_ZZ_ptr[m.i4_zz(I,J,K,N)]=RHOP_ptr[m.i3p1(I,J,K)]*RHO_ZZ_ptr[m.i4_zz(I,J,K,N)]/sum_rz;
                   } else {
                       RHO_ZZ_ptr[m.i4_zz(I,J,K,nmax)]=test;
                   }
               }
           }
       }
       end_check2:;

       // ZZ = ZZ / RHO (行 669-678)
       #pragma omp parallel for
       for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
           if(m.iss(I,J,K)) continue;
           for(int N=1;N<=n_total;++N)
               m.ZZ[m.i4_zz(I,J,K,N)]=m.ZZ[m.i4_zz(I,J,K,N)]/m.RHO[m.i3p1(I,J,K)];
       }

       // CLIP_PASSIVE_SCALARS (行 682-685)
       if(m.n_passive>0){
           int zi=m.zeta_idx;
           for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
               if(m.iss(I,J,K)) continue;
               double val=m.ZZ[m.i4_zz(I,J,K,zi)];
               m.ZZ[m.i4_zz(I,J,K,zi)]=std::max(0.0,std::min(1.0,val));
           }
       }

       // PBAR 更新 
       // zone I (1..N_ZONE) 偏移 K+(kbp1+1)*I (zone 0 起点), D_PBAR_DT_S[I-1]
       // Fortran PBAR(:,I)=0.5*(PBAR(:,I)+PBAR_S(:,I)+D_PBAR_DT_S(I)*DT) 为整列 0:KBP1
       for(int I=1;I<=n_zone;++I){
           for(int K=0;K<=kbp1;++K){
               double pbar_k = m.PBAR[K+(kbp1+1)*I];
               double pbar_s_k = m.PBAR_S?m.PBAR_S[K+(kbp1+1)*I]:pbar_k;
               m.PBAR[K+(kbp1+1)*I] = 0.5*(pbar_k + pbar_s_k + m.D_PBAR_DT_S[I-1]*DT);
           }
       }

       // RSUM 计算 
       #pragma omp parallel for
       for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
           if(m.iss(I,J,K)) continue;
           double zz_get_loc[512];
           for(int N=1;N<=nts;++N) zz_get_loc[N-1]=m.ZZ[m.i4_zz(I,J,K,N)];
           fds_get_specific_gas_constant(zz_get_loc,&m.RSUM[m.i3p1(I,J,K)]);
       }

       // TMP = PBAR/(RSUM*RHO) (行 712-721)
       #pragma omp parallel for
       for(int K=1;K<=kbar;++K) for(int J=1;J<=jbar;++J) for(int I=1;I<=ibar;++I){
           if(m.iss(I,J,K)) continue;
           int PZ=m.PRESSURE_ZONE[m.i3p1(I,J,K)];
           double pb=pbar_at(m.PBAR,kbp1+1,K,PZ);
           m.TMP[m.i3p1(I,J,K)]=pb/(m.RSUM[m.i3p1(I,J,K)]*m.RHO[m.i3p1(I,J,K)]);
       }
   }
}

} // anonymous namespace

// =============================================================================
// extern "C" 桥接函数 — 从 Fortran 包装器调用
// =============================================================================
extern "C" void mass_c_finite_differences(const MassMeshPointersC* c) {
   mass_fd_impl(c);
}
extern "C" void mass_c_density(double T, double DT, const MassMeshPointersC* c) {
   density_impl(T, DT, c);
}
