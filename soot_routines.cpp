#include "soot_routines.h"
#include <cmath>
#include <algorithm>
#include <numeric>

// =============================================================================
// 模块级全局变量定义
// =============================================================================
namespace soot_routines {

double MIN_AGGLOMERATION = 1.0e-4;

std::vector<double> BIN_S;
std::vector< std::vector<double> > BIN_M;
std::vector< std::vector<double> > BIN_X;
std::vector< std::vector<double> > MOBILITY_FAC;
std::vector< std::vector<double> > A_FAC;
std::vector< std::vector<double> > PARTICLE_RADIUS;

std::vector< std::vector< std::vector<double> > > PHI_B_FAC;
std::vector< std::vector< std::vector<double> > > PHI_G_FAC;
std::vector< std::vector< std::vector<double> > > PHI_S_FAC;
std::vector< std::vector< std::vector<double> > > PHI_I_FAC;
std::vector< std::vector< std::vector<double> > > FU1_FAC;
std::vector< std::vector< std::vector<double> > > FU2_FAC;
std::vector< std::vector< std::vector<double> > > PARTICLE_MASS;

std::vector< std::vector< std::vector< std::vector<double> > > > BIN_ETA;
std::vector< std::vector< std::vector< std::vector<int> > > > BIN_ETA_INDEX;

} // namespace

// =============================================================================
// Cunningham slip correction factor (matches Fortran PHYSICAL_FUNCTIONS)
// =============================================================================
static double cunningham(double KN) {
    constexpr double K1 = 1.257, K2 = 0.4, K3 = 1.1;
    return 1.0 + K1 * KN + K2 * KN * std::exp(-K3 / KN);
}

// =============================================================================
// PBAR lookup with PRESSURE_ZONE — matches Fortran PBAR(K,PRESSURE_ZONE(I,J,K))
// Fortran PBAR(0:KBAR, 1:N_ZONE) is column-major, flat array PBAR_RAW
// =============================================================================
static double pbar_at(double* PBAR_RAW, int PBAR_NK, int K, int PZ) {
    if (!PBAR_RAW) return 101325.0; // fallback
    // column-major: PBAR_RAW[(zone-1)*PBAR_NK + K]
    // Fortran is 1-based for zone, K is 0-based in C++
    return PBAR_RAW[(PZ - 1) * PBAR_NK + K];
}

// =============================================================================
// Fortran → C++ callbacks for physical functions
// These call the Fortran BIND(C) functions defined in soot.f90
// =============================================================================

// ---- callers that pack ZZ into a contiguous buffer ----
static void get_viscosity_cb(const double* ZZ_GET, int nts, double& MU_OUT, double TMP_G) {
    // ZZ_GET is already a contiguous array, pass C pointer to Fortran
    fds_get_viscosity(ZZ_GET, &MU_OUT, TMP_G);
}

static void get_conductivity_cb(const double* ZZ_GET, int nts, double& K_OUT, double TMP_G) {
    fds_get_conductivity(ZZ_GET, &K_OUT, TMP_G);
}

static void get_mass_fraction_cb(const double* ZZ_GET, int nts, int idx, double& Y_OUT) {
    fds_get_mass_fraction(ZZ_GET, idx, &Y_OUT);
}

static void get_molecular_weight_cb(const double* ZZ_GET, int nts, double& MW_OUT) {
    fds_get_molecular_weight(ZZ_GET, &MW_OUT);
}

static void get_specific_heat_cb(const double* ZZ_GET, int nts, double& CP_OUT, double TMP_G) {
    fds_get_specific_heat(ZZ_GET, &CP_OUT, TMP_G);
}

static void get_specific_gas_constant_cb(const double* ZZ_GET, int nts, double& RSUM_OUT) {
    fds_get_specific_gas_constant(ZZ_GET, &RSUM_OUT);
}

static void get_average_specific_heat_cb(const double* ZZ_GET, int nts, double& CPBAR_OUT, double TMP_G) {
    fds_get_average_specific_heat(ZZ_GET, &CPBAR_OUT, TMP_G);
}

static void d_z_lookup_cb(double TMP_VAL, int NS_VAL, double& D_OUT) {
    fds_d_z_lookup(TMP_VAL, NS_VAL, &D_OUT);
}

// =============================================================================
// INITIALIZE_AGGLOMERATION — matches Fortran INITIALIZE_AGGLOMERATION
// =============================================================================
void soot_routines::initialize_agglomeration(const GlobalConstants& gc) {
    using namespace soot_routines;
    int N_SPEC = gc.N_AGGLOMERATION_SPECIES;
    int MAX_BINS = 0;
    for (int n = 0; n < N_SPEC; ++n)
        if (gc.N_PARTICLE_BINS[n] > MAX_BINS) MAX_BINS = gc.N_PARTICLE_BINS[n];

    BIN_S.resize(N_SPEC);
    BIN_M.assign(N_SPEC, std::vector<double>(MAX_BINS + 1, 0.0));
    BIN_X.assign(N_SPEC, std::vector<double>(MAX_BINS, 0.0));
    BIN_ETA.assign(N_SPEC, std::vector< std::vector< std::vector<double> > >(
        MAX_BINS, std::vector< std::vector<double> >(MAX_BINS, std::vector<double>(2, 0.0))));
    BIN_ETA_INDEX.assign(N_SPEC, std::vector< std::vector< std::vector<int> > >(
        MAX_BINS, std::vector< std::vector<int> >(MAX_BINS, std::vector<int>(2, -1))));

    auto a3 = [&](auto& a){ a.assign(N_SPEC, std::vector< std::vector<double> >(
        MAX_BINS, std::vector<double>(MAX_BINS, 0.0))); };
    a3(PHI_B_FAC); a3(PARTICLE_MASS); a3(PHI_G_FAC); a3(PHI_S_FAC); a3(PHI_I_FAC); a3(FU1_FAC); a3(FU2_FAC);
    MOBILITY_FAC.assign(N_SPEC, std::vector<double>(MAX_BINS,0.0));
    A_FAC.assign(N_SPEC, std::vector<double>(MAX_BINS,0.0));
    PARTICLE_RADIUS.assign(N_SPEC, std::vector<double>(MAX_BINS,0.0));

    for (int N = 0; N < N_SPEC; ++N) {
        int NPB = gc.N_PARTICLE_BINS[N];
        double DENS_S = gc.SPECIES[gc.AGGLOMERATION_SPEC_INDEX[N]-1]->DENSITY_SOLID;
        double MIN_MASS = 0.125*FOTHPI*DENS_S*std::pow(gc.MIN_PARTICLE_DIAMETER[N],3);
        double MAX_MASS = 0.125*FOTHPI*DENS_S*std::pow(gc.MAX_PARTICLE_DIAMETER[N],3);
        BIN_S[N] = std::pow(MAX_MASS/MIN_MASS, 1.0/NPB);
        BIN_M[N][0] = MIN_MASS;
        for (int I=1; I<=NPB; ++I) BIN_M[N][I] = BIN_M[N][I-1]*BIN_S[N];
        for (int I=0; I<NPB; ++I) BIN_X[N][I] = 2.0*BIN_M[N][I+1]/(1.0+BIN_S[N]);
        for (int I=0; I<NPB; ++I) {
            PARTICLE_RADIUS[N][I] = std::pow(BIN_X[N][I]/FOTHPI/DENS_S, ONTH);
            MOBILITY_FAC[N][I] = 1.0/(6.0*PI*PARTICLE_RADIUS[N][I]);
            A_FAC[N][I] = std::sqrt(2.0*K_BOLTZMANN*BIN_X[N][I]/PI);
        }
        for (int I=0; I<NPB; ++I) for (int II=0; II<NPB; ++II) {
            PARTICLE_MASS[N][I][II] = BIN_X[N][I]+BIN_X[N][II];
            double r_min = std::min(PARTICLE_RADIUS[N][I],PARTICLE_RADIUS[N][II]);
            double r_sum = PARTICLE_RADIUS[N][I]+PARTICLE_RADIUS[N][II];
            double E_PK = r_min*r_min/(2.0*r_sum*r_sum);
            PHI_G_FAC[N][I][II] = gc.GRAVITATIONAL_SETTLING ? E_PK*r_sum*r_sum : 0.0;
            PHI_B_FAC[N][I][II] = 4.0*PI*K_BOLTZMANN*r_sum;
            PHI_S_FAC[N][I][II] = E_PK*r_sum*r_sum*r_sum*std::sqrt(8.0*PI/15.0);
            PHI_I_FAC[N][I][II] = E_PK*r_sum*r_sum*std::pow(512.0*PI*PI*PI/15.0,0.25);
            FU1_FAC[N][I][II] = r_sum/K_BOLTZMANN*std::sqrt(8.0*K_BOLTZMANN/PI*(1.0/BIN_X[N][I]+1.0/BIN_X[N][II]));
            FU2_FAC[N][I][II] = 2.0/r_sum;
            double PM=PARTICLE_MASS[N][I][II], BX_last=BIN_X[N][NPB-1];
            bool found=false;
            for (int III=1; III<NPB; ++III) {
                if (PM>BX_last) {
                    BIN_ETA_INDEX[N][I][II][0]=NPB; BIN_ETA_INDEX[N][I][II][1]=NPB;
                    BIN_ETA[N][I][II][0]=0.5*PM/BX_last; BIN_ETA[N][I][II][1]=0.5*PM/BX_last;
                    found=true; break;
                } else if (PM>BIN_X[N][III-1] && PM<BIN_X[N][III]) {
                    double bx_lo=BIN_X[N][III-1], bx_hi=BIN_X[N][III];
                    BIN_ETA_INDEX[N][I][II][0]=III; BIN_ETA[N][I][II][0]=(bx_hi-PM)/(bx_hi-bx_lo);
                    BIN_ETA_INDEX[N][I][II][1]=III+1; BIN_ETA[N][I][II][1]=(PM-bx_lo)/(bx_hi-bx_lo);
                    found=true; break;
                }
            }
            if (!found && PM<=BIN_X[N][0]) {
                BIN_ETA_INDEX[N][I][II][0]=1; BIN_ETA[N][I][II][0]=1.0;
                BIN_ETA_INDEX[N][I][II][1]=1; BIN_ETA[N][I][II][1]=0.0;
            }
        }
        int last=NPB-1;
        BIN_ETA_INDEX[N][last][last][0]=NPB; BIN_ETA[N][last][last][0]=1.0;
        BIN_ETA_INDEX[N][last][last][1]=NPB; BIN_ETA[N][last][last][1]=0.0;
    }
}

// =============================================================================
// SETTLING_VELOCITY — matches Fortran SETTLING_VELOCITY exactly
// =============================================================================
void soot_routines::settling_velocity(int, MeshPointers& mp, const GlobalConstants& gc) {
    const int& NTS = gc.N_TRACKED_SPECIES;
    const auto& SMV = gc.SPECIES_MIXTURE;
    const auto& GVEC = gc.GVEC;
    const int& NB = gc.NULL_BOUNDARY, &OB = gc.OPEN_BOUNDARY, &IB = gc.INTERPOLATED_BOUNDARY;
    const bool& GS = gc.GRAVITATIONAL_SETTLING, &TS = gc.THERMOPHORETIC_SETTLING;
    constexpr double CS=1.17, CT=2.2, CM=1.146, CM3=3.0*CM, CS2=CS*2.0, CT2=2.0*CT;
    int I,J,K,N,IW,IIG,JJG,KKG,IOR;

    auto* RHOP = mp.PREDICTOR ? mp.RHO : mp.RHOS;
    auto& US=*mp.WORK7, &VS=*mp.WORK8, &WS=*mp.WORK9;
    int IBAR=*mp.IBAR, JBAR=*mp.JBAR, KBAR=*mp.KBAR, NW=*mp.N_EXTERNAL_WALL_CELLS+*mp.N_INTERNAL_WALL_CELLS;

    for (N=1; N<=NTS; ++N) {
        const auto& SM = *SMV[N-1];
        if (!SM.DEPOSITING) continue;
        US.assign(0.0); VS.assign(0.0); WS.assign(0.0);
        double GF = GS ? SM.MEAN_DIAMETER*SM.MEAN_DIAMETER*SM.DENSITY_SOLID/18.0 : 0.0;

        for (K=0; K<=KBAR; ++K) for (J=0; J<=JBAR; ++J) for (I=0; I<=IBAR; ++I) {
            if ((*mp.CELL)[(*mp.CELL_INDEX)(I+1,J+1,K+1)-1].SOLID) continue;

            double TG,RG,PG,DTDN,MUG,KNF,KN,CN,ALPHA,KG,GV;
            int PZ;
            std::vector<double> ZZ_GET(NTS+1); // 1-based

            // ---- U direction (X face: I direction) ----
            auto* ZZP = mp.PREDICTOR ? mp.ZZ : mp.ZZS;
            TG = 0.5*((*mp.TMP)(I+1,J+1,K+1)+(*mp.TMP)(I+2,J+1,K+1));
            for (int ns=1; ns<=NTS; ++ns)
                ZZ_GET[ns] = 0.5*((*ZZP)(I+1,J+1,K+1,ns)+(*ZZP)(I+2,J+1,K+1,ns));
            RG = 0.5*((*RHOP)(I+1,J+1,K+1)+(*RHOP)(I+2,J+1,K+1));
            PZ = mp.PRESSURE_ZONE ? (*mp.PRESSURE_ZONE)(I+1,J+1,K+1) : 1;
            PG = 0.5*(pbar_at(mp.PBAR_RAW, mp.PBAR_NK, K, PZ) +
                      pbar_at(mp.PBAR_RAW, mp.PBAR_NK, K, mp.PRESSURE_ZONE ? (*mp.PRESSURE_ZONE)(I+2,J+1,K+1) : 1));
           DTDN = -((*mp.TMP)(I+2,J+1,K+1)-(*mp.TMP)(I+1,J+1,K+1))*(*mp.RDXN)[I];
            get_viscosity_cb(&ZZ_GET[1], NTS, MUG, TG);
            KNF = std::sqrt(2.0*PI/(PG*RG))*MUG;
            if (TS) {
                KN = KNF/SM.THERMOPHORETIC_DIAMETER; CN = cunningham(KN);
                get_conductivity_cb(&ZZ_GET[1], NTS, KG, TG);
                ALPHA = KG/SM.CONDUCTIVITY_SOLID;
                US(I+1,J+1,K+1) += CS2*(ALPHA+CT*KN)*CN/((1.0+CM3*KN)*(1.0+2.0*ALPHA+CT2*KN))*MUG/(TG*RG)*DTDN;
            }
            if (GS) {
                KN = KNF/SM.MEAN_DIAMETER; CN = cunningham(KN);
                GV = GF*CN/MUG;
                US(I+1,J+1,K+1) += GVEC[0]*GV;
            }

            // ---- V direction (Y face: J direction) ----
            TG = 0.5*((*mp.TMP)(I+1,J+1,K+1)+(*mp.TMP)(I+1,J+2,K+1));
            for (int ns=1; ns<=NTS; ++ns)
                ZZ_GET[ns] = 0.5*((*ZZP)(I+1,J+1,K+1,ns)+(*ZZP)(I+1,J+2,K+1,ns));
            RG = 0.5*((*RHOP)(I+1,J+1,K+1)+(*RHOP)(I+1,J+2,K+1));
            PZ = mp.PRESSURE_ZONE ? (*mp.PRESSURE_ZONE)(I+1,J+1,K+1) : 1;
            PG = 0.5*(pbar_at(mp.PBAR_RAW, mp.PBAR_NK, K, (*mp.PRESSURE_ZONE)(I+1,J+1,K+1)) +
                      pbar_at(mp.PBAR_RAW, mp.PBAR_NK, K, (*mp.PRESSURE_ZONE)(I+1,J+2,K+1)));
            DTDN = -((*mp.TMP)(I+1,J+2,K+1)-(*mp.TMP)(I+1,J+1,K+1))*(*mp.RDYN)[J];
            get_viscosity_cb(&ZZ_GET[1], NTS, MUG, TG);
            KNF = std::sqrt(2.0*PI/(PG*RG))*MUG;
            if (TS) {
                KN = KNF/SM.THERMOPHORETIC_DIAMETER; CN = cunningham(KN);
                get_conductivity_cb(&ZZ_GET[1], NTS, KG, TG);
                ALPHA = KG/SM.CONDUCTIVITY_SOLID;
                VS(I+1,J+1,K+1) += CS2*(ALPHA+CT*KN)*CN/((1.0+CM3*KN)*(1.0+2.0*ALPHA+CT2*KN))*MUG/(TG*RG)*DTDN;
            }
            if (GS) {
                KN = KNF/SM.MEAN_DIAMETER; CN = cunningham(KN);
                GV = GF*CN/MUG;
                VS(I+1,J+1,K+1) += GVEC[1]*GV;
            }

            // ---- W direction (Z face: K direction) ----
            TG = 0.5*((*mp.TMP)(I+1,J+1,K+1)+(*mp.TMP)(I+1,J+1,K+2));
            for (int ns=1; ns<=NTS; ++ns)
                ZZ_GET[ns] = 0.5*((*ZZP)(I+1,J+1,K+1,ns)+(*ZZP)(I+1,J+1,K+2,ns));
            RG = 0.5*((*RHOP)(I+1,J+1,K+1)+(*RHOP)(I+1,J+1,K+2));
            PZ = mp.PRESSURE_ZONE ? (*mp.PRESSURE_ZONE)(I+1,J+1,K+1) : 1;
            PG = 0.5*(pbar_at(mp.PBAR_RAW, mp.PBAR_NK, K, (*mp.PRESSURE_ZONE)(I+1,J+1,K+1)) +
                      pbar_at(mp.PBAR_RAW, mp.PBAR_NK, K, (*mp.PRESSURE_ZONE)(I+1,J+1,K+2)));
            DTDN = -((*mp.TMP)(I+1,J+1,K+2)-(*mp.TMP)(I+1,J+1,K+1))*(*mp.RDZN)[K];
            get_viscosity_cb(&ZZ_GET[1], NTS, MUG, TG);
            KNF = std::sqrt(2.0*PI/(PG*RG))*MUG;
            if (TS) {
                KN = KNF/SM.THERMOPHORETIC_DIAMETER; CN = cunningham(KN);
                get_conductivity_cb(&ZZ_GET[1], NTS, KG, TG);
                ALPHA = KG/SM.CONDUCTIVITY_SOLID;
                WS(I+1,J+1,K+1) += CS2*(ALPHA+CT*KN)*CN/((1.0+CM3*KN)*(1.0+2.0*ALPHA+CT2*KN))*MUG/(TG*RG)*DTDN;
            }
            if (GS) {
                KN = KNF/SM.MEAN_DIAMETER; CN = cunningham(KN);
                GV = GF*CN/MUG;
                WS(I+1,J+1,K+1) += GVEC[2]*GV;
            }
        }

        // Wall BC — zero settling velocity at walls
        for (IW=1; IW<=NW; ++IW) {
            auto& WC = (*mp.WALL)[IW-1]; auto& BC = (*mp.BOUNDARY_COORD)[WC.BC_INDEX-1];
            if (WC.BOUNDARY_TYPE==NB||WC.BOUNDARY_TYPE==OB||WC.BOUNDARY_TYPE==IB) continue;
            IIG=BC.IIG; JJG=BC.JJG; KKG=BC.KKG; IOR=BC.IOR;
            switch(IOR) {
                case -1: US(IIG,JJG,KKG)=0; break; case 1: US(IIG-1,JJG,KKG)=0; break;
                case -2: VS(IIG,JJG,KKG)=0; break; case 2: VS(IIG,JJG-1,KKG)=0; break;
                case -3: WS(IIG,JJG,KKG)=0; break; case 3: WS(IIG,JJG,KKG-1)=0; break;
            }
        }

        // Divergence term for settling
        for (K=1; K<=KBAR; ++K) for (J=1; J<=JBAR; ++J) for (I=1; I<=IBAR; ++I) {
            if ((*mp.CELL)[(*mp.CELL_INDEX)(I,J,K)-1].SOLID) continue;
            double RHS = ((*mp.R)[I]*(*mp.FX)(I,J,K,N)*US(I+1,J+1,K+1)-(*mp.R)[I-1]*(*mp.FX)(I-1,J,K,N)*US(I,J+1,K+1))*(*mp.RDX)[I]*(*mp.RRN)[I]
                + ((*mp.FY)(I,J,K,N)*VS(I+1,J+1,K+1)-(*mp.FY)(I,J-1,K,N)*VS(I+1,J,K+1))*(*mp.RDY)[J]
                + ((*mp.FZ)(I,J,K,N)*WS(I+1,J+1,K+1)-(*mp.FZ)(I,J,K-1,N)*WS(I+1,J+1,K))*(*mp.RDZ)[K];
            (*mp.DEL_RHO_D_DEL_Z)(I,J,K,N) -= RHS;
        }
    }
}

// =============================================================================
// CALC_AGGLOMERATION — matches Fortran CALC_AGGLOMERATION exactly
// =============================================================================
void soot_routines::calc_agglomeration(double DT, int, MeshPointers& mp, const GlobalConstants& gc) {
    const auto& NBINS = gc.N_PARTICLE_BINS;
    const auto& SMIX = gc.AGGLOMERATION_SMIX_INDEX;
    using namespace soot_routines;
    int I,J,K,NS,N,NN,IM1,IM2,JM1,JM2,KM1,KM2,IP1,JP1,KP1;

    for (NS=1; NS<=gc.N_AGGLOMERATION_SPECIES; ++NS) {
        for (K=1; K<=*mp.KBAR; ++K) for (J=1; J<=*mp.JBAR; ++J) for (I=1; I<=*mp.IBAR; ++I) {
            if ((*mp.CELL)[(*mp.CELL_INDEX)(I,J,K)-1].SOLID) continue;
            int NPB = NBINS[NS-1];
            std::vector<double> N0(NPB+1,0.0);
            for (int idx=1; idx<=NPB; ++idx) N0[idx]=(*mp.ZZ)(I,J,K,SMIX[NS-1]+idx-1);
            double RHO_G=(*mp.RHO)(I,J,K);
            for (int idx=1; idx<=NPB; ++idx) N0[idx]*=RHO_G/BIN_X[NS-1][idx-1];
            bool all_small=true;
            for (int idx=1; idx<=NPB; ++idx) if(N0[idx]>=MIN_AGGLOMERATION){all_small=false;break;}
            if(all_small) continue;
            double TMP_G=(*mp.TMP)(I,J,K), MUG;
            // ZZ_GET = ZZ(I,J,K,1:N_TRACKED_SPECIES) for property lookup
            std::vector<double> ZZ_GET(gc.N_TRACKED_SPECIES+1);
            for (int ns=1; ns<=gc.N_TRACKED_SPECIES; ++ns) ZZ_GET[ns]=(*mp.ZZ)(I,J,K,ns);
            get_viscosity_cb(&ZZ_GET[1], gc.N_TRACKED_SPECIES, MUG, TMP_G);
            // PBAR(K, PRESSURE_ZONE(I,J,K)) — zone-aware lookup
            double PBAR_VAL;
            if (mp.PRESSURE_ZONE && mp.PBAR_RAW) {
                int PZ = (*mp.PRESSURE_ZONE)(I,J,K);
                PBAR_VAL = pbar_at(mp.PBAR_RAW, mp.PBAR_NK, K-1, PZ);
            } else if (mp.PBAR) {
                // Single-zone fallback
                PBAR_VAL = (*mp.PBAR)[K];
            } else {
                PBAR_VAL = 101325.0;
            }
            double KN_FAC=MUG*std::sqrt(PI/(2.0*PBAR_VAL*RHO_G));
            IM1=std::max(0,I-1); JM1=std::max(0,J-1); KM1=std::max(0,K-1);
            IM2=std::max(1,I-1); JM2=std::max(1,J-1); KM2=std::max(1,K-1);
            IP1=std::min(*mp.IBAR,I+1); JP1=std::min(*mp.JBAR,J+1); KP1=std::min(*mp.KBAR,K+1);
            double DUDX=(*mp.RDX)[I]*((*mp.U)(I+1,J+1,K+1)-(*mp.U)(IM1+1,J+1,K+1));
            double DVDY=(*mp.RDY)[J]*((*mp.V)(I+1,J+1,K+1)-(*mp.V)(I+1,JM1+1,K+1));
            double DWDZ=(*mp.RDZ)[K]*((*mp.W)(I+1,J+1,K+1)-(*mp.W)(I+1,J+1,KM1+1));
            double OD=ONTH*(DUDX+DVDY+DWDZ), S11=DUDX-OD, S22=DVDY-OD, S33=DWDZ-OD;
            double DUDY=0.25*(*mp.RDY)[J]*((*mp.U)(I+1,JP1+1,K+1)-(*mp.U)(I+1,JM2+1,K+1)+(*mp.U)(IM1+1,JP1+1,K+1)-(*mp.U)(IM1+1,JM2+1,K+1));
            double DUDZ=0.25*(*mp.RDZ)[K]*((*mp.U)(I+1,J+1,KP1+1)-(*mp.U)(I+1,J+1,KM2+1)+(*mp.U)(IM1+1,J+1,KP1+1)-(*mp.U)(IM1+1,J+1,KM2+1));
            double DVDX=0.25*(*mp.RDX)[I]*((*mp.V)(IP1+1,J+1,K+1)-(*mp.V)(IM2+1,J+1,K+1)+(*mp.V)(IP1+1,JM1+1,K+1)-(*mp.V)(IM2+1,JM1+1,K+1));
            double DVDZ=0.25*(*mp.RDZ)[K]*((*mp.V)(I+1,J+1,KP1+1)-(*mp.V)(I+1,J+1,KM2+1)+(*mp.V)(I+1,JM1+1,KP1+1)-(*mp.V)(I+1,JM1+1,KM2+1));
            double DWDX=0.25*(*mp.RDX)[I]*((*mp.W)(IP1+1,J+1,K+1)-(*mp.W)(IM2+1,J+1,K+1)+(*mp.W)(IP1+1,J+1,KM1+1)-(*mp.W)(IM2+1,J+1,KM1+1));
            double DWDY=0.25*(*mp.RDY)[J]*((*mp.W)(I+1,JP1+1,K+1)-(*mp.W)(I+1,JM2+1,K+1)+(*mp.W)(I+1,JP1+1,KM1+1)-(*mp.W)(I+1,JM2+1,KM1+1));
            double SR=2.0*(S11*S11+S22*S22+S33*S33+2.0*(std::pow(0.5*(DUDY+DVDX),2)+std::pow(0.5*(DUDZ+DWDX),2)+std::pow(0.5*(DVDZ+DWDY),2)));
            double DR=(*mp.MU)(I,J,K)/RHO_G*SR;
            // Wait: Fortran uses MU(I,J,K)/RHO_G*STRAIN_RATE, correct
            auto NI=N0;
            std::vector<double> MOB(NPB+1), TERM(NPB+1), AMT(NPB+1);
            std::vector< std::vector<double> > PHI(NPB+1, std::vector<double>(NPB+1,0.0));
            for (N=1; N<=NPB; ++N) {
                double KN=KN_FAC/PARTICLE_RADIUS[NS-1][N-1];
                MOB[N]=cunningham(KN)*MOBILITY_FAC[NS-1][N-1]/MUG;
                TERM[N]=MOB[N]*GRAV*BIN_X[NS-1][N-1];
                double AM=A_FAC[NS-1][N-1]*std::sqrt(TMP_G)*MOB[N];
                AMT[N]=(std::pow(PARTICLE_RADIUS[NS-1][N-1]+AM,3)-std::pow(PARTICLE_RADIUS[NS-1][N-1]*PARTICLE_RADIUS[NS-1][N-1]+AM*AM,1.5))/(3.0*PARTICLE_RADIUS[NS-1][N-1]*AM)-PARTICLE_RADIUS[NS-1][N-1];
            }
            for (N=1; N<=NPB; ++N) for (NN=1; NN<=NPB; ++NN) {
                if (NN<N) continue;
                double FU1=FU1_FAC[NS-1][N-1][NN-1]/(std::sqrt(TMP_G)*(MOB[N]+MOB[NN]));
                double FU2=1.0+FU2_FAC[NS-1][N-1][NN-1]*std::sqrt(AMT[NN]*AMT[NN]+AMT[N]*AMT[N]);
                double FU=1.0/(1.0/FU1+1.0/FU2);
                double PB=PHI_B_FAC[NS-1][N-1][NN-1]*(MOB[N]+MOB[NN])*FU*TMP_G;
                double VR=std::abs(TERM[N]-TERM[NN]);
                double PG=PHI_G_FAC[NS-1][N-1][NN-1]*VR;
                double PS_loc=PHI_S_FAC[NS-1][N-1][NN-1]*std::sqrt(RHO_G/MUG*DR);
                double PII=(GRAV<=TWO_EPSILON_EB)?0.0:PHI_I_FAC[NS-1][N-1][NN-1]*std::pow(RHO_G/MUG*DR*DR*DR,0.25)*VR/GRAV;
                PHI[N][NN]=PB+PG+std::sqrt(PS_loc*PS_loc+PII*PII); PHI[NN][N]=PHI[N][NN];
            }
            double DTS=DT, DSUM=0.0;
            std::vector<double> N3f(NPB+1,0.0);
            while (DSUM<DT) {
                std::vector<double> N1(NPB+1,0.0), N2(NPB+1,0.0), N3(NPB+1,0.0);
                for (N=1; N<=NPB; ++N) for (NN=N; NN<=NPB; ++NN) {
                    if (N0[N]<MIN_AGGLOMERATION||N0[NN]<MIN_AGGLOMERATION) continue;
                    double DN=PHI[NN][N]*N0[N]*N0[NN]*DTS;
                    N1[N]-=DN; N1[NN]-=DN;
                    int i1=BIN_ETA_INDEX[NS-1][N-1][NN-1][0], i2=BIN_ETA_INDEX[NS-1][N-1][NN-1][1];
                    N2[i1]+=BIN_ETA[NS-1][N-1][NN-1][0]*DN; N2[i2]+=BIN_ETA[NS-1][N-1][NN-1][1]*DN;
                }
                for (int idx=1; idx<=NPB; ++idx) N3[idx]=N1[idx]+N2[idx]+N0[idx];
                double TOL=0.0;
                for (int idx=1; idx<=NPB; ++idx) TOL=std::max(TOL,std::abs((N0[idx]-N3[idx])/(N0[idx]+TINY_EB)));
                if (TOL>0.3) DTS*=0.3/TOL;
                else { DSUM+=DTS; DTS=std::min(DT-DSUM,1.5*DTS); N0=N3; }
                N3f=N3;
            }
            double sNI=0, sN3=0;
            for (int idx=1; idx<=NPB; ++idx) { sNI+=NI[idx]*BIN_X[NS-1][idx-1]; sN3+=N3f[idx]*BIN_X[NS-1][idx-1]; }
            if (sN3>0) { double f=sNI/sN3; for (int idx=1; idx<=NPB; ++idx) N3f[idx]*=f; }
            for (int idx=1; idx<=NPB; ++idx) (*mp.ZZ)(I,J,K,SMIX[NS-1]+idx-1)=N3f[idx]*BIN_X[NS-1][idx-1]/RHO_G;
        }
    }
}

// ---- Forward declarations for helpers used in soot_surface_oxidation ----
static double nu_soot_ox_val(const GlobalConstants&, int nts, int ns);
static double min_nu_soot_ox(const GlobalConstants&, int nts);
static bool is_agglomerating(const GlobalConstants&, int soot_index);
static int find_agglomeration_spec_index(const GlobalConstants&, int soot_index);

// =============================================================================
// SOOT_SURFACE_OXIDATION — matches Fortran SOOT_SURFACE_OXIDATION exactly
// Hartman, Beyler, Riahi, Beyler, Fire and Materials (2012), 36:177-184
// =============================================================================
void soot_routines::soot_surface_oxidation(double DT, int, MeshPointers& mp, const GlobalConstants& gc) {
    constexpr double R0=8314.0, E=-2.110e5*1000.0/R0, A=4.7e10;
    const int& NTS = gc.N_TRACKED_SPECIES;
    const auto& SMV = gc.SPECIES_MIXTURE;
    const int& SOOT_INDEX = gc.SOOT_INDEX;
    const int& O2_INDEX = gc.O2_INDEX;
    const double& MW_O2 = gc.MW_O2_VAL;
    const double& ZZ_MIN = gc.ZZ_MIN_GLOBAL_VAL;
    int IW, ICF, NS, NS2, IIG, JJG, KKG;
    // Fortran equivalent: SS => SPECIES(SOOT_INDEX); SS%AWM_INDEX (1-based)
    int SOOT_AWM = (SOOT_INDEX > 0 && SOOT_INDEX <= (int)gc.SPECIES.size())
                   ? gc.SPECIES[SOOT_INDEX-1]->AWM_INDEX : 1;

    // ---- WALL loop ----
    for (IW=1; IW<=*mp.N_EXTERNAL_WALL_CELLS+*mp.N_INTERNAL_WALL_CELLS; ++IW) {
        auto& WC=(*mp.WALL)[IW-1];
        auto& B1=(*mp.BOUNDARY_PROP1)[WC.B1_INDEX-1];
        auto& BC=(*mp.BOUNDARY_COORD)[WC.BC_INDEX-1];
        if (WC.BOUNDARY_TYPE!=1) continue;  // SOLID_BOUNDARY=1
        if (!B1.AWM_AEROSOL || B1.AWM_AEROSOL[SOOT_AWM-1] < ZZ_MIN) continue;
        double M_SOOT = B1.AWM_AEROSOL[SOOT_AWM-1] * B1.AREA;
        IIG=BC.IIG; JJG=BC.JJG; KKG=BC.KKG;
        double RHO_G = (*mp.RHO)(IIG,JJG,KKG);
        double TMP_G = (*mp.TMP)(IIG,JJG,KKG);
        // ZZ_GET = MAX(0, ZZ - M_DOT_PPP*DT/RHO_G)
        std::vector<double> ZZ_GET(NTS+1);
        for (int ns=1; ns<=NTS; ++ns)
            ZZ_GET[ns] = std::max(0.0, (*mp.ZZ)(IIG,JJG,KKG,ns) -
                (*mp.M_DOT_PPP)(IIG,JJG,KKG,ns)*DT/RHO_G);
        double Y_O2;
        get_mass_fraction_cb(&ZZ_GET[1], NTS, O2_INDEX, Y_O2);
        if (Y_O2 < ZZ_MIN) continue;
        double MW;
        get_molecular_weight_cb(&ZZ_GET[1], NTS, MW);
        double X_O2 = Y_O2*MW/MW_O2;
        std::vector<double> DZZ(NTS+1, 0.0);
        double DMDT = A * M_SOOT * X_O2 * std::exp(E / B1.TMP_F);
        double VOL = (*mp.DX)[IIG]*(*mp.RC)[IIG]*(*mp.DY)[JJG]*(*mp.DZ)[KKG];
        double DM = std::min({M_SOOT, DMDT*DT, -Y_O2*RHO_G*VOL/min_nu_soot_ox(gc, NTS)});
        // DZZ = NU_SOOT_OX * DM / VOL
        for (int ns=1; ns<=NTS; ++ns) DZZ[ns] = nu_soot_ox_val(gc, NTS, ns) * DM / VOL;
        // Q -= SUM(H_F * DZZ) / DT
        double H_F_sum = 0.0;
        for (int ns=1; ns<=NTS; ++ns) H_F_sum += SMV[ns-1]->H_F * DZZ[ns];
        (*mp.Q)(IIG,JJG,KKG) -= H_F_sum / DT;
        DM /= B1.AREA;
        // Agglomerating species redistribution
        if (is_agglomerating(gc, SOOT_INDEX)) {
            M_SOOT /= B1.AREA;
            NS2 = find_agglomeration_spec_index(gc, SOOT_INDEX);
            for (int bin_idx=0; bin_idx<gc.N_PARTICLE_BINS[NS2-1]; ++bin_idx) {
                int sm_idx = gc.AGGLOMERATION_SMIX_INDEX[NS2-1] + bin_idx;
                auto& SM_BIN = *SMV[sm_idx-1];
                if (SM_BIN.AWM_INDEX > 0)
                    B1.AWM_AEROSOL[SM_BIN.AWM_INDEX-1] -= DM * B1.AWM_AEROSOL[SM_BIN.AWM_INDEX-1] / M_SOOT;
            }
        }
        B1.AWM_AEROSOL[SOOT_AWM-1] -= DM;
        // Divergence term
        double CP, H_G, RSUM_LOC, CPBAR, CPBAR2, MW_RATIO, DELTA_H_G, M_DOT_PPP_SINGLE;
        get_specific_heat_cb(&ZZ_GET[1], NTS, CP, TMP_G);
        H_G = CP * TMP_G;
        get_specific_gas_constant_cb(&ZZ_GET[1], NTS, RSUM_LOC);
        for (NS=1; NS<=NTS; ++NS) {
            if (std::abs(DZZ[NS]) < TWO_EPSILON_EB) continue;
            std::vector<double> ZZ_PURE(NTS+1, 0.0);
            ZZ_PURE[NS] = 1.0;
            M_DOT_PPP_SINGLE = DZZ[NS] / DT;
            MW_RATIO = SMV[NS-1]->RCON / (*mp.RSUM)(IIG,JJG,KKG);
            if (DZZ[NS] < 0.0) {
                (*mp.D_SOURCE)(IIG,JJG,KKG) += M_DOT_PPP_SINGLE * MW_RATIO / RHO_G;
            } else {
                get_average_specific_heat_cb(&ZZ_PURE[1], NTS, CPBAR, TMP_G);
                get_average_specific_heat_cb(&ZZ_PURE[1], NTS, CPBAR2, B1.TMP_F);
                DELTA_H_G = CPBAR2*B1.TMP_F - CPBAR*TMP_G;
                (*mp.D_SOURCE)(IIG,JJG,KKG) += M_DOT_PPP_SINGLE * (MW_RATIO + DELTA_H_G/H_G) / RHO_G;
                (*mp.M_DOT_PPP)(IIG,JJG,KKG,NS) += M_DOT_PPP_SINGLE;
            }
        }
    }

    // ---- CFACE loop ----
    for (ICF=*mp.INTERNAL_CFACE_CELLS_LB+1; ICF<=*mp.INTERNAL_CFACE_CELLS_LB+*mp.N_INTERNAL_CFACE_CELLS; ++ICF) {
        int cf_idx = ICF - *mp.INTERNAL_CFACE_CELLS_LB - 1;  // 0-based
        if (cf_idx < 0 || cf_idx >= (int)mp.CFACE->size()) continue;
        auto& CFA = (*mp.CFACE)[cf_idx];
        auto& B1 = (*mp.BOUNDARY_PROP1)[CFA.B1_INDEX-1];
        auto& BC = (*mp.BOUNDARY_COORD)[CFA.BC_INDEX-1];
        if (CFA.BOUNDARY_TYPE!=1) continue;
        if (!B1.AWM_AEROSOL || B1.AWM_AEROSOL[SOOT_AWM-1] < ZZ_MIN) continue;
        double M_SOOT = B1.AWM_AEROSOL[SOOT_AWM-1] * B1.AREA;
        IIG=BC.IIG; JJG=BC.JJG; KKG=BC.KKG;
        double RHO_G = B1.RHO_G;  // CFACE uses B1%RHO_G
        double TMP_G = B1.TMP_G;  // CFACE uses B1%TMP_G
        std::vector<double> ZZ_GET(NTS+1);
        if (B1.ZZ_G) {
            for (int ns=1; ns<=NTS; ++ns)
                ZZ_GET[ns] = std::max(0.0, B1.ZZ_G[ns-1] -
                    (*mp.M_DOT_PPP)(IIG,JJG,KKG,ns)*DT/RHO_G);
        } else {
            for (int ns=1; ns<=NTS; ++ns)
                ZZ_GET[ns] = std::max(0.0, (*mp.ZZ)(IIG,JJG,KKG,ns) -
                    (*mp.M_DOT_PPP)(IIG,JJG,KKG,ns)*DT/RHO_G);
        }
        double Y_O2;
        get_mass_fraction_cb(&ZZ_GET[1], NTS, O2_INDEX, Y_O2);
        if (Y_O2 < ZZ_MIN) continue;
        double MW;
        get_molecular_weight_cb(&ZZ_GET[1], NTS, MW);
        double X_O2 = Y_O2*MW/MW_O2;
        std::vector<double> DZZ(NTS+1, 0.0);
        double DMDT = A * M_SOOT * X_O2 * std::exp(E / B1.TMP_F);
        double VOL = (*mp.DX)[IIG]*(*mp.RC)[IIG]*(*mp.DY)[JJG]*(*mp.DZ)[KKG];
        double DM = std::min({M_SOOT, DMDT*DT, -Y_O2*RHO_G*VOL/min_nu_soot_ox(gc, NTS)});
        for (int ns=1; ns<=NTS; ++ns) DZZ[ns] = nu_soot_ox_val(gc, NTS, ns) * DM / VOL;
        double H_F_sum = 0.0;
        for (int ns=1; ns<=NTS; ++ns) H_F_sum += SMV[ns-1]->H_F * DZZ[ns];
        (*mp.Q)(IIG,JJG,KKG) -= H_F_sum / DT;
        DM /= B1.AREA;
        if (is_agglomerating(gc, SOOT_INDEX)) {
            M_SOOT /= B1.AREA;
            NS2 = find_agglomeration_spec_index(gc, SOOT_INDEX);
            for (int bin_idx=0; bin_idx<gc.N_PARTICLE_BINS[NS2-1]; ++bin_idx) {
                int sm_idx = gc.AGGLOMERATION_SMIX_INDEX[NS2-1] + bin_idx;
                auto& SM_BIN = *SMV[sm_idx-1];
                if (SM_BIN.AWM_INDEX > 0)
                    B1.AWM_AEROSOL[SM_BIN.AWM_INDEX-1] -= DM * B1.AWM_AEROSOL[SM_BIN.AWM_INDEX-1] / M_SOOT;
            }
        }
        B1.AWM_AEROSOL[SOOT_AWM-1] -= DM;
        // Divergence term for CFACE
        {
        double CP, H_G, RSUM_LOC, CPBAR, CPBAR2, MW_RATIO, DELTA_H_G, M_DOT_PPP_SINGLE;
        get_specific_heat_cb(&ZZ_GET[1], NTS, CP, TMP_G);
        H_G = CP * TMP_G;
        get_specific_gas_constant_cb(&ZZ_GET[1], NTS, RSUM_LOC);
        for (NS=1; NS<=NTS; ++NS) {
            if (std::abs(DZZ[NS]) < TWO_EPSILON_EB) continue;
            std::vector<double> ZZ_PURE(NTS+1, 0.0);
            ZZ_PURE[NS] = 1.0;
            M_DOT_PPP_SINGLE = DZZ[NS] / DT;
            MW_RATIO = SMV[NS-1]->RCON / CFA.RSUM_G;  // CFACE uses CFA%RSUM_G
            if (DZZ[NS] < 0.0) {
                (*mp.D_SOURCE)(IIG,JJG,KKG) += M_DOT_PPP_SINGLE * MW_RATIO / RHO_G;
            } else {
                get_average_specific_heat_cb(&ZZ_PURE[1], NTS, CPBAR, TMP_G);
                get_average_specific_heat_cb(&ZZ_PURE[1], NTS, CPBAR2, B1.TMP_F);
                DELTA_H_G = CPBAR2*B1.TMP_F - CPBAR*TMP_G;
                (*mp.D_SOURCE)(IIG,JJG,KKG) += M_DOT_PPP_SINGLE * (MW_RATIO + DELTA_H_G/H_G) / RHO_G;
                (*mp.M_DOT_PPP)(IIG,JJG,KKG,NS) += M_DOT_PPP_SINGLE;
            }
            }
        }
    }
}

// =============================================================================
// DROPLET_SCRUBBING — matches Fortran DROPLET_SCRUBBING exactly
// =============================================================================
void soot_routines::droplet_scrubbing(int IP, int, double DT, double DT_P, MeshPointers& mp, const GlobalConstants& gc) {
    const int& NTS = gc.N_TRACKED_SPECIES;
    const auto& SMV = gc.SPECIES_MIXTURE;
    
    int IIG=1, JJG=1, KKG=1;
    
    double LP_U=0, LP_V=0, LP_W=0, LP_RAD=1e-4, LP_RVC=1.0, LP_PWT=1.0;
    
    if (!mp.BOUNDARY_COORD->empty()) {
        IIG = (*mp.BOUNDARY_COORD)[0].IIG;
        JJG = (*mp.BOUNDARY_COORD)[0].JJG;
        KKG = (*mp.BOUNDARY_COORD)[0].KKG;
    }
    if (!mp.BOUNDARY_PROP1->empty()) {
        auto& B1 = (*mp.BOUNDARY_PROP1)[0];
        LP_RAD = 0.5e-3; 
        
    }
    double VEL = std::sqrt(LP_U*LP_U + LP_V*LP_V + LP_W*LP_W);
    double R_D = LP_RAD;
    double VOL = LP_RVC > 0 ? 1.0/LP_RVC : 0.0;
    double FRAC = std::min(1.0, VEL*DT_P*std::min(std::pow(VOL,TWTH), LP_PWT*PI*R_D*R_D) / (VOL > 0 ? VOL : 1.0));
    double VREL = std::sqrt(std::pow(LP_U-(*mp.U)(IIG,JJG,KKG),2) +
                            std::pow(LP_V-(*mp.V)(IIG,JJG,KKG),2) +
                            std::pow(LP_W-(*mp.W)(IIG,JJG,KKG),2));
    double MU_G;
    // ZZ_GET = ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES)
    std::vector<double> ZZ_GET(NTS+1);
    for (int ns=1; ns<=NTS; ++ns) ZZ_GET[ns] = (*mp.ZZ)(IIG,JJG,KKG,ns);
    get_viscosity_cb(&ZZ_GET[1], NTS, MU_G, (*mp.TMP)(IIG,JJG,KKG));
    double RE = (*mp.RHO)(IIG,JJG,KKG) * VREL * 2.0 * R_D / MU_G;
    for (int NS=1; NS<=NTS; ++NS) {
        auto& SM = *SMV[NS-1];
        if (!SM.DEPOSITING) continue;
        double R_RATIO = 0.5 * SM.MEAN_DIAMETER / R_D;
        double EFF_IN_VIS = std::pow(1.0+R_RATIO,2) * (1.0 - 3.0/(2.0*(1.0+R_RATIO)) + 1.0/(2.0*std::pow(1.0+R_RATIO,3)));
        double EFF_IN_POT = std::pow(1.0+R_RATIO,2) - (1.0+R_RATIO);
        double EFF_IN = (EFF_IN_VIS + EFF_IN_POT * RE/60.0) / (1.0 + RE/60.0);
        double STK = 0.5 * SM.MEAN_DIAMETER * SM.MEAN_DIAMETER * SM.DENSITY_SOLID * VREL / (9.0 * MU_G * R_D);
        double EFF_IM_POT;
        if (STK <= 0.0834) {
            EFF_IM_POT = 0.0;
        } else {
            EFF_IM_POT = std::pow(STK/(STK+0.5), 2);
            if (STK < 0.2) EFF_IM_POT = (STK-0.0834)/(0.2-0.0834) * EFF_IM_POT;
        }
        double EFF_IM_VIS;
        if (STK > 1.214) {
            EFF_IM_VIS = std::pow(1.0 + 0.75*std::log(2.0*STK)/(STK-1.214), -2);
        } else {
            EFF_IM_VIS = 0.0;
        }
        double EFF_IM = (EFF_IM_VIS + EFF_IM_POT * RE/60.0) / (1.0 + RE/60.0);
        double PE, D_OUT;
        d_z_lookup_cb((*mp.TMP)(IIG,JJG,KKG), NS, D_OUT);
        PE = 2.0 * R_D * VREL / D_OUT;
        double EFF = 1.0 - (1.0 - EFF_IN) * (1.0 - EFF_IM);
        double RATE = EFF * FRAC * ZZ_GET[NS] * (*mp.RHO)(IIG,JJG,KKG) / DT;
        (*mp.D_SOURCE)(IIG,JJG,KKG) -= RATE * SMV[0]->MW / SM.MW / (*mp.RHO)(IIG,JJG,KKG);
        (*mp.M_DOT_PPP)(IIG,JJG,KKG,NS) -= RATE;
    }
}


static double nu_soot_ox_val(const GlobalConstants& gc, int nts, int ns) {
    (void)nts;
    if ((int)gc.NU_SOOT_OX.size() >= ns)
        return gc.NU_SOOT_OX[ns-1];
    return -1.0;
}
static double min_nu_soot_ox(const GlobalConstants& gc, int nts) {
    double mn = 1e10;
    for (int i=1; i<=nts; ++i) {
        double v = nu_soot_ox_val(gc, nts, i);
        if (v < mn) mn = v;
    }
    return mn;
}
static bool is_agglomerating(const GlobalConstants& gc, int soot_index) {
    if (soot_index > 0 && soot_index <= (int)gc.SPECIES.size())
        return gc.SPECIES[soot_index-1]->AGGLOMERATING;
    return false;
}
static int find_agglomeration_spec_index(const GlobalConstants& gc, int soot_index) {
    for (int i=0; i<gc.N_AGGLOMERATION_SPECIES; ++i)
        if (gc.AGGLOMERATION_SPEC_INDEX[i] == soot_index) return i+1;
    return 0;
}
// =============================================================================
// 桥接函数: 从 C 指针结构体构建 C++ 对象后转发
// =============================================================================
namespace {

void bridge_build_mp(MeshPointers& mp, std::vector<Array3D<double>*>& m3,
                     std::vector<Array4D<double>*>& m4, const SootMeshPointersC* c)
{
    int ib=c->IBAR, jb=c->JBAR, kb=c->KBAR, nt=c->N_TRACKED_SPECIES;
    auto w3=[&](double* p, int ni, int nj, int nk)->Array3D<double>*{
        if(!p)return nullptr; auto*a=new Array3D<double>(); a->wrap(p,ni,nj,nk); m3.push_back(a); return a;};
    auto w4=[&](double* p, int ni, int nj, int nk, int nl)->Array4D<double>*{
        if(!p)return nullptr; auto*a=new Array4D<double>(); a->wrap(p,ni,nj,nk,nl); m4.push_back(a); return a;};


    mp.U  = w3(c->U,  ib+1, jb+1, kb+1);
    mp.V  = w3(c->V,  ib+1, jb+1, kb+1);
    mp.W  = w3(c->W,  ib+1, jb+1, kb+1);
    mp.RHO    = w3(c->RHO,    ib+1, jb+1, kb+1);
    mp.RHOS   = w3(c->RHOS,   ib+1, jb+1, kb+1);
    mp.TMP    = w3(c->TMP,    ib+1, jb+1, kb+1);
    mp.MU     = w3(c->MU,     ib+1, jb+1, kb+1);
    mp.Q      = w3(c->Q,      ib+1, jb+1, kb+1);
    mp.D_SOURCE = w3(c->D_SOURCE, ib+1, jb+1, kb+1);
    mp.RSUM   = w3(c->RSUM,   ib+1, jb+1, kb+1);
    mp.WORK7  = w3(c->WORK7,  ib+1, jb+1, kb+1);
    mp.WORK8  = w3(c->WORK8,  ib+1, jb+1, kb+1);
    mp.WORK9  = w3(c->WORK9,  ib+1, jb+1, kb+1);

    mp.ZZ           = w4(c->ZZ,           ib+1, jb+1, kb+1, nt);
    mp.ZZS          = w4(c->ZZS,          ib+1, jb+1, kb+1, nt);
    mp.DEL_RHO_D_DEL_Z = w4(c->DEL_RHO_D_DEL_Z, ib+1, jb+1, kb+1, nt);
    mp.FX           = w4(c->FX,           ib+1, jb+1, kb+1, nt);
    mp.FY           = w4(c->FY,           ib+1, jb+1, kb+1, nt);
    mp.FZ           = w4(c->FZ,           ib+1, jb+1, kb+1, nt);
    mp.M_DOT_PPP    = w4(c->M_DOT_PPP,    ib+1, jb+1, kb+1, nt);

    static auto vecwrap = [](double* p, int n){ auto*v=new std::vector<double>(p,p+n); return v; };
    mp.R    = vecwrap(c->R,    ib+1);
    mp.RC   = vecwrap(c->RC,   ib+1);
    mp.RRN  = vecwrap(c->RRN,  ib+1);
    mp.DX   = vecwrap(c->DX,   ib+1);
    mp.RDX  = vecwrap(c->RDX,  ib+1);
    mp.RDXN = vecwrap(c->RDXN, ib+1);
    mp.DY   = vecwrap(c->DY,   jb+1);
    mp.RDY  = vecwrap(c->RDY,  jb+1);
    mp.DYN  = vecwrap(c->DYN,  jb+1);
    mp.RDYN = vecwrap(c->RDYN, jb+1);
    mp.DZ   = vecwrap(c->DZ,   kb+1);
    mp.RDZ  = vecwrap(c->RDZ,  kb+1);
    mp.DZN  = vecwrap(c->DZN,  kb+1);
    mp.RDZN = vecwrap(c->RDZN, kb+1);


    int pbar_nk = kb + 1;  // KBAR+1
    mp.PBAR_RAW = c->PBAR; // direct access to Fortran 2D PBAR
    mp.PBAR_NK  = pbar_nk;

    if (c->PBAR) mp.PBAR = new std::vector<double>(c->PBAR, c->PBAR + pbar_nk);
    else         mp.PBAR = nullptr;


    if (c->PRESSURE_ZONE) {
        auto* pz = new Array3D<int>();
        pz->wrap(c->PRESSURE_ZONE, ib+1, jb+1, kb+1);
        mp.PRESSURE_ZONE = pz;
    } else {
        mp.PRESSURE_ZONE = nullptr;
    }


    static int ibar_s=c->IBAR, jbar_s=c->JBAR, kbar_s=c->KBAR;
    static int n_ext_s=c->N_EXTERNAL_WALL_CELLS, n_int_s=c->N_INTERNAL_WALL_CELLS;
    static int icf_lo_s=c->INTERNAL_CFACE_CELLS_LB, n_icf_s=c->N_INTERNAL_CFACE_CELLS;
    mp.IBAR = &ibar_s; mp.JBAR = &jbar_s; mp.KBAR = &kbar_s;
    mp.N_EXTERNAL_WALL_CELLS = &n_ext_s; mp.N_INTERNAL_WALL_CELLS = &n_int_s;
    mp.INTERNAL_CFACE_CELLS_LB = &icf_lo_s; mp.N_INTERNAL_CFACE_CELLS = &n_icf_s;


    mp.NU_SOOT_OX = c->NU_SOOT_OX;
    mp.N_SURFACE_DENSITY_SPECIES = c->N_SURFACE_DENSITY_SPECIES;
    mp.I_MAX_TEMP_VAL = c->I_MAX_TEMP_VAL;
    mp.N_TOTAL_SCALARS_VAL = c->N_TOTAL_SCALARS_VAL;
    mp.D_Z_PTR = c->D_Z_PTR;


    int nw = c->N_EXTERNAL_WALL_CELLS + c->N_INTERNAL_WALL_CELLS;
    int nsds = c->N_SURFACE_DENSITY_SPECIES > 0 ? c->N_SURFACE_DENSITY_SPECIES : 1;
    auto* wvec = new std::vector<WallType>(nw);
    auto* bcvec = new std::vector<BoundaryCoordType>(nw);
    auto* bp1vec = new std::vector<BoundaryProp1Type>(nw);
    auto* bp2vec = new std::vector<BoundaryProp2Type>(nw);
    for (int i=0; i<nw; ++i) {
        (*wvec)[i].BOUNDARY_TYPE = c->WALL_BC_TYPE[i];
        (*wvec)[i].BC_INDEX = c->WALL_BC_INDEX[i];
        (*wvec)[i].B1_INDEX = c->WALL_B1_INDEX[i];
        (*bcvec)[i].IIG = c->BC_IIG[i]; (*bcvec)[i].JJG = c->BC_JJG[i];
        (*bcvec)[i].KKG = c->BC_KKG[i]; (*bcvec)[i].IOR = c->BC_IOR[i];
        (*bp1vec)[i].AREA = c->BP1_AREA[i];
        (*bp1vec)[i].TMP_F = c->BP1_TMP_F[i];
        (*bp1vec)[i].RHO_G = c->BP1_RHO_G[i];
        (*bp1vec)[i].TMP_G = c->BP1_TMP_G[i];
        (*bp1vec)[i].NSDS = nsds;

        (*bp1vec)[i].AWM_AEROSOL = c->BP1_AWM_AEROSOL ? &c->BP1_AWM_AEROSOL[i * nsds] : nullptr;

        (*bp1vec)[i].ZZ_G = c->BP1_ZZ_G ? &c->BP1_ZZ_G[i * nt] : nullptr;
    }
    mp.WALL = wvec; mp.BOUNDARY_COORD = bcvec;
    mp.BOUNDARY_PROP1 = bp1vec; mp.BOUNDARY_PROP2 = bp2vec;


    int nc = (ib+1)*(jb+1)*(kb+1);
    auto* ci = new Array3D<int>();
    ci->wrap(c->CELL_INDEX, ib+1, jb+1, kb+1);
    mp.CELL_INDEX = ci;
    auto* cv = new std::vector<CellType>(nc);
    if (c->CELL_SOLID) for (int i=0; i<nc; ++i) (*cv)[i].SOLID = (c->CELL_SOLID[i] != 0.0);
    mp.CELL = cv;

    int nicf = c->N_INTERNAL_CFACE_CELLS;
    auto* cfvec = new std::vector<CfaceType>(nicf);
    for (int i=0; i<nicf; ++i) {
        (*cfvec)[i].BOUNDARY_TYPE = c->CFACE_BC_TYPE[i];
        (*cfvec)[i].B1_INDEX = c->CFACE_B1_INDEX[i];
        (*cfvec)[i].BC_INDEX = c->CFACE_BC_INDEX[i];
        (*cfvec)[i].RSUM_G = c->CFACE_RSUM_G[i];
    }
    mp.CFACE = cfvec;

    mp.PREDICTOR = (c->PREDICTOR_FLAG != 0);
}

void bridge_build_gc(GlobalConstants& gc, const SootMeshPointersC* c) {
    gc.N_TRACKED_SPECIES = c->N_TRACKED_SPECIES;
    gc.GVEC = {c->GVEC_X, c->GVEC_Y, c->GVEC_Z};
    gc.NULL_BOUNDARY = c->NULL_BOUNDARY_VAL;
    gc.OPEN_BOUNDARY = c->OPEN_BOUNDARY_VAL;
    gc.INTERPOLATED_BOUNDARY = c->INTERPOLATED_BOUNDARY_VAL;
    gc.GRAVITATIONAL_SETTLING = true;
    gc.THERMOPHORETIC_SETTLING = true;
    gc.N_AGGLOMERATION_SPECIES = c->N_AGGLOMERATION_SPECIES;
    gc.SOOT_INDEX = c->SOOT_INDEX;
    gc.O2_INDEX = c->O2_INDEX;
    gc.MW_O2_VAL = c->MW_O2_VAL;
    gc.ZZ_MIN_GLOBAL_VAL = c->ZZ_MIN_GLOBAL_VAL;

    int na = c->N_AGGLOMERATION_SPECIES;
    gc.N_PARTICLE_BINS.assign(c->N_PARTICLE_BINS, c->N_PARTICLE_BINS + na);
    gc.MIN_PARTICLE_DIAMETER.assign(c->MIN_PARTICLE_DIAMETER, c->MIN_PARTICLE_DIAMETER + na);
    gc.MAX_PARTICLE_DIAMETER.assign(c->MAX_PARTICLE_DIAMETER, c->MAX_PARTICLE_DIAMETER + na);
    gc.AGGLOMERATION_SPEC_INDEX.assign(c->AGGLOMERATION_SPEC_INDEX, c->AGGLOMERATION_SPEC_INDEX + na);
    gc.AGGLOMERATION_SMIX_INDEX.assign(c->AGGLOMERATION_SMIX_INDEX, c->AGGLOMERATION_SMIX_INDEX + na);

    int nts = c->N_TRACKED_SPECIES;
    gc.SPECIES_MIXTURE.clear();
    for (int i=0; i<nts; ++i) {
        auto* sm = new SpeciesMixtureType();
        sm->DEPOSITING = c->SM_DEPOSITING[i] != 0.0;
        sm->MEAN_DIAMETER = c->SM_MEAN_DIAMETER[i];
        sm->DENSITY_SOLID = c->SM_DENSITY_SOLID[i];
        sm->THERMOPHORETIC_DIAMETER = c->SM_THERMOPHORETIC_DIAMETER[i];
        sm->CONDUCTIVITY_SOLID = c->SM_CONDUCTIVITY_SOLID[i];
        sm->H_F = c->SM_H_F[i]; sm->RCON = c->SM_RCON[i]; sm->MW = c->SM_MW[i];
        sm->AWM_INDEX = c->SM_AWM_INDEX[i];
        sm->AGGLOMERATION_INDEX = c->SM_AGGLOMERATION_INDEX[i];
        gc.SPECIES_MIXTURE.push_back(sm);
    }

    // NU_SOOT_OX
    gc.NU_SOOT_OX.assign(c->NU_SOOT_OX, c->NU_SOOT_OX + nts);

    int nsp = c->N_SPECIES;
    gc.SPECIES.clear();
    for (int i=0; i<nsp; ++i) {
        auto* sp = new SpecType();
        sp->DENSITY_SOLID = c->SP_DENSITY_SOLID[i];
        sp->CONDUCTIVITY_SOLID = c->SP_CONDUCTIVITY_SOLID[i];
        sp->AWM_INDEX = c->SP_AWM_INDEX[i];
        sp->AGGLOMERATING = c->SP_AGGLOMERATING[i] != 0;
        gc.SPECIES.push_back(sp);
    }
}

void cleanup_mp_allocations(std::vector<Array3D<double>*>& m3,
                            std::vector<Array4D<double>*>& m4,
                            MeshPointers& mp) {
    for (auto* p : m3) delete p;
    for (auto* p : m4) delete p;
    m3.clear(); m4.clear();
    if (mp.WALL)           { delete mp.WALL;           mp.WALL = nullptr; }
    if (mp.BOUNDARY_COORD) { delete mp.BOUNDARY_COORD; mp.BOUNDARY_COORD = nullptr; }
    if (mp.BOUNDARY_PROP1) { delete mp.BOUNDARY_PROP1; mp.BOUNDARY_PROP1 = nullptr; }
    if (mp.BOUNDARY_PROP2) { delete mp.BOUNDARY_PROP2; mp.BOUNDARY_PROP2 = nullptr; }
    if (mp.CELL)           { delete mp.CELL;           mp.CELL = nullptr; }
    if (mp.CELL_INDEX)     { delete mp.CELL_INDEX;     mp.CELL_INDEX = nullptr; }
    if (mp.CFACE)          { delete mp.CFACE;          mp.CFACE = nullptr; }
    if (mp.PRESSURE_ZONE)  { delete mp.PRESSURE_ZONE;  mp.PRESSURE_ZONE = nullptr; }
    if (mp.PBAR)           { delete mp.PBAR;           mp.PBAR = nullptr; }
}

} // anonymous namespace

// =============================================================================
// extern "C" 桥接函数实现
// =============================================================================

void soot_c_initialize_agglomeration(const SootMeshPointersC* c) {
    GlobalConstants gc; bridge_build_gc(gc, c);
    soot_routines::initialize_agglomeration(gc);
}

void soot_c_settling_velocity(const SootMeshPointersC* c) {
    MeshPointers mp;
    std::vector<Array3D<double>*> m3; std::vector<Array4D<double>*> m4;
    bridge_build_mp(mp, m3, m4, c);
    GlobalConstants gc; bridge_build_gc(gc, c);
    soot_routines::settling_velocity(0, mp, gc);
    cleanup_mp_allocations(m3, m4, mp);
}

void soot_c_calc_agglomeration(double DT, const SootMeshPointersC* c) {
    MeshPointers mp;
    std::vector<Array3D<double>*> m3; std::vector<Array4D<double>*> m4;
    bridge_build_mp(mp, m3, m4, c);
    GlobalConstants gc; bridge_build_gc(gc, c);
    soot_routines::calc_agglomeration(DT, 0, mp, gc);
    cleanup_mp_allocations(m3, m4, mp);
}

void soot_c_soot_surface_oxidation(double DT, const SootMeshPointersC* c) {
    MeshPointers mp;
    std::vector<Array3D<double>*> m3; std::vector<Array4D<double>*> m4;
    bridge_build_mp(mp, m3, m4, c);
    GlobalConstants gc; bridge_build_gc(gc, c);
    soot_routines::soot_surface_oxidation(DT, 0, mp, gc);
    cleanup_mp_allocations(m3, m4, mp);
}

void soot_c_droplet_scrubbing(int IP, double DT, double DT_P, const SootMeshPointersC* c) {
    MeshPointers mp;
    std::vector<Array3D<double>*> m3; std::vector<Array4D<double>*> m4;
    bridge_build_mp(mp, m3, m4, c);
    GlobalConstants gc; bridge_build_gc(gc, c);
    // Override Lagrangian particle data from scalar values
    LagrangianParticleType lp;
    lp.U = c->LP_U_VAL; lp.V = c->LP_V_VAL; lp.W = c->LP_W_VAL;
    lp.RADIUS = c->LP_RADIUS_VAL; lp.RVC = c->LP_RVC_VAL; lp.PWT = c->LP_PWT_VAL;
    lp.BC_INDEX = c->LP_BC_INDEX_VAL; lp.B1_INDEX = c->LP_B1_INDEX_VAL;
    // Inject into boundary structures for the C++ routine to find
    if (c->LP_BC_INDEX_VAL > 0 && c->LP_BC_INDEX_VAL <= (int)mp.BOUNDARY_COORD->size()) {
        auto& bc = (*mp.BOUNDARY_COORD)[c->LP_BC_INDEX_VAL-1];
        bc.IIG = c->LP_IIG_VAL; bc.JJG = c->LP_JJG_VAL; bc.KKG = c->LP_KKG_VAL;
        // Override first boundary prop with LP data
        if (!mp.BOUNDARY_PROP1->empty()) {
            auto& b1 = (*mp.BOUNDARY_PROP1)[0];
            b1.AREA = c->LP_BP1_AREA_VAL;
            b1.TMP_F = c->LP_BP1_TMP_F_VAL;
            b1.RHO_G = c->LP_BP1_RHO_G_VAL;
            b1.TMP_G = c->LP_BP1_TMP_G_VAL;
        }
    }
    soot_routines::droplet_scrubbing(IP, 0, DT, DT_P, mp, gc);
    cleanup_mp_allocations(m3, m4, mp);
}

int soot_c_n_agglomeration_species(void) {
    return (int)soot_routines::BIN_S.size();
}

int soot_c_max_particle_bins(void) {
    int mx = 0;
    for (const auto& v : soot_routines::BIN_X)
        if ((int)v.size() > mx) mx = (int)v.size();
    return mx;
}

double soot_c_particle_radius(int sp, int bin) {
    if (sp < 0 || sp >= (int)soot_routines::PARTICLE_RADIUS.size()) return 0.0;
    if (bin < 0 || bin >= (int)soot_routines::PARTICLE_RADIUS[sp].size()) return 0.0;
    return soot_routines::PARTICLE_RADIUS[sp][bin];
}
