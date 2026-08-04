#ifndef MASS_ROUTINES_H
#define MASS_ROUTINES_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// =============================================================================
// mass 模块 C++ 转换 — 类型定义和函数声明

// =============================================================================


constexpr int NULL_BOUNDARY        = 0;
constexpr int SOLID_BOUNDARY       = 1;
constexpr int OPEN_BOUNDARY        = 2;
constexpr int MIRROR_BOUNDARY      = 3;
constexpr int INTERPOLATED_BOUNDARY= 6;
constexpr int PERIODIC_BOUNDARY    = 7;


constexpr double ONTH           = 1.0 / 3.0;
constexpr double FOTH           = 4.0 / 3.0;
constexpr double TWTH           = 2.0 / 3.0;
constexpr double PI             = 3.14159265358979323846;
constexpr double FOTHPI         = (4.0 / 3.0) * PI;
constexpr double TWO_EPSILON_EB = 2.0 * 2.2204460492503131e-16;
constexpr double TINY_EB        = 2.2250738585072014e-308;


struct SpecType {
    double DENSITY_SOLID{1800.0};
    double CONDUCTIVITY_SOLID{0.1};
    int AWM_INDEX{0};
    bool AGGLOMERATING{false};
};

struct SpeciesMixtureType {
    bool DEPOSITING{false};
    double MEAN_DIAMETER{1e-6};
    double DENSITY_SOLID{1800.0};
    double THERMOPHORETIC_DIAMETER{1e-6};
    double CONDUCTIVITY_SOLID{0.1};
    double H_F{0.0};
    double RCON{287.0};
    double MW{29.0};
    int AWM_INDEX{0};
    int AGGLOMERATION_INDEX{0};
};

struct WallType {
    int BOUNDARY_TYPE{0};
    int BC_INDEX{0};
    int B1_INDEX{0};
};

struct BoundaryCoordType {
    int II{0}, JJ{0}, KK{0};
    int IIG{0}, JJG{0}, KKG{0};
    int IOR{0};
};

struct BoundaryProp1Type {
    double AREA{0.0};
    double TMP_F{0.0}, RHO_F{0.0};
    double TMP_G{0.0}, RHO_G{0.0};
    double* ZZ_F{nullptr};
};

struct ExternalWallType {
    int BOUNDARY_TYPE_PREVIOUS{0};
};

struct CellType {
    bool SOLID{false};
    bool EXTERIOR{false};
    int WALL_INDEX[7] = {0}; // -3..+3, 索引从 0 开始
    inline int wi(int dir) const { return WALL_INDEX[dir+3]; }
};

typedef struct {

    void* U, *V, *W;
    void* US, *VS, *WS;

    void* RHO, *RHOS;
    void* TMP, *RSUM, *MU, *Q, *D_SOURCE;

    void* ZZ, *ZZS;              // (0:IBP1, 0:JBP1, 0:KBP1, N_TRACKED_SPECIES)
    void* DEL_RHO_D_DEL_Z;       // (0:IBP1, 0:JBP1, 0:KBP1, N_TRACKED_SPECIES)
    void* FX, *FY, *FZ;          // (0:IBAR, 0:JBAR, 0:KBAR, 0:N_TOTAL_SCALARS) 含索引 0
    void* M_DOT_PPP;             // (0:IBP1, 0:JBP1, 0:KBP1, N_TRACKED_SPECIES)
    void* ADV_FX, *ADV_FY, *ADV_FZ; // (0:IBP1, 0:JBP1, 0:KBP1, N_TRACKED_SPECIES)

    void* R, *RC, *RRN;
    void* DX, *RDX, *RDXN;
    void* DY, *RDY, *DYN, *RDYN;
    void* DZ, *RDZ, *DZN, *RDZN;
    void* XC, *ZC;

    void* WORK1, *WORK2, *WORK3, *WORK4, *WORK5;
    void* SWORK4;  
    void* PBAR;           
    void* PBAR_S;         
    void* PRESSURE_ZONE;  
    void* D_PBAR_DT;     
    void* D_PBAR_DT_S;   
  
    void* CELL_INDEX;    
    void* CELL_SOLID;     
    void* CELL_WALL_INDEX; 
    void* CELL_EXTERIOR;  
    void* UVW_SAVE;       
    
    void* WALL_BC_TYPE;   // int[]
    void* WALL_B1_INDEX;  // int[]
    void* BC_IIG, *BC_JJG, *BC_KKG, *BC_IOR; // int[]
    void* BC_II, *BC_JJ, *BC_KK;            
    void* BP1_RHO_F;      // double[]
    void* BP1_ZZ_F;      
    // --- EXTERNAL_WALL ---
    void* EWC_BOUNDARY_TYPE_PREVIOUS; // int[]
    
    int IBAR, JBAR, KBAR;
    int IBP1, JBP1, KBP1;  // IBAR+1, JBAR+1, KBAR+1
    int N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS;
    int N_ZONE;
    int N_TOTAL_SCALARS, N_TRACKED_SPECIES;
    int N_PASSIVE_SCALARS, ZETA_INDEX;

    int PREDICTOR_FLAG;
    int FIRST_PASS_FLAG;
    int CC_IBM_FLAG;
    int STORE_SPECIES_FLUX_FLAG;
    int I_FLUX_LIMITER;
    int FLUX_LIMITER_MW_CORRECTION_FLAG;
    int PERIODIC_TEST_VAL;
    int ICYC_VAL;
    int SOLID_PHASE_ONLY_FLAG;
    int GRAVITATIONAL_SETTLING_FLAG;
    int THERMOPHORETIC_SETTLING_FLAG;
    double RHOMIN_VAL, RHOMAX_VAL;

    void* SM_MW;           // double[]
    void* SM_DEPOSITING;   // double[] (0.0/1.0 as flag)
    // --- SPECIES (长度 N_SPECIES) ---
    void* SP_DENSITY_SOLID; // double[]
} MassMeshPointersC;


struct IX3 {
    int ni, nj, nk;
    IX3() : ni(1), nj(1), nk(1) {}
    IX3(int ni_, int nj_, int nk_) : ni(ni_), nj(nj_), nk(nk_) {}
    inline int operator()(int i, int j, int k) const {
        return i + ni * (j + nj * k);
    }
};

struct IX4 {
    int ni, nj, nk, nl, l0;
    IX4() : ni(1), nj(1), nk(1), nl(1), l0(0) {}
    IX4(int ni_, int nj_, int nk_, int nl_) : ni(ni_), nj(nj_), nk(nk_), nl(nl_), l0(0) {}
    IX4(int ni_, int nj_, int nk_, int nl_, int l0_) : ni(ni_), nj(nj_), nk(nk_), nl(nl_), l0(l0_) {}
    inline int operator()(int i, int j, int k, int l) const {
        return i + ni * (j + nj * (k + nk * (l - l0)));
    }
};

// =============================================================================
// extern "C" 桥接函数 — 从 Fortran 包装器调用
// =============================================================================
#ifdef __cplusplus
extern "C" {
#endif


    void mass_c_finite_differences(const MassMeshPointersC* c);
    void mass_c_density(double T, double DT, const MassMeshPointersC* c);


    void fds_get_scalar_face_value(const double* uu, const double* scalar,
                                   double* face,
                                   int ilo, int ihi, int jlo, int jhi,
                                   int klo, int khi, int dir,
                                   int i_flux_limiter,
                                   int ni, int nj, int nk);
    void fds_get_specific_gas_constant(const double* zz_get,
                                       double* rsum_out);
    void fds_get_molecular_weight(const double* zz_get,
                                  double* mw_out);

#ifdef __cplusplus
}
#endif

#endif // MASS_ROUTINES_H
