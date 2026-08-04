#ifndef SOOT_ROUTINES_H
#define SOOT_ROUTINES_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

constexpr double ONTH           = 1.0 / 3.0;
constexpr double FOTH           = 4.0 / 3.0;
constexpr double TWTH           = 2.0 / 3.0;
constexpr double PI             = 3.14159265358979323846;
constexpr double FOTHPI         = (4.0 / 3.0) * PI;
constexpr double TWO_EPSILON_EB = 2.0 * 2.2204460492503131e-16;
constexpr double TINY_EB        = 2.2250738585072014e-308;
constexpr double K_BOLTZMANN    = 1.3806488e-23;
constexpr double GRAV           = 9.80665;

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

struct CellType {
    bool SOLID{false};
};

struct WallType {
    int BOUNDARY_TYPE{0};
    int BC_INDEX{1};
    int B1_INDEX{1};
};

struct BoundaryCoordType {
    int IIG{1}, JJG{1}, KKG{1}, IOR{0};
};

struct BoundaryProp1Type {
    double AREA{1.0};
    double TMP_F{300.0};
    double RHO_G{1.2};
    double TMP_G{300.0};
    double* AWM_AEROSOL{nullptr};  
    double* ZZ_G{nullptr};          
    int NSDS{0};                   
};

struct BoundaryProp2Type {
    double V_DEP{0.0};
};

struct CfaceType {
    int BOUNDARY_TYPE{0};
    int B1_INDEX{1};
    int BC_INDEX{1};
    double RSUM_G{287.0};
};

struct LagrangianParticleType {
    double U{0.0}, V{0.0}, W{0.0};
    double RADIUS{1e-6};
    double RVC{1.0};
    double PWT{1.0};
    int BC_INDEX{1};
    int B1_INDEX{1};
};


template<typename T>
class Array3D {
public:
    int ni{0}, nj{0}, nk{0};
    std::vector<T> data;
    T* ext_ptr{nullptr};  // non-owning external memory (Fortran)

    Array3D() = default;
    Array3D(int i, int j, int k) : ni(i), nj(j), nk(k), data(i*j*k, T{}) {}

    void resize(int i, int j, int k) { ni = i; nj = j; nk = k; data.assign(i*j*k, T{}); ext_ptr = nullptr; }
    void assign(T val) {
        if (ext_ptr) { std::fill_n(ext_ptr, ni*nj*nk, val); }
        else         { std::fill(data.begin(), data.end(), val); }
    }

    void wrap(T* ext, int i, int j, int k) { ext_ptr = ext; ni = i; nj = j; nk = k; data.clear(); }

    T* raw() { return ext_ptr ? ext_ptr : data.data(); }

          T& operator()(int i, int j, int k)       { return idx(i,j,k); }
    const T& operator()(int i, int j, int k) const { return const_cast<Array3D*>(this)->idx(i,j,k); }

private:
    T& idx(int i, int j, int k) {
        T* base = ext_ptr ? ext_ptr : data.data();
        return base[((k-1)*nj + (j-1))*ni + (i-1)];
    }
};

template<typename T>
class Array4D {
public:
    int ni{0}, nj{0}, nk{0}, nl{0};
    std::vector<T> data;
    T* ext_ptr{nullptr};  // non-owning external memory (Fortran)

    Array4D() = default;
    Array4D(int i, int j, int k, int l) : ni(i), nj(j), nk(k), nl(l), data(i*j*k*l, T{}) {}

    void resize(int i, int j, int k, int l) { ni = i; nj = j; nk = k; nl = l; data.assign(i*j*k*l, T{}); ext_ptr = nullptr; }
    void assign(T val) {
        if (ext_ptr) { std::fill_n(ext_ptr, ni*nj*nk*nl, val); }
        else         { std::fill(data.begin(), data.end(), val); }
    }

    void wrap(T* ext, int i, int j, int k, int l) { ext_ptr = ext; ni = i; nj = j; nk = k; nl = l; data.clear(); }

          T& operator()(int i, int j, int k, int l)       { return idx(i,j,k,l); }
    const T& operator()(int i, int j, int k, int l) const { return const_cast<Array4D*>(this)->idx(i,j,k,l); }

private:
    T& idx(int i, int j, int k, int l) {
        T* base = ext_ptr ? ext_ptr : data.data();
        return base[(((l-1)*nk + (k-1))*nj + (j-1))*ni + (i-1)];
    }
};

struct MeshPointers {
    Array3D<double>  *U{nullptr}, *V{nullptr}, *W{nullptr};
    Array3D<double>  *RHO{nullptr}, *RHOS{nullptr};
    Array3D<double>  *TMP{nullptr}, *MU{nullptr};
    Array3D<double>  *Q{nullptr}, *D_SOURCE{nullptr}, *RSUM{nullptr};
    Array4D<double>  *ZZ{nullptr}, *ZZS{nullptr};
    Array4D<double>  *DEL_RHO_D_DEL_Z{nullptr};
    Array4D<double>  *FX{nullptr}, *FY{nullptr}, *FZ{nullptr};
    Array4D<double>  *M_DOT_PPP{nullptr};

    std::vector<double> *R{nullptr}, *RC{nullptr}, *RRN{nullptr};
    std::vector<double> *DX{nullptr}, *RDX{nullptr}, *RDXN{nullptr};
    std::vector<double> *DY{nullptr}, *RDY{nullptr}, *DYN{nullptr}, *RDYN{nullptr};
    std::vector<double> *DZ{nullptr}, *RDZ{nullptr}, *DZN{nullptr}, *RDZN{nullptr};
    double           *PBAR_RAW{nullptr};   // Fortran 2D PBAR(0:KBAR,1:N_ZONE), column-major
    int               PBAR_NK{0};         
    std::vector<double> *PBAR{nullptr};  

    double *NU_SOOT_OX{nullptr};
    int    N_SURFACE_DENSITY_SPECIES{0};
    int    I_MAX_TEMP_VAL{0}, N_TOTAL_SCALARS_VAL{0};
    double *D_Z_PTR{nullptr};
    Array3D<int>     *PRESSURE_ZONE{nullptr};

    int *IBAR{nullptr}, *JBAR{nullptr}, *KBAR{nullptr};
    int *N_EXTERNAL_WALL_CELLS{nullptr}, *N_INTERNAL_WALL_CELLS{nullptr};
    int *INTERNAL_CFACE_CELLS_LB{nullptr}, *N_INTERNAL_CFACE_CELLS{nullptr};

    std::vector<WallType>              *WALL{nullptr};
    std::vector<CfaceType>             *CFACE{nullptr};
    std::vector<BoundaryCoordType>     *BOUNDARY_COORD{nullptr};
    std::vector<BoundaryProp1Type>     *BOUNDARY_PROP1{nullptr};
    std::vector<BoundaryProp2Type>     *BOUNDARY_PROP2{nullptr};
    std::vector<CellType>              *CELL{nullptr};
    Array3D<int>     *CELL_INDEX{nullptr};
    Array3D<double>  *WORK7{nullptr}, *WORK8{nullptr}, *WORK9{nullptr};

    bool PREDICTOR{true};
};

// ---------------------------------------------------------------------------
// Fortran 全局常量上下文 — 传递给 C++ 函数的参数集合
// ---------------------------------------------------------------------------
struct GlobalConstants {
    int N_TRACKED_SPECIES{0};
    int SOOT_INDEX{1};
    int O2_INDEX{2};
    double MW_O2_VAL{32.0};
    double ZZ_MIN_GLOBAL_VAL{1e-10};
    std::vector<double> NU_SOOT_OX;
    std::vector<double> GVEC{0.0, 0.0, -9.80665};
    int NULL_BOUNDARY{0};
    int OPEN_BOUNDARY{2};
    int INTERPOLATED_BOUNDARY{6};
    bool GRAVITATIONAL_SETTLING{true};
    bool THERMOPHORETIC_SETTLING{true};

    int N_AGGLOMERATION_SPECIES{0};
    std::vector<int> N_PARTICLE_BINS;
    std::vector<double> MIN_PARTICLE_DIAMETER;
    std::vector<double> MAX_PARTICLE_DIAMETER;
    std::vector<int> AGGLOMERATION_SPEC_INDEX;
    std::vector<int> AGGLOMERATION_SMIX_INDEX;
    std::vector<SpecType*> SPECIES;
    std::vector<SpeciesMixtureType*> SPECIES_MIXTURE;
};

// =============================================================================
// soot_routines 命名空间 — 模块级全局变量和子程序
// =============================================================================
namespace soot_routines {

extern double MIN_AGGLOMERATION;
extern std::vector<double> BIN_S;
extern std::vector< std::vector<double> > BIN_M;
extern std::vector< std::vector<double> > BIN_X;
extern std::vector< std::vector<double> > MOBILITY_FAC;
extern std::vector< std::vector<double> > A_FAC;
extern std::vector< std::vector<double> > PARTICLE_RADIUS;
extern std::vector< std::vector< std::vector<double> > > PHI_B_FAC;
extern std::vector< std::vector< std::vector<double> > > PHI_G_FAC;
extern std::vector< std::vector< std::vector<double> > > PHI_S_FAC;
extern std::vector< std::vector< std::vector<double> > > PHI_I_FAC;
extern std::vector< std::vector< std::vector<double> > > FU1_FAC;
extern std::vector< std::vector< std::vector<double> > > FU2_FAC;
extern std::vector< std::vector< std::vector<double> > > PARTICLE_MASS;
extern std::vector< std::vector< std::vector< std::vector<double> > > > BIN_ETA;
extern std::vector< std::vector< std::vector< std::vector<int> > > > BIN_ETA_INDEX;

/// 重力沉降和热泳 (气体相) — Fortran SUBROUTINE SETTLING_VELOCITY(NM)
void settling_velocity(int nm, MeshPointers& mp, const GlobalConstants& gc);

/// 初始化粒子箱和团聚核预因子 — Fortran SUBROUTINE INITIALIZE_AGGLOMERATION
void initialize_agglomeration(const GlobalConstants& gc);

/// 团聚主计算 (含自适应子步进) — Fortran SUBROUTINE CALC_AGGLOMERATION(DT,NM)
void calc_agglomeration(double DT, int nm, MeshPointers& mp, const GlobalConstants& gc);

/// 壁面积碳氧化 — Fortran SUBROUTINE SOOT_SURFACE_OXIDATION(DT,NM)
void soot_surface_oxidation(double DT, int nm, MeshPointers& mp, const GlobalConstants& gc);

/// 液滴洗涤 — Fortran SUBROUTINE DROPLET_SCRUBBING(IP,NM,DT,DT_P)
void droplet_scrubbing(int IP, int nm, double DT, double DT_P, MeshPointers& mp, const GlobalConstants& gc);

} // namespace soot_routines

// =============================================================================
// C 兼容桥接结构体 + extern "C" 函数 — 供 Fortran -> C++ 混合编译
// =============================================================================
#ifdef __cplusplus
extern "C" {
#endif


typedef struct {

    double *U, *V, *W;
    double *RHO, *RHOS;
    double *TMP, *MU;
    double *Q, *D_SOURCE, *RSUM;

    double *ZZ, *ZZS;
    double *DEL_RHO_D_DEL_Z;
    double *FX, *FY, *FZ;
    double *M_DOT_PPP;

    double *R, *RC, *RRN;
    double *DX, *RDX, *RDXN;
    double *DY, *RDY, *DYN, *RDYN;
    double *DZ, *RDZ, *DZN, *RDZN;
    double *PBAR;
    int    *PRESSURE_ZONE;
    int    *CELL_INDEX;
    int     IBAR, JBAR, KBAR;
    int     N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS;
    int     INTERNAL_CFACE_CELLS_LB, N_INTERNAL_CFACE_CELLS;
    int     N_SPECIES;

    int    *WALL_BC_TYPE;
    int    *WALL_BC_INDEX;
    int    *WALL_B1_INDEX;
    int    *BC_IIG, *BC_JJG, *BC_KKG, *BC_IOR;
    double *BP1_AREA, *BP1_TMP_F;
    double *BP1_RHO_G, *BP1_TMP_G;
    double *BP1_AWM_AEROSOL;
    double *BP1_ZZ_G;

    int    *CFACE_BC_TYPE;
    int    *CFACE_B1_INDEX;
    int    *CFACE_BC_INDEX;
    double *CFACE_RSUM_G;

    double *LP_U, *LP_V, *LP_W;
    double *LP_RADIUS, *LP_RVC, *LP_PWT;
    int    *LP_BC_INDEX, *LP_B1_INDEX;

    double  LP_U_VAL, LP_V_VAL, LP_W_VAL;
    double  LP_RADIUS_VAL, LP_RVC_VAL, LP_PWT_VAL;
    int     LP_BC_INDEX_VAL, LP_B1_INDEX_VAL;

    int     LP_IIG_VAL, LP_JJG_VAL, LP_KKG_VAL;
    double  LP_BP1_AREA_VAL, LP_BP1_TMP_F_VAL;
    double  LP_BP1_RHO_G_VAL, LP_BP1_TMP_G_VAL;

    double *CELL_SOLID;

    double *WORK7, *WORK8, *WORK9;

    int     PREDICTOR_FLAG;
    int     N_TRACKED_SPECIES;
    double  GVEC_X, GVEC_Y, GVEC_Z;
    int     NULL_BOUNDARY_VAL, OPEN_BOUNDARY_VAL, INTERPOLATED_BOUNDARY_VAL;
    int     N_AGGLOMERATION_SPECIES;
    int     SOOT_INDEX, O2_INDEX;
    int     N_SURFACE_DENSITY_SPECIES;
    double  R0_VAL, MW_O2_VAL, ZZ_MIN_GLOBAL_VAL;
    double *NU_SOOT_OX;            
    int     I_MAX_TEMP_VAL, N_TOTAL_SCALARS_VAL;
    double *D_Z_PTR;               
    
    double *SM_DEPOSITING;
    double *SM_MEAN_DIAMETER, *SM_DENSITY_SOLID;
    double *SM_THERMOPHORETIC_DIAMETER, *SM_CONDUCTIVITY_SOLID;
    double *SM_H_F, *SM_RCON, *SM_MW;
    int    *SM_AWM_INDEX, *SM_AGGLOMERATION_INDEX;
   
    double *SP_DENSITY_SOLID, *SP_CONDUCTIVITY_SOLID;
    int    *SP_AWM_INDEX, *SP_AGGLOMERATING;
   
    int    *N_PARTICLE_BINS;
    int    *AGGLOMERATION_SPEC_INDEX, *AGGLOMERATION_SMIX_INDEX;
    double *MIN_PARTICLE_DIAMETER, *MAX_PARTICLE_DIAMETER;
} SootMeshPointersC;

// ---- Fortran 回调函数: C++ → Fortran 调用物理函数 ----
void fds_get_viscosity(const double* zz_arr, double* mu_out, double tmpg);
void fds_get_conductivity(const double* zz_arr, double* k_out, double tmpg);
void fds_get_mass_fraction(const double* zz_arr, int idx, double* y_out);
void fds_get_molecular_weight(const double* zz_arr, double* mw_out);
void fds_get_specific_heat(const double* zz_arr, double* cp_out, double tmpg);
void fds_get_specific_gas_constant(const double* zz_arr, double* rsum_out);
void fds_get_average_specific_heat(const double* zz_arr, double* cpbar_out, double tmpg);
void fds_d_z_lookup(double tmp_val, int ns_val, double* d_out);

// ---- 桥接函数: Fortran 可调用 (extern "C", 无 C++ name mangling) ----
void soot_c_initialize_agglomeration(const SootMeshPointersC* mp);
void soot_c_settling_velocity(const SootMeshPointersC* mp);
void soot_c_calc_agglomeration(double DT, const SootMeshPointersC* mp);
void soot_c_soot_surface_oxidation(double DT, const SootMeshPointersC* mp);
void soot_c_droplet_scrubbing(int IP, double DT, double DT_P, const SootMeshPointersC* mp);

// ---- 模块变量访问函数 ----
int    soot_c_n_agglomeration_species(void);
int    soot_c_max_particle_bins(void);
double soot_c_particle_radius(int sp, int bin);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // SOOT_ROUTINES_H
