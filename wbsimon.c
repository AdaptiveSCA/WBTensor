#include "wbsimon.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// =============================================================
// Helper Functions
// =============================================================

/**
 * Flips a specific bit in a Boolean matrix.
 * Used for applying masking noise.
 */
void mzd_flip_bit(mzd_t *M, int row, int col) {
    int val = mzd_read_bit(M, row, col);
    mzd_write_bit(M, row, col, val ^ 1);
}

/**
 * Sets all bits in a matrix to zero.
 * Ensures a clean state for matrix generation.
 */
void wb_mzd_set_zero(mzd_t *M) {
    for(int r = 0; r < M->nrows; r++) {
        for(int c = 0; c < M->ncols; c++) {
            mzd_write_bit(M, r, c, 0);
        }
    }
}

/**
 * Maps 2D indices (i, j) to the 1D index of the flattened upper-triangular quadratic vector.
 * White-box rows store quadratic forms x_i*x_j as a single vector.
 * The mapping order is Row-Major Upper-Triangular: (0,0), (0,1)...(0,N-1), (1,1)...
 */
int get_quad_idx(int i, int j, int N) {
    // Ensure i <= j for upper triangular form
    if (i > j) { int tmp = i; i = j; j = tmp; }
    // Formula derived from arithmetic progression sum of row lengths
    return (i * N) - (i * (i + 1)) / 2 + j;
}

/**
 * Circular right shift for 16-bit integers.
 * Essential for the SIMON key schedule and round function.
 */
static inline uint16_t rotr16(uint16_t x, int n) {
    return (x >> n) | (x << (16 - n));
}

// =============================================================
// Linear Algebra Helpers for White-Box Encodings
// =============================================================

/**
 * Generates a random N x N invertible binary matrix and its inverse.
 * Uses Gaussian elimination on an augmented matrix [P | I] to ensure robustness.
 * * @param N Dimension of the matrix.
 * @param out_P Pointer to store the generated matrix P.
 * @param out_P_inv Pointer to store the calculated inverse P^-1.
 */
void wb_gen_invertible_matrix(int N, mzd_t **out_P, mzd_t **out_P_inv) {
    mzd_t *P = mzd_init(N, N);
    mzd_t *Aug = mzd_init(N, 2 * N);
    mzd_t *Identity = mzd_init(N, N);
    
    // Prepare Identity matrix
    wb_mzd_set_zero(Identity);
    for(int i=0; i<N; i++) mzd_write_bit(Identity, i, i, 1);

    while (1) {
        mzd_randomize(P);
        
        // Construct Augmented Matrix [P | Identity]
        wb_mzd_set_zero(Aug);
        // Copy P to the left half
        for(int r=0; r<N; r++) {
            for(int c=0; c<N; c++) {
                mzd_write_bit(Aug, r, c, mzd_read_bit(P, r, c));
            }
        }
        // Copy Identity to the right half
        for(int r=0; r<N; r++) {
            for(int c=0; c<N; c++) {
                mzd_write_bit(Aug, r, N + c, mzd_read_bit(Identity, r, c));
            }
        }

        // Perform Gaussian Elimination to reach Row Echelon Form
        mzd_echelonize(Aug, 1);

        // Check for singularity: The left half must become an Identity matrix.
        // A simple check is ensuring the main diagonal is all 1s.
        int singular = 0;
        for(int i=0; i<N; i++) {
            if (mzd_read_bit(Aug, i, i) == 0) {
                singular = 1;
                break;
            }
        }

        if (!singular) {
            // Success: The right half of the augmented matrix is now P^-1
            *out_P_inv = mzd_init(N, N);
            for(int r=0; r<N; r++) {
                for(int c=0; c<N; c++) {
                    mzd_write_bit(*out_P_inv, r, c, mzd_read_bit(Aug, r, N + c));
                }
            }
            break; // Valid P and P_inv found
        }
        // If singular, loop repeats with a new random P
    }
    
    mzd_free(Identity);
    mzd_free(Aug);
    *out_P = P;
}

/**
 * Applies a linear Input Encoding P_in to the quadratic round matrix M.
 * * @param M The round matrix containing flattened quadratic rows.
 * @param P_inv The inverse of the input encoding matrix.
 * @param N Dimension of the linear state.
 */
void wb_apply_input_encoding(mzd_t *M, mzd_t *P_inv, int N) {
    // Allocation for temporary calculation matrices
    mzd_t *Q = mzd_init(N, N);
    mzd_t *Q_new = mzd_init(N, N);
    mzd_t *P_inv_T = mzd_transpose(NULL, P_inv);
    mzd_t *Temp = mzd_init(N, N);
    
    int N_quad = M->ncols;
    mzd_t *M_new = mzd_init(N, N_quad);
    wb_mzd_set_zero(M_new);

    // Iterate over each output bit (row) of the white-box matrix
    for (int r = 0; r < N; r++) {
        wb_mzd_set_zero(Q);
        
        // 1. Reconstruct the quadratic coefficient matrix Q from the flattened row.
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                int val = mzd_read_bit(M, r, get_quad_idx(i, j, N));
                if (val) {
                    mzd_write_bit(Q, i, j, 1);
                }
            }
        }

        // 2. Apply Linear Transformation: Q_new = P_inv^T * Q * P_inv
        // Temp = Q * P_inv
        mzd_mul(Temp, Q, P_inv, 0);
        // Q_new = P_inv_T * Temp
        mzd_mul(Q_new, P_inv_T, Temp, 0);

        // 3. Fold and Flatten back to Upper-Triangular format
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                int val;
                if (i == j) {
                    // Diagonal terms (Linear terms in GF2): x_i * x_i
                    val = mzd_read_bit(Q_new, i, i);
                } else {
                    // Off-diagonal terms (Cross terms): x_i * x_j
                    // The coefficient is the sum of the symmetric positions
                    val = mzd_read_bit(Q_new, i, j) ^ mzd_read_bit(Q_new, j, i);
                }
                
                if (val) {
                    mzd_write_bit(M_new, r, get_quad_idx(i, j, N), 1);
                }
            }
        }
    }

    // Update the matrix M with the encoded version
    mzd_copy(M, M_new);

    // Cleanup
    mzd_free(Q); mzd_free(Q_new); mzd_free(P_inv_T); mzd_free(Temp); mzd_free(M_new);
}

// =============================================================
// Standard SIMON Component Generators
// =============================================================

uint16_t* simon32_key_schedule(const uint8_t *master_key, int T) {
    uint16_t *rk = (uint16_t*)malloc(sizeof(uint16_t) * T);
    uint16_t k[4];
    // Load key as Big Endian words to match standard test vectors
    for(int i=0; i<4; i++) {
        int base = (3 - i) * 2; 
        k[i] = ((uint16_t)master_key[base] << 8) | (uint16_t)master_key[base+1];
    }
    for(int i=0; i<4; i++) rk[i] = k[i];
    
    // Standard Z0 sequence
    const uint8_t z0_seq[] = {
        1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,
        1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0
    };
    for(int i=4; i<T; i++) {
        uint16_t tmp = rotr16(rk[i-1], 3);
        tmp ^= rk[i-3];
        tmp ^= rotr16(tmp, 1);
        uint16_t z_bit = z0_seq[(i-4) % 62];
        rk[i] = ~rk[i-4] ^ tmp ^ z_bit ^ 3;
    }
    return rk;
}

// =============================================================
// Context Management
// =============================================================

void wb_ctx_init(WB_Context *ctx, SimonMode mode) {
    ctx->mode = mode;
    srand(time(NULL));
    if (mode == SIMON_32_64) {
        ctx->n = 16;
        ctx->T = 32;
    }
    ctx->d = ctx->n;
    // Total state size: Constants(2) + Shares(4n) + Dummies(2d)
    ctx->N = 2 + 4 * ctx->n + 2 * ctx->d;
    ctx->N_quad = ctx->N * (ctx->N + 1) / 2;

    ctx->round_matrices = (mzd_t **)malloc(sizeof(mzd_t *) * ctx->T);
    for (int i = 0; i < ctx->T; i++) {
        ctx->round_matrices[i] = mzd_init(ctx->N, ctx->N_quad);
    }
    ctx->dummy_consts = (uint16_t*)malloc(sizeof(uint16_t) * ctx->T);
    ctx->ext_encoding_in = NULL;
    ctx->ext_encoding_out_inv = NULL;
}

void wb_ctx_free(WB_Context *ctx) {
    for (int i = 0; i < ctx->T; i++) {
        mzd_free(ctx->round_matrices[i]);
    }
    free(ctx->round_matrices);
    free(ctx->dummy_consts);
    if(ctx->ext_encoding_in) mzd_free(ctx->ext_encoding_in);
    if(ctx->ext_encoding_out_inv) mzd_free(ctx->ext_encoding_out_inv);
}

// =============================================================
// Layer Generators (Logic, Masking, Dummy, QAM)
// =============================================================

mzd_t* wb_gen_layer_constants(WB_Context *ctx) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int idx_c0c0 = get_quad_idx(0, 0, ctx->N);
    int idx_c1c1 = get_quad_idx(1, 1, ctx->N);
    int idx_c0c1 = get_quad_idx(0, 1, ctx->N);
    mzd_write_bit(M, 1, idx_c0c0, 1);
    mzd_write_bit(M, 1, idx_c1c1, 1);
    mzd_write_bit(M, 1, idx_c0c1, 1);
    return M;
}

mzd_t* wb_gen_layer_simon_lrmix(WB_Context *ctx) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int n = ctx->n; int N = ctx->N;
    // Offsets for Share 0 (U0, V0) and Share 1 (U1, V1)
    int off_u0 = 2; int off_v0 = 2 + n;
    int off_u1 = 2 + 2 * n; int off_v1 = 2 + 3 * n;
    
    // Generate masked logic for SIMON Round Function
    for (int i = 0; i < n; i++) {
        int row_u0 = off_u0 + i; int row_v0 = off_v0 + i;
        int row_u1 = off_u1 + i; int row_v1 = off_v1 + i;
        int idx_u0 = get_quad_idx(off_u0+i, off_u0+i, N);
        int idx_v0 = get_quad_idx(off_v0+i, off_v0+i, N);
        int idx_u1 = get_quad_idx(off_u1+i, off_u1+i, N);
        int idx_v1 = get_quad_idx(off_v1+i, off_v1+i, N);
        
        // Circular shifts
        int r1 = (i - 1 + n) % n; int r8 = (i - 8 + n) % n; int r2 = (i - 2 + n) % n;
        int u0_r1 = off_u0 + r1; int u0_r8 = off_u0 + r8; int u0_r2 = off_u0 + r2;
        int u1_r1 = off_u1 + r1; int u1_r8 = off_u1 + r8; int u1_r2 = off_u1 + r2;
        
        int f_lin_0 = get_quad_idx(u0_r2, u0_r2, N);
        int f_lin_1 = get_quad_idx(u1_r2, u1_r2, N);
        
        // Expanded AND gate: (A0+A1)(B0+B1) = A0B0 + A0B1 + A1B0 + A1B1
        int term_00 = get_quad_idx(u0_r1, u0_r8, N);
        int term_11 = get_quad_idx(u1_r1, u1_r8, N);
        int term_01 = get_quad_idx(u0_r1, u1_r8, N);
        int term_10 = get_quad_idx(u1_r1, u0_r8, N);
        
        // Share 0 logic
        mzd_write_bit(M, row_u0, idx_u0, 1); mzd_write_bit(M, row_u0, idx_v0, 1);
        mzd_write_bit(M, row_u0, f_lin_0, 1); mzd_write_bit(M, row_u0, term_00, 1);
        mzd_write_bit(M, row_u0, term_01, 1); mzd_write_bit(M, row_u0, term_10, 1);
        mzd_write_bit(M, row_v0, idx_v0, 1); mzd_write_bit(M, row_v0, f_lin_0, 1);
        mzd_write_bit(M, row_v0, term_00, 1); mzd_write_bit(M, row_v0, term_01, 1);
        mzd_write_bit(M, row_v0, term_10, 1);
        
        // Share 1 logic
        mzd_write_bit(M, row_u1, idx_u1, 1); mzd_write_bit(M, row_u1, idx_v1, 1);
        mzd_write_bit(M, row_u1, f_lin_1, 1); mzd_write_bit(M, row_u1, term_11, 1);
        mzd_write_bit(M, row_v1, idx_v1, 1); mzd_write_bit(M, row_v1, f_lin_1, 1);
        mzd_write_bit(M, row_v1, term_11, 1);
    }
    return M;
}

mzd_t* wb_gen_layer_key(WB_Context *ctx, uint16_t round_key) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int n = ctx->n; int N = ctx->N;
    int idx_c0c0 = get_quad_idx(0, 0, N);
    int idx_c1c1 = get_quad_idx(1, 1, N);
    int idx_c0c1 = get_quad_idx(0, 1, N);
    
    int off_u0 = 2; int off_v0 = 2 + n;
    int off_u1 = 2 + 2 * n; int off_v1 = 2 + 3 * n;
    
    // Split key into K0 (random) and K1 (rest)
    uint16_t K0 = rand() & 0xFFFF;
    uint16_t K1 = round_key ^ K0;
    
    for (int i = 0; i < n; i++) {
        if ((K0 >> i) & 1) {
            mzd_write_bit(M, off_u0 + i, idx_c0c0, 1); mzd_write_bit(M, off_u0 + i, idx_c1c1, 1); mzd_write_bit(M, off_u0 + i, idx_c0c1, 1);
            mzd_write_bit(M, off_v0 + i, idx_c0c0, 1); mzd_write_bit(M, off_v0 + i, idx_c1c1, 1); mzd_write_bit(M, off_v0 + i, idx_c0c1, 1);
        }
        if ((K1 >> i) & 1) {
            mzd_write_bit(M, off_u1 + i, idx_c0c0, 1); mzd_write_bit(M, off_u1 + i, idx_c1c1, 1); mzd_write_bit(M, off_u1 + i, idx_c0c1, 1);
            mzd_write_bit(M, off_v1 + i, idx_c0c0, 1); mzd_write_bit(M, off_v1 + i, idx_c1c1, 1); mzd_write_bit(M, off_v1 + i, idx_c0c1, 1);
        }
    }
    return M;
}

mzd_t* wb_gen_layer_dummy(WB_Context *ctx, int r) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int n = ctx->n; int N = ctx->N; int d = ctx->d;
    uint16_t kd0 = rand() & 0xFFFF;
    uint16_t kd1 = rand() & 0xFFFF;
    
    // Store expected difference for QAM check in next round
    ctx->dummy_consts[r] = kd0 ^ kd1;
    
    int off_d0 = 2 + 4 * n; int off_d1 = 2 + 4 * n + d; 
    int idx_c0c0 = get_quad_idx(0, 0, N);
    int idx_c1c1 = get_quad_idx(1, 1, N);
    int idx_c0c1 = get_quad_idx(0, 1, N);
    
    for (int i = 0; i < d; i++) {
        int row_d0 = off_d0 + i; int row_d1 = off_d1 + i;
        int r1 = (i - 1 + d) % d; int r8 = (i - 8 + d) % d; int r2 = (i - 2 + d) % d;
        int idx_d0_r1 = off_d0 + r1; int idx_d0_r8 = off_d0 + r8; int idx_d0_r2 = off_d0 + r2;
        int quad_f0 = get_quad_idx(idx_d0_r1, idx_d0_r8, N); int lin_f0  = get_quad_idx(idx_d0_r2, idx_d0_r2, N);
        int quad_f1 = get_quad_idx(off_d0 + r1, off_d0 + r8, N); // Same logic source (D0) for both branches
        int lin_f1  = get_quad_idx(off_d0 + r2, off_d0 + r2, N);
        
        int idx_R1 = get_quad_idx(off_d1 + i, off_d1 + i, N); // R term comes from D1
        
        // Unified Logic: Y = D1 ^ F(D0)
        // Both D0_next and D1_next use this Y
        mzd_write_bit(M, row_d0, idx_R1, 1); mzd_write_bit(M, row_d0, quad_f0, 1); mzd_write_bit(M, row_d0, lin_f0, 1);
        mzd_write_bit(M, row_d1, idx_R1, 1); mzd_write_bit(M, row_d1, quad_f0, 1); mzd_write_bit(M, row_d1, lin_f0, 1);
        
        // Add distinct dummy keys
        if ((kd0 >> i) & 1) {
            mzd_write_bit(M, row_d0, idx_c0c0, 1); mzd_write_bit(M, row_d0, idx_c1c1, 1); mzd_write_bit(M, row_d0, idx_c0c1, 1);
        }
        if ((kd1 >> i) & 1) {
            mzd_write_bit(M, row_d1, idx_c0c0, 1); mzd_write_bit(M, row_d1, idx_c1c1, 1); mzd_write_bit(M, row_d1, idx_c0c1, 1);
        }
    }
    return M;
}

mzd_t* wb_gen_layer_dense_mask(WB_Context *ctx) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int n = ctx->n; int d = ctx->d; int N_quad = ctx->N_quad;
    // Apply random noise that cancels out in XOR pairs (Share0^Share1, D0^D1, c0^c1)
    for (int col = 0; col < N_quad; col++) {
        if (rand() & 1) { mzd_write_bit(M, 0, col, 1); mzd_write_bit(M, 1, col, 1); }
    }
    int off_s0 = 2; int off_s1 = 2 + 2 * n;
    for (int i = 0; i < 2 * n; i++) {
        int row_s0 = off_s0 + i; int row_s1 = off_s1 + i;
        for (int col = 0; col < N_quad; col++) {
            if (rand() & 1) { mzd_write_bit(M, row_s0, col, 1); mzd_write_bit(M, row_s1, col, 1); }
        }
    }
    int off_d0 = 2 + 4 * n; int off_d1 = 2 + 4 * n + d;
    for (int i = 0; i < d; i++) {
        int row_d0 = off_d0 + i; int row_d1 = off_d1 + i;
        for (int col = 0; col < N_quad; col++) {
            if (rand() & 1) { mzd_write_bit(M, row_d0, col, 1); mzd_write_bit(M, row_d1, col, 1); }
        }
    }
    return M;
}

mzd_t* wb_gen_layer_qam(WB_Context *ctx, uint16_t target_diff) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int N = ctx->N; int n = ctx->n; int d = ctx->d;
    int idx_c0 = 0; int idx_c1 = 1;
    int off_d0 = 2 + 4 * n; int off_d1 = 2 + 4 * n + d;
    int *lambda = (int*)malloc(sizeof(int) * (N + 1));
    int *rho    = (int*)malloc(sizeof(int) * N);
    
    // Create quadratic annihilators L(x)*R(x) where L(x)=0 for valid inputs
    for (int row = 0; row < N; row++) {
        for(int k=0; k<=N; k++) lambda[k] = 0;
        
        // L_ctrl: c0 + c1 + 1 = 0
        if (rand() & 1) {
            lambda[idx_c0] ^= 1; lambda[idx_c1] ^= 1; lambda[N] ^= 1;
        }
        
        // L_dummy: D0_i + D1_i + C_i = 0
        for (int l = 0; l < d; l++) {
            if (rand() & 1) {
                lambda[off_d0 + l] ^= 1; lambda[off_d1 + l] ^= 1;
                if ((target_diff >> l) & 1) { lambda[idx_c0] ^= 1; lambda[idx_c1] ^= 1; }
            }
        }
        
        // R(x): Random linear form
        for(int q=0; q<N; q++) rho[q] = rand() & 1;
        
        // Q(x) = L(x)*R(x)
        for (int q = 0; q < N; q++) {
            if (rho[q] == 0) continue;
            // Linear terms from const part of L(x)
            if (lambda[N] == 1) mzd_flip_bit(M, row, get_quad_idx(q, q, N));
            // Quadratic terms from vars in L(x)
            for (int p = 0; p < N; p++) {
                if (lambda[p] == 1) mzd_flip_bit(M, row, get_quad_idx(p, q, N));
            }
        }
    }
    free(lambda); free(rho);
    return M;
}

// =============================================================
// Main Generation Procedure
// =============================================================

void wb_generate_matrix(WB_Context *ctx, const uint8_t *master_key) {
    printf("[WB-Gen] Generating White-Box Matrix...\n");
    uint16_t *rks = simon32_key_schedule(master_key, ctx->T);

    int N = ctx->N;
    
    // P_in_inv tracks the inverse of the input encoding for the current round
    mzd_t *P_in_inv = NULL; 

    for (int r = 0; r < ctx->T; r++) {
        mzd_t *M_final = ctx->round_matrices[r];
        
        // 1. Generate Core Layers
        mzd_t *M_const = wb_gen_layer_constants(ctx);
        mzd_t *M_logic = wb_gen_layer_simon_lrmix(ctx);
        mzd_t *M_key   = wb_gen_layer_key(ctx, rks[r]);
        mzd_t *M_dummy = wb_gen_layer_dummy(ctx, r);
        mzd_t *M_drm   = wb_gen_layer_dense_mask(ctx);
        uint16_t constraint = (r == 0) ? 0 : ctx->dummy_consts[r - 1];
        mzd_t *M_qam   = wb_gen_layer_qam(ctx, constraint);

        mzd_add(M_final, M_final, M_const);
        mzd_add(M_final, M_final, M_logic);
        mzd_add(M_final, M_final, M_key);
        mzd_add(M_final, M_final, M_dummy);
        mzd_add(M_final, M_final, M_drm);
        mzd_add(M_final, M_final, M_qam);

        mzd_free(M_const); mzd_free(M_logic); mzd_free(M_key);
        mzd_free(M_dummy); mzd_free(M_drm); mzd_free(M_qam);

        // 2. Setup Input Encoding P_in
        if (r == 0) {
            // First round uses External Input Encoding
            mzd_t *P_ext_in = NULL;
            wb_gen_invertible_matrix(N, &P_ext_in, &P_in_inv);
            ctx->ext_encoding_in = P_ext_in; 
        } 
        // Else: P_in_inv is inherited from previous round's P_out_inv
        
        // 3. Setup Output Encoding P_out
        mzd_t *P_out = NULL;
        mzd_t *P_out_inv = NULL;

        if (r == ctx->T - 1) {
            // Last round uses External Output Encoding
            wb_gen_invertible_matrix(N, &P_out, &P_out_inv);
            ctx->ext_encoding_out_inv = P_out_inv;
        } else {
            // Internal rounds use random bijections
            wb_gen_invertible_matrix(N, &P_out, &P_out_inv);
        }

        // 4. Apply Encodings to Core Matrix
        // Transform Input: Q' = P_in_inv^T * Q * P_in_inv
        wb_apply_input_encoding(M_final, P_in_inv, N);

        // Transform Output: M' = P_out * M
        mzd_t *M_temp = mzd_init(N, ctx->N_quad);
        mzd_mul(M_temp, P_out, M_final, 0);
        mzd_copy(M_final, M_temp);
        mzd_free(M_temp);

        // 5. Transfer encoding to next round
        mzd_free(P_in_inv); // Done with current input inverse
        P_in_inv = P_out_inv; // Next round's input inv is current output inv
        mzd_free(P_out); // Matrix P_out no longer needed
    }
    free(rks);
    printf("[WB-Gen] Matrix generation complete.\n");
}

// =============================================================
// Encryption Execution
// =============================================================

mzd_t* wb_map_input(WB_Context *ctx, const uint8_t *in) {
    int N = ctx->N;
    int n = ctx->n;
    int d = ctx->d;

    // 1. Decode bytes to words
    uint16_t L = ((uint16_t)in[0] << 8) | (uint16_t)in[1];
    uint16_t R = ((uint16_t)in[2] << 8) | (uint16_t)in[3];
    uint16_t U = L; 
    uint16_t V = L ^ R;
    
    mzd_t *x = mzd_init(N, 1);
    
    // 2. Initialize Constants (c0 != c1)
    int c0 = rand() & 1; 
    int c1 = 1 ^ c0;
    mzd_write_bit(x, 0, 0, c0); 
    mzd_write_bit(x, 1, 0, c1);
    
    // 3. Initialize Real Shares (U0^U1=U, V0^V1=V)
    uint16_t U0 = rand() & 0xFFFF; uint16_t U1 = U ^ U0;
    uint16_t V0 = rand() & 0xFFFF; uint16_t V1 = V ^ V0;
    
    for(int i=0; i<n; i++) {
        mzd_write_bit(x, 2 + i, 0, (U0 >> i) & 1);
        mzd_write_bit(x, 2 + n + i, 0, (V0 >> i) & 1);
        mzd_write_bit(x, 2 + 2*n + i, 0, (U1 >> i) & 1);
        mzd_write_bit(x, 2 + 3*n + i, 0, (V1 >> i) & 1);
    }
    
    // 4. Initialize Dummy Shares (Initialized Equal)
    uint16_t initial_dummy = rand() & 0xFFFF;
    int off_d0 = 2 + 4 * n; 
    int off_d1 = 2 + 4 * n + d;
    for(int i=0; i<d; i++) {
        int bit = (initial_dummy >> i) & 1;
        mzd_write_bit(x, off_d0 + i, 0, bit);
        mzd_write_bit(x, off_d1 + i, 0, bit);
    }

    // 5. Apply External Input Encoding: x_enc = P_ext_in * x
    mzd_t *x_enc = mzd_init(N, 1);
    mzd_mul(x_enc, ctx->ext_encoding_in, x, 0);
    mzd_copy(x, x_enc);
    mzd_free(x_enc);

    return x;
}

void wb_execute_rounds(WB_Context *ctx, mzd_t *state) {
    int N = ctx->N;
    int N_quad = ctx->N_quad;
    
    mzd_t *x_quad = mzd_init(N_quad, 1);
    mzd_t *y = mzd_init(N, 1);

    for (int r = 0; r < ctx->T; r++) {
        // Expand state to quadratic form
        int q_idx = 0;
        for (int i = 0; i < N; i++) {
            int bit_i = mzd_read_bit(state, i, 0);
            for (int j = i; j < N; j++) {
                int bit_j = mzd_read_bit(state, j, 0);
                mzd_write_bit(x_quad, q_idx, 0, bit_i & bit_j);
                q_idx++;
            }
        }
        
        // Execute Round: y = M_r * x_quad
        mzd_mul(y, ctx->round_matrices[r], x_quad, 0);
        
        // Update state
        mzd_copy(state, y);
    }

    mzd_free(x_quad);
    mzd_free(y);
}

void wb_map_output(WB_Context *ctx, mzd_t *state, uint8_t *out) {
    int N = ctx->N;
    int n = ctx->n;

    // 1. Apply External Output Decoding: x_raw = P_ext_out_inv * state
    mzd_t *x_raw = mzd_init(N, 1);
    mzd_mul(x_raw, ctx->ext_encoding_out_inv, state, 0);
    mzd_copy(state, x_raw);
    mzd_free(x_raw);

    // 2. Decode Output (Recombine Shares)
    uint16_t U_final = 0, V_final = 0;
    for(int i=0; i<n; i++) {
        int u_bit = mzd_read_bit(state, 2 + i, 0) ^ mzd_read_bit(state, 2 + 2*n + i, 0);
        int v_bit = mzd_read_bit(state, 2 + n + i, 0) ^ mzd_read_bit(state, 2 + 3*n + i, 0);
        if(u_bit) U_final |= (1 << i);
        if(v_bit) V_final |= (1 << i);
    }

    uint16_t L = U_final; 
    uint16_t R = U_final ^ V_final;
    
    // 3. Store Bytes (Big Endian)
    out[0] = (L >> 8) & 0xFF; out[1] = L & 0xFF;
    out[2] = (R >> 8) & 0xFF; out[3] = R & 0xFF;

    // 4. Free the state vector
    mzd_free(state);
}