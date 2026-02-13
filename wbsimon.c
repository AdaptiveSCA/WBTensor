#include "wbsimon.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// --- Helper Functions ---

/**
 * Calculates the index in the flattened quadratic vector for the term x_i * x_j.
 * Assumes a Row-Major Upper Triangular mapping.
 * Ensures i <= j before calculation.
 */
int get_quad_idx(int i, int j, int N) {
    if (i > j) { int tmp = i; i = j; j = tmp; }
    // Formula: Index = (Start of Row i) + (Offset in Row i)
    // Start of Row i = Sum(k=0 to i-1) of (N-k)
    return (i * N) - (i * (i + 1)) / 2 + j;
}

/**
 * Performs a circular right shift on a 16-bit word.
 */
static inline uint16_t rotr16(uint16_t x, int n) {
    return (x >> n) | (x << (16 - n));
}

// --- Key Schedule ---

/**
 * Generates round keys for SIMON 32/64.
 * Handles endianness conversion to match standard test vectors.
 */
uint16_t* simon32_key_schedule(const uint8_t *master_key, int T) {
    uint16_t *rk = (uint16_t*)malloc(sizeof(uint16_t) * T);
    uint16_t k[4];
    
    // Load Master Key words in reverse order to match NIST vector representation.
    // Bytes are loaded as Big Endian words (High Byte, Low Byte).
    for(int i=0; i<4; i++) {
        int base = (3 - i) * 2; 
        k[i] = ((uint16_t)master_key[base] << 8) | (uint16_t)master_key[base+1];
    }
    
    // Initialize first 4 round keys
    for(int i=0; i<4; i++) rk[i] = k[i];
    
    // SIMON Z0 Sequence (Period 62)
    const uint8_t z0_seq[] = {
        1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,
        1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0
    };

    // Generate remaining round keys
    for(int i=4; i<T; i++) {
        uint16_t tmp = rotr16(rk[i-1], 3);
        tmp ^= rk[i-3];
        tmp ^= rotr16(tmp, 1);
        uint16_t z_bit = z0_seq[(i-4) % 62];
        rk[i] = ~rk[i-4] ^ tmp ^ z_bit ^ 3;
    }
    return rk;
}

// --- Context Management ---

void wb_ctx_init(WB_Context *ctx, SimonMode mode) {
    ctx->mode = mode;
    srand(time(NULL)); // Initialize RNG

    if (mode == SIMON_32_64) {
        ctx->n = 16;
        ctx->T = 32;
    }
    
    // Configuration: d (Dummy size) equals n
    ctx->d = ctx->n;

    // Linear State Size N:
    // 2 (Constants) + 2n (Share 0) + 2n (Share 1) + 2d (Dummies)
    ctx->N = 2 + 4 * ctx->n + 2 * ctx->d;
    
    // Quadratic State Size
    ctx->N_quad = ctx->N * (ctx->N + 1) / 2;

    // Allocate matrices
    ctx->round_matrices = (mzd_t **)malloc(sizeof(mzd_t *) * ctx->T);
    for (int i = 0; i < ctx->T; i++) {
        ctx->round_matrices[i] = mzd_init(ctx->N, ctx->N_quad);
    }
    
    ctx->dummy_consts = (uint16_t*)malloc(sizeof(uint16_t) * ctx->T);
}

void wb_ctx_free(WB_Context *ctx) {
    for (int i = 0; i < ctx->T; i++) {
        mzd_free(ctx->round_matrices[i]);
    }
    free(ctx->round_matrices);
    free(ctx->dummy_consts);
}

// --- Layer Generators ---

/**
 * LAYER: Constants Preservation
 * Ensures that the constant bits c0 and c1 evolve correctly.
 * Logic: c0' = 0, c1' = c0 + c1 = 1 (Identity sum).
 */
mzd_t* wb_gen_layer_constants(WB_Context *ctx) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    
    // Indices for quadratic terms c0^2, c1^2, c0*c1
    int idx_c0c0 = get_quad_idx(0, 0, ctx->N);
    int idx_c1c1 = get_quad_idx(1, 1, ctx->N);
    int idx_c0c1 = get_quad_idx(0, 1, ctx->N);

    // Row 1 (Target c1) accumulates c0 + c1 + c0c1
    mzd_write_bit(M, 1, idx_c0c0, 1);
    mzd_write_bit(M, 1, idx_c1c1, 1);
    mzd_write_bit(M, 1, idx_c0c1, 1);

    return M;
}

/**
 * LAYER: SIMON Logic with L-R Mixing and Boolean Masking
 * Implements the masked SIMON round function distributed across Share 0 and Share 1.
 * * Masking Strategy for AND gate A & B:
 * (A0 + A1)(B0 + B1) = A0B0 + A0B1 + A1B0 + A1B1
 * Share 0 computes: A0B0 + A0B1 + A1B0
 * Share 1 computes: A1B1
 */
mzd_t* wb_gen_layer_simon_lrmix(WB_Context *ctx) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int n = ctx->n;
    int N = ctx->N;
    
    // Row Offsets for Share 0 and Share 1
    int off_u0 = 2;
    int off_v0 = 2 + n;
    int off_u1 = 2 + 2 * n;
    int off_v1 = 2 + 3 * n;

    for (int i = 0; i < n; i++) {
        // Output Rows
        int row_u0 = off_u0 + i;
        int row_v0 = off_v0 + i;
        int row_u1 = off_u1 + i;
        int row_v1 = off_v1 + i;

        // Linear input indices (Identity terms)
        int idx_u0 = get_quad_idx(off_u0+i, off_u0+i, N);
        int idx_v0 = get_quad_idx(off_v0+i, off_v0+i, N);
        int idx_u1 = get_quad_idx(off_u1+i, off_u1+i, N);
        int idx_v1 = get_quad_idx(off_v1+i, off_v1+i, N);

        // Rotation indices for F(U) calculation
        int r1 = (i - 1 + n) % n;
        int r8 = (i - 8 + n) % n;
        int r2 = (i - 2 + n) % n;

        // Share 0 Rotations
        int u0_r1 = off_u0 + r1;
        int u0_r8 = off_u0 + r8;
        int u0_r2 = off_u0 + r2;

        // Share 1 Rotations
        int u1_r1 = off_u1 + r1;
        int u1_r8 = off_u1 + r8;
        int u1_r2 = off_u1 + r2;

        // --- F(U) Construction ---
        
        // Linear part: U <<< 2 (Distributed linearly)
        int f_lin_0 = get_quad_idx(u0_r2, u0_r2, N);
        int f_lin_1 = get_quad_idx(u1_r2, u1_r2, N);

        // Quadratic part: (U <<< 1) & (U <<< 8)
        // Expanded into 4 cross-terms
        int term_00 = get_quad_idx(u0_r1, u0_r8, N); // A0 & B0
        int term_11 = get_quad_idx(u1_r1, u1_r8, N); // A1 & B1
        int term_01 = get_quad_idx(u0_r1, u1_r8, N); // A0 & B1
        int term_10 = get_quad_idx(u1_r1, u0_r8, N); // A1 & B0

        // --- Matrix Population ---

        // Share 0 Update: Includes 3 terms of the expansion (00, 01, 10)
        // U0' = U0 + V0 + F_lin(U0) + QuadTerms(0)
        mzd_write_bit(M, row_u0, idx_u0, 1);
        mzd_write_bit(M, row_u0, idx_v0, 1);
        mzd_write_bit(M, row_u0, f_lin_0, 1);
        mzd_write_bit(M, row_u0, term_00, 1);
        mzd_write_bit(M, row_u0, term_01, 1);
        mzd_write_bit(M, row_u0, term_10, 1);

        // V0' = V0 + F_lin(U0) + QuadTerms(0)
        mzd_write_bit(M, row_v0, idx_v0, 1);
        mzd_write_bit(M, row_v0, f_lin_0, 1);
        mzd_write_bit(M, row_v0, term_00, 1);
        mzd_write_bit(M, row_v0, term_01, 1);
        mzd_write_bit(M, row_v0, term_10, 1);

        // Share 1 Update: Includes 1 term of the expansion (11)
        // U1' = U1 + V1 + F_lin(U1) + QuadTerms(1)
        mzd_write_bit(M, row_u1, idx_u1, 1);
        mzd_write_bit(M, row_u1, idx_v1, 1);
        mzd_write_bit(M, row_u1, f_lin_1, 1);
        mzd_write_bit(M, row_u1, term_11, 1);

        // V1' = V1 + F_lin(U1) + QuadTerms(1)
        mzd_write_bit(M, row_v1, idx_v1, 1);
        mzd_write_bit(M, row_v1, f_lin_1, 1);
        mzd_write_bit(M, row_v1, term_11, 1);
    }
    return M;
}

/**
 * LAYER: Key Addition with Splitting
 * Splits the round key K into K0 and K1.
 * K0 is added to Share 0 rows.
 * K1 is added to Share 1 rows.
 * Keys are embedded via the constant columns (c0, c1).
 */
mzd_t* wb_gen_layer_key(WB_Context *ctx, uint16_t round_key) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int n = ctx->n;
    int N = ctx->N;

    int idx_c0c0 = get_quad_idx(0, 0, N);
    int idx_c1c1 = get_quad_idx(1, 1, N);
    int idx_c0c1 = get_quad_idx(0, 1, N);

    int off_u0 = 2;
    int off_v0 = 2 + n;
    int off_u1 = 2 + 2 * n;
    int off_v1 = 2 + 3 * n;

    // Generate random key mask
    uint16_t K0 = rand() & 0xFFFF;
    uint16_t K1 = round_key ^ K0;

    for (int i = 0; i < n; i++) {
        // Apply K0 to Share 0
        int k0_bit = (K0 >> i) & 1;
        if (k0_bit) {
            // Add to U0
            mzd_write_bit(M, off_u0 + i, idx_c0c0, 1);
            mzd_write_bit(M, off_u0 + i, idx_c1c1, 1);
            mzd_write_bit(M, off_u0 + i, idx_c0c1, 1);
            // Add to V0
            mzd_write_bit(M, off_v0 + i, idx_c0c0, 1);
            mzd_write_bit(M, off_v0 + i, idx_c1c1, 1);
            mzd_write_bit(M, off_v0 + i, idx_c0c1, 1);
        }

        // Apply K1 to Share 1
        int k1_bit = (K1 >> i) & 1;
        if (k1_bit) {
            // Add to U1
            mzd_write_bit(M, off_u1 + i, idx_c0c0, 1);
            mzd_write_bit(M, off_u1 + i, idx_c1c1, 1);
            mzd_write_bit(M, off_u1 + i, idx_c0c1, 1);
            // Add to V1
            mzd_write_bit(M, off_v1 + i, idx_c0c0, 1);
            mzd_write_bit(M, off_v1 + i, idx_c1c1, 1);
            mzd_write_bit(M, off_v1 + i, idx_c0c1, 1);
        }
    }
    return M;
}

/**
 * LAYER: Dummy Inputs (Pairwise Redundancy)
 * Creates two independent SIMON-like execution paths (D0 and D1).
 * D0 and D1 are initialized identically but receive different random dummy keys.
 * The XOR difference (D0 ^ D1) remains a predictable constant (C) for QAM verification.
 */
mzd_t* wb_gen_layer_dummy(WB_Context *ctx, int r) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int n = ctx->n;
    int N = ctx->N;
    int d = ctx->d;

    // Generate random dummy keys
    uint16_t kd0 = rand() & 0xFFFF;
    uint16_t kd1 = rand() & 0xFFFF;
    
    // Store expected difference for verification
    ctx->dummy_consts[r] = kd0 ^ kd1;

    int off_d0 = 2 + 4 * n;     
    int off_d1 = 2 + 4 * n + d; 

    int idx_c0c0 = get_quad_idx(0, 0, N);
    int idx_c1c1 = get_quad_idx(1, 1, N);
    int idx_c0c1 = get_quad_idx(0, 1, N);

    for (int i = 0; i < d; i++) {
        int row_d0 = off_d0 + i;
        int row_d1 = off_d1 + i;
        
        // Rotation indices for dummy logic (Simulating SIMON round)
        int r1 = (i - 1 + d) % d;
        int r8 = (i - 8 + d) % d;
        int r2 = (i - 2 + d) % d;

        // D0 Logic
        int idx_d0_r1 = off_d0 + r1;
        int idx_d0_r8 = off_d0 + r8;
        int idx_d0_r2 = off_d0 + r2;
        int quad_f0 = get_quad_idx(idx_d0_r1, idx_d0_r8, N);
        int lin_f0  = get_quad_idx(idx_d0_r2, idx_d0_r2, N);

        // D1 Logic
        int idx_d1_r1 = off_d1 + r1;
        int idx_d1_r8 = off_d1 + r8;
        int idx_d1_r2 = off_d1 + r2;
        int quad_f1 = get_quad_idx(idx_d1_r1, idx_d1_r8, N);
        int lin_f1  = get_quad_idx(idx_d1_r2, idx_d1_r2, N);

        // Apply Logic: D_new = D + F(D)
        mzd_write_bit(M, row_d0, quad_f0, 1);
        mzd_write_bit(M, row_d0, lin_f0, 1);
        mzd_write_bit(M, row_d0, get_quad_idx(off_d0+i, off_d0+i, N), 1); 
        
        mzd_write_bit(M, row_d1, quad_f1, 1);
        mzd_write_bit(M, row_d1, lin_f1, 1);
        mzd_write_bit(M, row_d1, get_quad_idx(off_d1+i, off_d1+i, N), 1); 

        // Apply Dummy Keys
        if ((kd0 >> i) & 1) {
            mzd_write_bit(M, row_d0, idx_c0c0, 1);
            mzd_write_bit(M, row_d0, idx_c1c1, 1);
            mzd_write_bit(M, row_d0, idx_c0c1, 1);
        }
        if ((kd1 >> i) & 1) {
            mzd_write_bit(M, row_d1, idx_c0c0, 1);
            mzd_write_bit(M, row_d1, idx_c1c1, 1);
            mzd_write_bit(M, row_d1, idx_c0c1, 1);
        }
    }
    return M;
}

/**
 * LAYER: Dense Random Masking (QAM)
 * Injects high-entropy random noise into the system.
 * This noise is generated per complementary row pair such that it cancels out
 * when the pair is combined (XORed).
 * * Pairs targeted:
 * 1. (c0, c1)
 * 2. (Share0_Row_i, Share1_Row_i)
 * 3. (Dummy0_Row_i, Dummy1_Row_i)
 */
mzd_t* wb_gen_layer_dense_mask(WB_Context *ctx) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int n = ctx->n;
    int d = ctx->d;
    int N_quad = ctx->N_quad;

    // 1. Mask Constants (c0, c1)
    for (int col = 0; col < N_quad; col++) {
        if (rand() & 1) {
            mzd_write_bit(M, 0, col, 1);
            mzd_write_bit(M, 1, col, 1);
        }
    }

    // 2. Mask Real Shares (Share 0, Share 1)
    int off_s0 = 2;
    int off_s1 = 2 + 2 * n;
    
    // Covers both U and V regions
    for (int i = 0; i < 2 * n; i++) {
        int row_s0 = off_s0 + i;
        int row_s1 = off_s1 + i;

        for (int col = 0; col < N_quad; col++) {
            if (rand() & 1) {
                mzd_write_bit(M, row_s0, col, 1);
                mzd_write_bit(M, row_s1, col, 1);
            }
        }
    }

    // 3. Mask Dummy Shares (Dummy 0, Dummy 1)
    int off_d0 = 2 + 4 * n;
    int off_d1 = 2 + 4 * n + d;

    for (int i = 0; i < d; i++) {
        int row_d0 = off_d0 + i;
        int row_d1 = off_d1 + i;

        for (int col = 0; col < N_quad; col++) {
            if (rand() & 1) {
                mzd_write_bit(M, row_d0, col, 1);
                mzd_write_bit(M, row_d1, col, 1);
            }
        }
    }

    return M;
}

// --- Main Generator ---

void wb_generate_matrix(WB_Context *ctx, const uint8_t *master_key) {
    printf("[WB-Gen] Generating White-Box Matrices...\n");
    uint16_t *rks = simon32_key_schedule(master_key, ctx->T);

    for (int r = 0; r < ctx->T; r++) {
        mzd_t *M_final = ctx->round_matrices[r];
        
        // Generate independent layers
        mzd_t *M_const = wb_gen_layer_constants(ctx);
        mzd_t *M_logic = wb_gen_layer_simon_lrmix(ctx);
        mzd_t *M_key   = wb_gen_layer_key(ctx, rks[r]);
        mzd_t *M_dummy = wb_gen_layer_dummy(ctx, r);
        mzd_t *M_mask  = wb_gen_layer_dense_mask(ctx);

        // Fuse layers (XOR accumulation)
        mzd_add(M_final, M_final, M_const);
        mzd_add(M_final, M_final, M_logic);
        mzd_add(M_final, M_final, M_key);
        mzd_add(M_final, M_final, M_dummy);
        mzd_add(M_final, M_final, M_mask);

        // Cleanup
        mzd_free(M_const);
        mzd_free(M_logic);
        mzd_free(M_key);
        mzd_free(M_dummy);
        mzd_free(M_mask);
    }
    free(rks);
    printf("[WB-Gen] Matrix generation complete.\n");
}

// --- Execution ---

void wb_encrypt(WB_Context *ctx, const uint8_t *in, uint8_t *out) {
    int N = ctx->N;
    int N_quad = ctx->N_quad;
    int n = ctx->n;
    int d = ctx->d;

    // 1. Input Loading (Big Endian Load to match standard vectors)
    uint16_t L = ((uint16_t)in[0] << 8) | (uint16_t)in[1];
    uint16_t R = ((uint16_t)in[2] << 8) | (uint16_t)in[3];

    // 2. Map to Internal L-R Mixed State
    uint16_t U = L;
    uint16_t V = L ^ R;

    // 3. Initialize White-Box State Vector
    mzd_t *x = mzd_init(N, 1);
    
    // A. Randomize Constants (c0 + c1 = 1)
    int c0 = rand() & 1;
    int c1 = 1 ^ c0;
    mzd_write_bit(x, 0, 0, c0);
    mzd_write_bit(x, 1, 0, c1);

    // B. Create Real Shares (Masking)
    uint16_t U0 = rand() & 0xFFFF;
    uint16_t U1 = U ^ U0;
    uint16_t V0 = rand() & 0xFFFF;
    uint16_t V1 = V ^ V0;

    // Write Share 0
    for(int i=0; i<n; i++) {
        mzd_write_bit(x, 2 + i, 0, (U0 >> i) & 1);
        mzd_write_bit(x, 2 + n + i, 0, (V0 >> i) & 1);
    }
    // Write Share 1
    for(int i=0; i<n; i++) {
        mzd_write_bit(x, 2 + 2*n + i, 0, (U1 >> i) & 1);
        mzd_write_bit(x, 2 + 3*n + i, 0, (V1 >> i) & 1);
    }
    
    // C. Initialize Dummy Shares
    uint16_t initial_dummy = rand() & 0xFFFF;
    int off_d0 = 2 + 4 * n;
    int off_d1 = 2 + 4 * n + d;
    
    for(int i=0; i<d; i++) {
        int bit = (initial_dummy >> i) & 1;
        mzd_write_bit(x, off_d0 + i, 0, bit);
        mzd_write_bit(x, off_d1 + i, 0, bit);
    }

    // 4. Execution Loop
    mzd_t *x_quad = mzd_init(N_quad, 1);
    mzd_t *y = mzd_init(N, 1);

    for (int r = 0; r < ctx->T; r++) {
        // Expand: x -> x (tensor) x
        int q_idx = 0;
        for (int i = 0; i < N; i++) {
            int bit_i = mzd_read_bit(x, i, 0);
            for (int j = i; j < N; j++) {
                int bit_j = mzd_read_bit(x, j, 0);
                mzd_write_bit(x_quad, q_idx, 0, bit_i & bit_j);
                q_idx++;
            }
        }
        // Multiply: y = M * x_quad
        mzd_mul(y, ctx->round_matrices[r], x_quad, 0);
        
        // Update: x = y
        mzd_copy(x, y);
    }

    // 5. Output Decoding and Unmapping
    // Recombine Shares: U = U0 ^ U1, V = V0 ^ V1
    uint16_t U_final = 0, V_final = 0;
    
    for(int i=0; i<n; i++) {
        int u_bit = mzd_read_bit(x, 2 + i, 0) ^ mzd_read_bit(x, 2 + 2*n + i, 0);
        int v_bit = mzd_read_bit(x, 2 + n + i, 0) ^ mzd_read_bit(x, 2 + 3*n + i, 0);
        
        if(u_bit) U_final |= (1 << i);
        if(v_bit) V_final |= (1 << i);
    }

    L = U_final;
    R = U_final ^ V_final;

    // Big Endian Store
    out[0] = (L >> 8) & 0xFF;
    out[1] = L & 0xFF;
    out[2] = (R >> 8) & 0xFF;
    out[3] = R & 0xFF;

    mzd_free(x);
    mzd_free(x_quad);
    mzd_free(y);
}