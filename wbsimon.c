#include "wbsimon.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// --- Helpers ---
int get_quad_idx(int i, int j, int N) {
    if (i > j) { int tmp = i; i = j; j = tmp; }
    return (i * N) - (i * (i + 1)) / 2 + j;
}

static inline uint16_t rotr16(uint16_t x, int n) {
    return (x >> n) | (x << (16 - n));
}

// --- Key Schedule ---
uint16_t* simon32_key_schedule(const uint8_t *master_key, int T) {
    uint16_t *rk = (uint16_t*)malloc(sizeof(uint16_t) * T);
    uint16_t k[4];
    for(int i=0; i<4; i++) {
        int base = (3 - i) * 2; 
        k[i] = ((uint16_t)master_key[base] << 8) | (uint16_t)master_key[base+1];
    }
    for(int i=0; i<4; i++) rk[i] = k[i];
    
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

// --- Context ---

void wb_ctx_init(WB_Context *ctx, SimonMode mode) {
    ctx->mode = mode;
    srand(time(NULL));

    if (mode == SIMON_32_64) {
        ctx->n = 16;
        ctx->T = 32;
    }
    
    // Config: d = n
    ctx->d = ctx->n;

    // N = 2(c0,c1) + 4n(Real) + 2d(Dummy)
    // Layout: 
    // [0..1]: c0, c1
    // [2..2+n-1]: U
    // [2+n..2+2n-1]: V
    // [2+2n..2+4n-1]: Empty (Padding)
    // [2+4n..2+4n+d-1]: Dummy 0
    // [2+4n+d..2+4n+2d-1]: Dummy 1
    ctx->N = 2 + 4 * ctx->n + 2 * ctx->d;
    ctx->N_quad = ctx->N * (ctx->N + 1) / 2;

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

mzd_t* wb_gen_layer_constants(WB_Context *ctx) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    // c1' = c0 + c1
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
    int n = ctx->n;
    int N = ctx->N;
    int off_u = 2;
    int off_v = 2 + n;

    for (int i = 0; i < n; i++) {
        int row_u = off_u + i;
        int row_v = off_v + i;
        int lin_u = get_quad_idx(off_u+i, off_u+i, N);
        int lin_v = get_quad_idx(off_v+i, off_v+i, N);

        int r1 = (i - 1 + n) % n;
        int r8 = (i - 8 + n) % n;
        int r2 = (i - 2 + n) % n;
        int quad_f = get_quad_idx(off_u + r1, off_u + r8, N); 
        int lin_f  = get_quad_idx(off_u + r2, off_u + r2, N); 

        // U' += U + V + F(U)
        mzd_write_bit(M, row_u, lin_u, 1);
        mzd_write_bit(M, row_u, lin_v, 1);
        mzd_write_bit(M, row_u, lin_f, 1);
        mzd_write_bit(M, row_u, quad_f, 1);

        // V' += V + F(U)
        mzd_write_bit(M, row_v, lin_v, 1);
        mzd_write_bit(M, row_v, lin_f, 1);
        mzd_write_bit(M, row_v, quad_f, 1);
    }
    return M;
}

mzd_t* wb_gen_layer_key(WB_Context *ctx, uint16_t round_key) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int n = ctx->n;
    int N = ctx->N;
    int idx_c0c0 = get_quad_idx(0, 0, N);
    int idx_c1c1 = get_quad_idx(1, 1, N);
    int idx_c0c1 = get_quad_idx(0, 1, N);
    int off_u = 2;
    int off_v = 2 + n;

    for (int i = 0; i < n; i++) {
        int k_bit = (round_key >> i) & 1;
        if (k_bit) {
            mzd_write_bit(M, off_u + i, idx_c0c0, 1);
            mzd_write_bit(M, off_u + i, idx_c1c1, 1);
            mzd_write_bit(M, off_u + i, idx_c0c1, 1);
            mzd_write_bit(M, off_v + i, idx_c0c0, 1);
            mzd_write_bit(M, off_v + i, idx_c1c1, 1);
            mzd_write_bit(M, off_v + i, idx_c0c1, 1);
        }
    }
    return M;
}

/**
 * 4. Dummy Layer
 * Logic:
 * D0' = F_simon(D0) ^ K_dummy0
 * D1' = F_simon(D1) ^ K_dummy1
 * * Since D0_in = D1_in (initialized equal), F(D0) = F(D1).
 * Thus, D0' ^ D1' = K_dummy0 ^ K_dummy1 = Constant C^{(r)}
 */
mzd_t* wb_gen_layer_dummy(WB_Context *ctx, int r) {
    mzd_t *M = mzd_init(ctx->N, ctx->N_quad);
    int n = ctx->n;
    int N = ctx->N;
    int d = ctx->d;

    // 
    uint16_t kd0 = rand() & 0xFFFF;
    uint16_t kd1 = rand() & 0xFFFF;
    
    // 
    ctx->dummy_consts[r] = kd0 ^ kd1;

    // Dummy 
    int off_d0 = 2 + 4 * n;     
    int off_d1 = 2 + 4 * n + d; 

    // Constants indices for Key addition
    int idx_c0c0 = get_quad_idx(0, 0, N);
    int idx_c1c1 = get_quad_idx(1, 1, N);
    int idx_c0c1 = get_quad_idx(0, 1, N);

    // 对每一位 dummy bit 构造逻辑
    for (int i = 0; i < d; i++) {
        int row_d0 = off_d0 + i;
        int row_d1 = off_d1 + i;
        
        // --- 1. Cipher Logic (SIMON F(x)) ---
        
        // Indices relative to the specific dummy block
        int r1 = (i - 1 + d) % d;
        int r8 = (i - 8 + d) % d;
        int r2 = (i - 2 + d) % d;

        // D0 Logic
        int idx_d0_r1 = off_d0 + r1;
        int idx_d0_r8 = off_d0 + r8;
        int idx_d0_r2 = off_d0 + r2;
        int quad_f0 = get_quad_idx(idx_d0_r1, idx_d0_r8, N);
        int lin_f0  = get_quad_idx(idx_d0_r2, idx_d0_r2, N);

        // D1 Logic (Exactly same structure, just shifted offsets)
        int idx_d1_r1 = off_d1 + r1;
        int idx_d1_r8 = off_d1 + r8;
        int idx_d1_r2 = off_d1 + r2;
        int quad_f1 = get_quad_idx(idx_d1_r1, idx_d1_r8, N);
        int lin_f1  = get_quad_idx(idx_d1_r2, idx_d1_r2, N);

        // Apply Logic: D_new = F(D) 
        // D0' += F(D0)
        mzd_write_bit(M, row_d0, quad_f0, 1);
        mzd_write_bit(M, row_d0, lin_f0, 1);
        
        mzd_write_bit(M, row_d0, get_quad_idx(off_d0+i, off_d0+i, N), 1); // Add Linear D
        
        // D1' += F(D1)
        mzd_write_bit(M, row_d1, quad_f1, 1);
        mzd_write_bit(M, row_d1, lin_f1, 1);
        mzd_write_bit(M, row_d1, get_quad_idx(off_d1+i, off_d1+i, N), 1); // Add Linear D

        // --- 2. Key Addition (Masked on c0, c1) ---
        
        // Add Kd0 to D0 branch
        int k0_bit = (kd0 >> i) & 1;
        if (k0_bit) {
            mzd_write_bit(M, row_d0, idx_c0c0, 1);
            mzd_write_bit(M, row_d0, idx_c1c1, 1);
            mzd_write_bit(M, row_d0, idx_c0c1, 1);
        }

        // Add Kd1 to D1 branch
        int k1_bit = (kd1 >> i) & 1;
        if (k1_bit) {
            mzd_write_bit(M, row_d1, idx_c0c0, 1);
            mzd_write_bit(M, row_d1, idx_c1c1, 1);
            mzd_write_bit(M, row_d1, idx_c0c1, 1);
        }
    }
    return M;
}

// --- Main Generator ---

void wb_generate_tables(WB_Context *ctx, const uint8_t *master_key) {
    printf("[WB-Gen] Tables with Dummy Lines (d=%d)...\n", ctx->d);
    uint16_t *rks = simon32_key_schedule(master_key, ctx->T);

    for (int r = 0; r < ctx->T; r++) {
        mzd_t *M_final = ctx->round_matrices[r];
        
        mzd_t *M_const = wb_gen_layer_constants(ctx);
        mzd_t *M_logic = wb_gen_layer_simon_lrmix(ctx);
        mzd_t *M_key   = wb_gen_layer_key(ctx, rks[r]);
        mzd_t *M_dummy = wb_gen_layer_dummy(ctx, r); // New Layer!

        mzd_add(M_final, M_final, M_const);
        mzd_add(M_final, M_final, M_logic);
        mzd_add(M_final, M_final, M_key);
        mzd_add(M_final, M_final, M_dummy);

        mzd_free(M_const);
        mzd_free(M_logic);
        mzd_free(M_key);
        mzd_free(M_dummy);
    }
    free(rks);
}

// --- Encrypt ---

void wb_encrypt(WB_Context *ctx, const uint8_t *in, uint8_t *out) {
    int N = ctx->N;
    int N_quad = ctx->N_quad;
    int n = ctx->n;
    int d = ctx->d;

    // 1. Input Loading
    uint16_t L = ((uint16_t)in[0] << 8) | (uint16_t)in[1];
    uint16_t R = ((uint16_t)in[2] << 8) | (uint16_t)in[3];

    uint16_t U = L;
    uint16_t V = L ^ R;

    mzd_t *x = mzd_init(N, 1);
    
    // Constants
    int c0 = rand() & 1;
    int c1 = 1 ^ c0;
    mzd_write_bit(x, 0, 0, c0);
    mzd_write_bit(x, 1, 0, c1);

    // Real State
    for(int i=0; i<n; i++) {
        mzd_write_bit(x, 2 + i, 0, (U >> i) & 1);
        mzd_write_bit(x, 2 + n + i, 0, (V >> i) & 1);
    }
    
    // Dummy State Initialization
    // We initialize D0 = D1 = Random (or 0)
    // To verify QAM, they must start equal.
    uint16_t initial_dummy = rand() & 0xFFFF;
    int off_d0 = 2 + 4 * n;
    int off_d1 = 2 + 4 * n + d;
    
    for(int i=0; i<d; i++) {
        int bit = (initial_dummy >> i) & 1;
        mzd_write_bit(x, off_d0 + i, 0, bit); // D0
        mzd_write_bit(x, off_d1 + i, 0, bit); // D1
    }

    // Rounds
    mzd_t *x_quad = mzd_init(N_quad, 1);
    mzd_t *y = mzd_init(N, 1);

    for (int r = 0; r < ctx->T; r++) {
        int q_idx = 0;
        for (int i = 0; i < N; i++) {
            int bit_i = mzd_read_bit(x, i, 0);
            // Optimization: if bit_i is 0, skip inner loop writes? 
            // For logic simplicity, keep full loop.
            for (int j = i; j < N; j++) {
                int bit_j = mzd_read_bit(x, j, 0);
                mzd_write_bit(x_quad, q_idx, 0, bit_i & bit_j);
                q_idx++;
            }
        }
        mzd_mul(y, ctx->round_matrices[r], x_quad, 0);
        mzd_copy(x, y);
    }

    // Output
    U = 0; V = 0;
    for(int i=0; i<n; i++) {
        if(mzd_read_bit(x, 2 + i, 0)) U |= (1 << i);
        if(mzd_read_bit(x, 2 + n + i, 0)) V |= (1 << i);
    }

    L = U;
    R = U ^ V;

    out[0] = (L >> 8) & 0xFF;
    out[1] = L & 0xFF;
    out[2] = (R >> 8) & 0xFF;
    out[3] = R & 0xFF;

    mzd_free(x);
    mzd_free(x_quad);
    mzd_free(y);
}