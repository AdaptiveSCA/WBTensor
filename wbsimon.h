#ifndef WB_SIMON_H
#define WB_SIMON_H

#include <stdint.h>
#include <m4ri/m4ri.h>

typedef enum { SIMON_32_64 } SimonMode;

typedef struct {
    SimonMode mode;
    int n;          // 16
    int T;          // 32
    int d;          // Dummy size (Set to n=16)
    
    // Layout: [c0, c1, U(n), V(n), Empty(2n), Dummy0(d), Dummy1(d)]
    int N;          
    int N_quad;     
    mzd_t **round_matrices; 

    // QAM Constants (Stored for debugging or future QAM check verification)
    // C[r] = K_dummy0[r] ^ K_dummy1[r]
    uint16_t *dummy_consts; 

} WB_Context;

void wb_ctx_init(WB_Context *ctx, SimonMode mode); // d is implied as n
void wb_ctx_free(WB_Context *ctx);
void wb_generate_tables(WB_Context *ctx, const uint8_t *master_key);
void wb_encrypt(WB_Context *ctx, const uint8_t *in, uint8_t *out);

#endif // WB_SIMON_H