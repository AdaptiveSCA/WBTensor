#ifndef WB_SIMON_H
#define WB_SIMON_H

#include <stdint.h>
#include <m4ri/m4ri.h>

/**
 * SIMON Mode Configuration
 * Currently supports SIMON 32/64.
 */
typedef enum {
    SIMON_32_64
} SimonMode;

/**
 * White-Box Context Structure
 * Holds the configuration and the generated round matrices.
 */
typedef struct {
    SimonMode mode;
    int n;          // Word size (16 for SIMON 32/64)
    int T;          // Number of rounds (32 for SIMON 32/64)
    int d;          // Number of dummy bits per branch (equals n)
    
    // Dimension Definitions
    // N: Linear state size. 
    // Layout: [c0, c1, U_share0, V_share0, U_share1, V_share1, Dummy0, Dummy1]
    int N;          
    
    // N_quad: Quadratic state size (tensor product dimension)
    int N_quad;     

    // Round Matrices
    // An array of pointers to the compressed quadratic matrices for each round.
    mzd_t **round_matrices; 

    // QAM Debugging Constants
    // Stores the expected XOR difference (C) between Dummy0 and Dummy1 for each round.
    uint16_t *dummy_consts; 

} WB_Context;

// --- Context Management ---

/**
 * Initializes the white-box context, calculating dimensions and allocating memory.
 * @param ctx Pointer to the context structure.
 * @param mode The SIMON algorithm mode.
 */
void wb_ctx_init(WB_Context *ctx, SimonMode mode);

/**
 * Frees all memory associated with the white-box context.
 * @param ctx Pointer to the context structure.
 */
void wb_ctx_free(WB_Context *ctx);

// --- Core Generation ---

/**
 * Generates the complete white-box matrix for all rounds.
 * This function orchestrates the generation of logic layers (Constants, Logic, Key, Dummy, Masking)
 * and fuses them into a single quadratic matrix per round.
 * * @param ctx Pointer to the context structure.
 * @param master_key The 64-bit master key (byte array).
 */
void wb_generate_matrix(WB_Context *ctx, const uint8_t *master_key);

// --- Execution ---

/**
 * Performs white-box encryption on a plaintext block.
 * Handles input encoding (splitting into shares), state randomization,
 * round execution, and output decoding.
 * * @param ctx Pointer to the initialized context with generated matrices.
 * @param in Input plaintext (4 bytes).
 * @param out Output ciphertext (4 bytes).
 */
void wb_encrypt(WB_Context *ctx, const uint8_t *in, uint8_t *out);

#endif // WB_SIMON_H