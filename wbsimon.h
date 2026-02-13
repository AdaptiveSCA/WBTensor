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
    
    // External Encodings (Linear N x N matrices)
    // Used to bridge the external world (plaintext/ciphertext) with the white-box internal state.
    mzd_t *ext_encoding_in;      // Applied to Plaintext -> State
    mzd_t *ext_encoding_out_inv; // Applied to State -> Ciphertext (Decoding)

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
 * STEP 1: Input Mapping
 * Converts the 4-byte plaintext into the internal White-Box state vector.
 * This function handles:
 * 1. Loading bytes into words.
 * 2. Injecting randomness (Masking shares, Dummy values, Constants).
 * 3. Applying the External Input Linear Encoding (P_in).
 * * @param ctx The white-box context.
 * @param in The 4-byte plaintext input.
 * @return A pointer to the allocated M4RI state vector (N x 1).
 */
mzd_t* wb_map_input(WB_Context *ctx, const uint8_t *in);

/**
 * STEP 2: Core Execution
 * Executes the white-box round functions on the state vector.
 * This performs the quadratic evaluation and matrix multiplication for all rounds.
 * The state vector is updated in-place.
 * * @param ctx The white-box context.
 * @param state The state vector returned by wb_map_input.
 */
void wb_execute_rounds(WB_Context *ctx, mzd_t *state);

/**
 * STEP 3: Output Mapping
 * Converts the final internal state vector back to ciphertext bytes.
 * This function handles:
 * 1. Applying the External Output Linear Decoding (P_out^-1).
 * 2. Recombining the secret shares (U0^U1, V0^V1).
 * 3. Serializing the result to bytes.
 * 4. Frees the state vector memory.
 * * @param ctx The white-box context.
 * @param state The state vector (will be freed by this function).
 * @param out Buffer to store the 4-byte ciphertext.
 */
void wb_map_output(WB_Context *ctx, mzd_t *state, uint8_t *out);

#endif // WB_SIMON_H