#include <stdio.h>
#include "wbsimon.h"

void print_hex(const char *label, const uint8_t *data, int len) {
    printf("%s: ", label);
    for (int i = 0; i < len; i++) {
        printf("%02X ", data[i]);
    }
    printf("\n");
}

int main() {
    WB_Context ctx;
    
    printf("=== SIMON 32/64 White-Box Implementation ===\n");

    // Standard Test Vector
    uint8_t key[8] = {0x19, 0x18, 0x11, 0x10, 0x09, 0x08, 0x01, 0x00};
    uint8_t pt[4] = {0x65, 0x65, 0x68, 0x77};
    uint8_t ct[4] = {0};

    // 1. Initialize Context
    wb_ctx_init(&ctx, SIMON_32_64);

    // 2. Generate White-Box Matrices (Offline Phase)
    wb_generate_matrix(&ctx, key);

    // 3. Perform Encryption (Online Phase)
    // Step A: Map Input (External)
    mzd_t *state = wb_map_input(&ctx, pt);

    // Step B: Execute Core Rounds (Main WB Execution)
    wb_execute_rounds(&ctx, state);

    // Step C: Map Output (External)
    wb_map_output(&ctx, state, ct);
    
    // Verify
    print_hex("Ciphertext", ct, 4);
    if (ct[0] == 0xC6 && ct[1] == 0x9B && ct[2] == 0xE9 && ct[3] == 0xBB) {
        printf("[SUCCESS] Ciphertext matches standard test vector.\n");
    } else {
        printf("[FAILURE] Ciphertext mismatch.\n");
    }

    // Cleanup
    wb_ctx_free(&ctx);
    return 0;
}