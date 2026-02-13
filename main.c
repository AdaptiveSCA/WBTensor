#include <stdio.h>
#include "wbsimon.h"

// Utility to print hex output
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
    printf("Features: Tensor-based, L-R Mixing, Masking (Splitting), Dummies, QAM\n");

    // Standard Test Vector for SIMON 32/64
    // Key: 19 18 11 10 09 08 01 00
    uint8_t key[8] = {0x19, 0x18, 0x11, 0x10, 0x09, 0x08, 0x01, 0x00};
    // Plaintext: 65 65 68 77
    uint8_t pt[4] = {0x65, 0x65, 0x68, 0x77};
    // Expected Ciphertext: c6 9b e9 bb
    uint8_t ct[4] = {0};

    // Initialize context
    wb_ctx_init(&ctx, SIMON_32_64);

    // Generate White-Box Matrices
    wb_generate_matrix(&ctx, key);

    // Perform Encryption
    wb_encrypt(&ctx, pt, ct);
    
    print_hex("Ciphertext", ct, 4);
    
    // Verification
    if (ct[0] == 0xC6 && ct[1] == 0x9B && ct[2] == 0xE9 && ct[3] == 0xBB) {
        printf("[SUCCESS] Ciphertext matches standard test vector.\n");
    } else {
        printf("[FAILURE] Ciphertext mismatch.\n");
    }

    // Cleanup
    wb_ctx_free(&ctx);
    return 0;
}