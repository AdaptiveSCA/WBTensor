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
    
    // Key: 19 18 11 10 09 08 01 00
    uint8_t key[8] = {0x19, 0x18, 0x11, 0x10, 0x09, 0x08, 0x01, 0x00};
    uint8_t pt[4] = {0x65, 0x65, 0x68, 0x77};
    uint8_t ct[4] = {0};

    // Initialize with d=n 
    wb_ctx_init(&ctx, SIMON_32_64);

    wb_generate_tables(&ctx, key);

    wb_encrypt(&ctx, pt, ct);
    
    print_hex("Ciphertext", ct, 4);
    printf("Exp: C6 9B E9 BB\n");

    wb_ctx_free(&ctx);
    return 0;
}