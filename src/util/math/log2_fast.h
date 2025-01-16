#include <cstdint>


/**
 * Fast approximation of log2 (doesn't work when x<2)
 * @param x
 * @return log2 of x
 */
static inline float log2_approximate(float x)
{
    union { float f; uint32_t i; } z = { x };
    float log_2 = ((z.i >> 23) & 255) - 128;
    z.i &= ~(255 << 23);
    z.i += 127 << 23;
    log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
    return log_2;
}