#pragma once

#ifdef WITH_AVX2
#define HAVE_AVX2(x) x
#else
#define HAVE_AVX2(x)
#endif

#ifdef WITH_SSE4_1
#define HAVE_SSE4_1(x) x
#else
#define HAVE_SSE4_1(x)
#endif

#if ARCH_ID == 0

#define DISPATCH_0V(name)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { void name(); })\
HAVE_AVX2(namespace ARCH_AVX2 { void name(); })\
void name() {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: ARCH_AVX2::name(); break;)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: ARCH_SSE4_1::name(); break;)\
HAVE_SSE4_1(default:)\
ARCH_GENERIC::name();\
HAVE_SSE4_1(})\
}

#define DISPATCH_2(ret, name, t1, n1, t2, n2)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2); })\
ret name(t1 n1, t2 n2) {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2);)\
HAVE_SSE4_1(default:)\
return ARCH_GENERIC::name(n1, n2);\
HAVE_SSE4_1(})\
}

#define DISPATCH_3(ret, name, t1, n1, t2, n2, t3, n3)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3); })\
ret name(t1 n1, t2 n2, t3 n3) {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3);)\
HAVE_SSE4_1(default:)\
return ARCH_GENERIC::name(n1, n2, n3);\
HAVE_SSE4_1(})\
}

#define DISPATCH_3V(name, t1, n1, t2, n2, t3, n3)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { void name(t1 n1, t2 n2, t3 n3); })\
HAVE_AVX2(namespace ARCH_AVX2 { void name(t1 n1, t2 n2, t3 n3); })\
void name(t1 n1, t2 n2, t3 n3) {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: ARCH_AVX2::name(n1, n2, n3); break;)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: ARCH_SSE4_1::name(n1, n2, n3); break;)\
HAVE_SSE4_1(default:)\
ARCH_GENERIC::name(n1, n2, n3);\
HAVE_SSE4_1(})\
}

#define DISPATCH_4(ret, name, t1, n1, t2, n2, t3, n3, t4, n4)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3, t4 n4); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3, t4 n4); })\
ret name(t1 n1, t2 n2, t3 n3, t4 n4) {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3, n4);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3, n4);)\
HAVE_SSE4_1(default:)\
return ARCH_GENERIC::name(n1, n2, n3, n4);\
HAVE_SSE4_1(})\
}

#define DISPATCH_5(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5); })\
ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5) {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3, n4, n5);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3, n4, n5);)\
HAVE_SSE4_1(default:)\
return ARCH_GENERIC::name(n1, n2, n3, n4, n5);\
HAVE_SSE4_1(})\
}

#define DISPATCH_6(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6); })\
ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6) {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3, n4, n5, n6);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3, n4, n5, n6);)\
HAVE_SSE4_1(default:)\
return ARCH_GENERIC::name(n1, n2, n3, n4, n5, n6);\
HAVE_SSE4_1(})\
}

#define DISPATCH_6V(name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6); })\
HAVE_AVX2(namespace ARCH_AVX2 { void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6); })\
void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6) {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: ARCH_AVX2::name(n1, n2, n3, n4, n5, n6); break;)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: ARCH_SSE4_1::name(n1, n2, n3, n4, n5, n6); break;)\
HAVE_SSE4_1(default:)\
ARCH_GENERIC::name(n1, n2, n3, n4, n5, n6);\
HAVE_SSE4_1(})\
}

#define DISPATCH_7(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7); })\
ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7) {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3, n4, n5, n6, n7);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3, n4, n5, n6, n7);)\
HAVE_SSE4_1(default:)\
return ARCH_GENERIC::name(n1, n2, n3, n4, n5, n6, n7);\
HAVE_SSE4_1(})\
}

#define DISPATCH_7V(name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7); })\
HAVE_AVX2(namespace ARCH_AVX2 { void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7); })\
void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7) {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: ARCH_AVX2::name(n1, n2, n3, n4, n5, n6, n7); break;)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: ARCH_SSE4_1::name(n1, n2, n3, n4, n5, n6, n7); break;)\
HAVE_SSE4_1(default:)\
ARCH_GENERIC::name(n1, n2, n3, n4, n5, n6, n7);\
HAVE_SSE4_1(})\
}

#define DISPATCH_8(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7, t8, n8)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7, t8 n8); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7, t8 n8); })\
ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7, t8 n8) {\
HAVE_SSE4_1(switch(::SIMD::arch()) {)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3, n4, n5, n6, n7, n8);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3, n4, n5, n6, n7, n8);)\
HAVE_SSE4_1(default:)\
return ARCH_GENERIC::name(n1, n2, n3, n4, n5, n6, n7, n8);\
HAVE_SSE4_1(})\
}

#else

#define DISPATCH_0V(name)
#define DISPATCH_2(ret, name, t1, n1, t2, n2)
#define DISPATCH_3(ret, name, t1, n1, t2, n2, t3, n3)
#define DISPATCH_4(ret, name, t1, n1, t2, n2, t3, n3, t4, n4)
#define DISPATCH_5(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5)
#define DISPATCH_6(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6)
#define DISPATCH_7(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7)
#define DISPATCH_3V(name, t1, n1, t2, n2, t3, n3)
#define DISPATCH_6V(name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6)
#define DISPATCH_7V(name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7)
#define DISPATCH_8(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7, t8, n8)

#endif