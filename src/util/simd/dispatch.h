/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

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

#ifdef WITH_NEON
#define HAVE_NEON(x) x
#else
#define HAVE_NEON(x)
#endif

#if defined(WITH_NEON) | defined(WITH_AVX2) | defined(WITH_SSE4_1)
#define HAVE_SIMD(x) x
#else
#define HAVE_SIMD(x)
#endif

#if ARCH_ID == 0

#define DISPATCH_0V(name)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { void name(); })\
HAVE_AVX2(namespace ARCH_AVX2 { void name(); })\
HAVE_NEON(namespace ARCH_NEON { void name(); })\
void name() {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: ARCH_NEON::name(); break;)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: ARCH_AVX2::name(); break;)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: ARCH_SSE4_1::name(); break;)\
HAVE_SIMD(default:)\
ARCH_GENERIC::name();\
HAVE_SIMD(})\
}

#define DISPATCH_1V(name, t1, n1)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { void name(t1 n1); })\
HAVE_AVX2(namespace ARCH_AVX2 { void name(t1 n1); })\
HAVE_NEON(namespace ARCH_NEON { void name(t1 n1); })\
void name(t1 n1) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: ARCH_NEON::name(n1); break;)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: ARCH_AVX2::name(n1); break;)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: ARCH_SSE4_1::name(n1); break;)\
HAVE_SIMD(default:)\
ARCH_GENERIC::name(n1);\
HAVE_SIMD(})\
}

#define DISPATCH_1(ret, name, t1, n1)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1); })\
HAVE_NEON(namespace ARCH_NEON { ret name(t1 n1); })\
ret name(t1 n1) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: return ARCH_NEON::name(n1);)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1);)\
HAVE_SIMD(default:)\
return ARCH_GENERIC::name(n1);\
HAVE_SIMD(})\
}

#define DISPATCH_2(ret, name, t1, n1, t2, n2)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2); })\
HAVE_NEON(namespace ARCH_NEON { ret name(t1 n1, t2 n2); })\
ret name(t1 n1, t2 n2) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: return ARCH_NEON::name(n1, n2);)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2);)\
HAVE_SIMD(default:)\
return ARCH_GENERIC::name(n1, n2);\
HAVE_SIMD(})\
}

#define DISPATCH_3(ret, name, t1, n1, t2, n2, t3, n3)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3); })\
HAVE_NEON(namespace ARCH_NEON { ret name(t1 n1, t2 n2, t3 n3); })\
ret name(t1 n1, t2 n2, t3 n3) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: return ARCH_NEON::name(n1, n2, n3);)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3);)\
HAVE_SIMD(default:)\
return ARCH_GENERIC::name(n1, n2, n3);\
HAVE_SIMD(})\
}

#define DISPATCH_3V(name, t1, n1, t2, n2, t3, n3)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { void name(t1 n1, t2 n2, t3 n3); })\
HAVE_AVX2(namespace ARCH_AVX2 { void name(t1 n1, t2 n2, t3 n3); })\
HAVE_NEON(namespace ARCH_NEON { void name(t1 n1, t2 n2, t3 n3); })\
void name(t1 n1, t2 n2, t3 n3) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: ARCH_NEON::name(n1, n2, n3); break;)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: ARCH_AVX2::name(n1, n2, n3); break;)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: ARCH_SSE4_1::name(n1, n2, n3); break;)\
HAVE_SIMD(default:)\
ARCH_GENERIC::name(n1, n2, n3);\
HAVE_SIMD(})\
}

#define DISPATCH_4(ret, name, t1, n1, t2, n2, t3, n3, t4, n4)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3, t4 n4); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3, t4 n4); })\
HAVE_NEON(namespace ARCH_NEON { ret name(t1 n1, t2 n2, t3 n3, t4 n4); })\
ret name(t1 n1, t2 n2, t3 n3, t4 n4) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: return ARCH_NEON::name(n1, n2, n3, n4);)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3, n4);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3, n4);)\
HAVE_SIMD(default:)\
return ARCH_GENERIC::name(n1, n2, n3, n4);\
HAVE_SIMD(})\
}

#define DISPATCH_5(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5); })\
HAVE_NEON(namespace ARCH_NEON { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5); })\
ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: return ARCH_NEON::name(n1, n2, n3, n4, n5);)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3, n4, n5);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3, n4, n5);)\
HAVE_SIMD(default:)\
return ARCH_GENERIC::name(n1, n2, n3, n4, n5);\
HAVE_SIMD(})\
}

#define DISPATCH_6(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6); })\
HAVE_NEON(namespace ARCH_NEON { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6); })\
ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: return ARCH_NEON::name(n1, n2, n3, n4, n5, n6);)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3, n4, n5, n6);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3, n4, n5, n6);)\
HAVE_SIMD(default:)\
return ARCH_GENERIC::name(n1, n2, n3, n4, n5, n6);\
HAVE_SIMD(})\
}

#define DISPATCH_6V(name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6); })\
HAVE_AVX2(namespace ARCH_AVX2 { void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6); })\
HAVE_NEON(namespace ARCH_NEON { void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6); })\
void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: ARCH_NEON::name(n1, n2, n3, n4, n5, n6); break;)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: ARCH_AVX2::name(n1, n2, n3, n4, n5, n6); break;)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: ARCH_SSE4_1::name(n1, n2, n3, n4, n5, n6); break;)\
HAVE_SIMD(default:)\
ARCH_GENERIC::name(n1, n2, n3, n4, n5, n6);\
HAVE_SIMD(})\
}

#define DISPATCH_7(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7); })\
HAVE_NEON(namespace ARCH_NEON { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7); })\
ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: return ARCH_NEON::name(n1, n2, n3, n4, n5, n6, n7);)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3, n4, n5, n6, n7);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3, n4, n5, n6, n7);)\
HAVE_SIMD(default:)\
return ARCH_GENERIC::name(n1, n2, n3, n4, n5, n6, n7);\
HAVE_SIMD(})\
}

#define DISPATCH_7V(name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7); })\
HAVE_AVX2(namespace ARCH_AVX2 { void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7); })\
HAVE_NEON(namespace ARCH_NEON { void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7); })\
void name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: ARCH_NEON::name(n1, n2, n3, n4, n5, n6, n7); break;)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: ARCH_AVX2::name(n1, n2, n3, n4, n5, n6, n7); break;)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: ARCH_SSE4_1::name(n1, n2, n3, n4, n5, n6, n7); break;)\
HAVE_SIMD(default:)\
ARCH_GENERIC::name(n1, n2, n3, n4, n5, n6, n7);\
HAVE_SIMD(})\
}

#define DISPATCH_8(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7, t8, n8)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7, t8 n8); })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7, t8 n8); })\
HAVE_NEON(namespace ARCH_NEON { ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7, t8 n8); })\
ret name(t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7, t8 n8) {\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: return ARCH_NEON::name(n1, n2, n3, n4, n5, n6, n7, n8);)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name(n1, n2, n3, n4, n5, n6, n7, n8);)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name(n1, n2, n3, n4, n5, n6, n7, n8);)\
HAVE_SIMD(default:)\
return ARCH_GENERIC::name(n1, n2, n3, n4, n5, n6, n7, n8);\
HAVE_SIMD(})\
}

#else

#define DISPATCH_0V(name)
#define DISPATCH_1V(name, t1, n1)
#define DISPATCH_1(ret, name, t1, n1)
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