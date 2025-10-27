/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <string_view>
#include <utility>
#include "util/io/mmap.h"
#include "phr.h"

using std::string;
using std::runtime_error;
using std::out_of_range;

// --------------------- Util ---------------------

static inline void die(const char* msg) {
    std::fprintf(stderr, "fatal: %s\n", msg);
    std::exit(1);
}

static inline uint32_t be32(const uint8_t* p) {
    return (uint32_t)p[0] << 24 | (uint32_t)p[1] << 16 | (uint32_t)p[2] << 8 | (uint32_t)p[3];
}
static inline uint64_t le64(const uint8_t* p) {
    return (uint64_t)p[0] | (uint64_t)p[1] << 8 | (uint64_t)p[2] << 16 | (uint64_t)p[3] << 24 |
        (uint64_t)p[4] << 32 | (uint64_t)p[5] << 40 | (uint64_t)p[6] << 48 | (uint64_t)p[7] << 56;
}

struct Span {
    const uint8_t* p{ nullptr };
    const uint8_t* end{ nullptr };
    bool empty() const { return p >= end; }
    size_t size() const { return (size_t)(end - p); }
    void debug_print(OId i) const {
        size_t sn = std::min<size_t>(32, (size_t)(end - p));
        std::fprintf(stderr, "dbg i=%zu first-bytes:", i);
        for (size_t k = 0; k < sn; ++k) std::fprintf(stderr, " %02X", p[k]);
        std::fprintf(stderr, "\n");
    }
};

enum TagClass : uint8_t { TC_Universal = 0, TC_Application = 1, TC_Context = 2, TC_Private = 3 };

struct TLV {
    TagClass cls;
    bool constructed;
    uint8_t tagnum;      // assume < 31
    bool indefinite;
    size_t len;          // valid iff !indefinite
    const uint8_t* val;  // value start
};

static inline bool parse_tag(Span& s, TLV& out) {
    if (s.p >= s.end) return false;
    uint8_t b = *s.p++;
    out.cls = (TagClass)((b >> 6) & 0x3);
    out.constructed = (b & 0x20) != 0;
    uint8_t tn = (uint8_t)(b & 0x1F);
    if (tn == 0x1F) return false;
    out.tagnum = tn;
    return true;
}
static inline bool parse_len(Span& s, TLV& out) {
    if (s.p >= s.end) return false;
    uint8_t b = *s.p++;
    if (b == 0x80) { out.indefinite = true; out.len = 0; out.val = s.p; return true; }
    out.indefinite = false;
    if ((b & 0x80) == 0) { out.len = b; }
    else {
        int n = b & 0x7F; if (n <= 0 || n > 8) return false;
        if ((size_t)(s.end - s.p) < (size_t)n) return false;
        size_t L = 0; for (int i = 0; i < n; ++i) L = (L << 8) | s.p[i];
        s.p += n; out.len = L;
    }
    out.val = s.p;
    if ((size_t)(s.end - s.p) < out.len) return false;
    return true;
}
static inline bool parse_tlv_header(Span& s, TLV& out) { return parse_tag(s, out) && parse_len(s, out); }
static inline bool is_eoc(const Span& s) { return (s.end - s.p) >= 2 && s.p[0] == 0x00 && s.p[1] == 0x00; }
static inline Span limit_span(const TLV& t) { Span sub{ t.val, t.val + t.len }; return sub; }
static bool skip_value(Span& s, const TLV& t);
static bool skip_constructed(Span& s, const TLV& t) {
    if (!t.constructed) return false;
    if (!t.indefinite) {
        Span sub = limit_span(t);
        while (!sub.empty()) {
            TLV c; if (!parse_tlv_header(sub, c)) return false;
            if (!skip_value(sub, c)) return false;
        }
        s.p = t.val + t.len; return true;
    }
    for (;;) {
        if (is_eoc(s)) { s.p += 2; break; }
        TLV c; if (!parse_tlv_header(s, c)) return false;
        if (!skip_value(s, c)) return false;
    }
    return true;
}
static bool skip_value(Span& s, const TLV& t) {
    if (!t.constructed) { s.p = t.val + (t.indefinite ? 0 : t.len); return true; }
    return skip_constructed(s, t);
}

static inline bool is_string_primitive(const TLV& t) {
    if (t.cls != TC_Universal || t.constructed) return false;
    switch (t.tagnum) { // UTF8(12), Printable(19), IA5(22), Visible(26), General(27)
    case 0x0C: case 0x13: case 0x16: case 0x1A: case 0x1B: return true;
    default: return false;
    }
}
static inline bool is_integer_primitive(const TLV& t) {
    return (t.cls == TC_Universal && !t.constructed && t.tagnum == 0x02);
}
static bool parse_integer(const TLV& t, uint64_t& out) {
    if (t.indefinite || t.constructed) return false;
    const uint8_t* v = t.val; size_t n = t.len;
    if (n == 0) return false; if ((v[0] & 0x80) != 0) return false;
    uint64_t x = 0; for (size_t i = 0; i < n; ++i) x = (x << 8) | v[i]; out = x; return true;
}

// Context helpers: [n] EXPLICIT wrappers are common in .phr
static inline bool peek_ctx_constructed(Span s, uint8_t tagnum, TLV& out) {
    if (!parse_tlv_header(s, out)) return false;
    return (out.cls == TC_Context && out.constructed && out.tagnum == tagnum);
}

static bool read_ctx_string(Span& s, uint8_t tagnum,
    const uint8_t*& ptr, size_t& len) {
    TLV w;
    if (!peek_ctx_constructed(s, tagnum, w)) return false;

    if (!w.indefinite) {
        // Definite wrapper: read from a local span, then advance s past wrapper
        Span body = limit_span(w);
        TLV inner;
        if (!parse_tlv_header(body, inner) || !is_string_primitive(inner)) return false;
        ptr = inner.val; len = inner.len;
        s.p = w.val + w.len;
        return true;
    }
    else {
        // IMPORTANT: actually consume the wrapper header on 's'
        TLV w2; if (!parse_tlv_header(s, w2)) return false; // s -> wrapper body
        TLV inner;
        if (!parse_tlv_header(s, inner) || !is_string_primitive(inner)) return false;
        ptr = inner.val; len = inner.len;
        s.p = inner.val + inner.len;          // after string bytes
        if (is_eoc(s)) s.p += 2;              // EOC for wrapper body
        if (is_eoc(s)) s.p += 2;              // (some encoders double-EOC)
        return true;
    }
}

static bool read_ctx_sequence(Span& s, uint8_t tagnum, Span& out_body, bool& indefinite) {
    TLV w;
    if (!peek_ctx_constructed(s, tagnum, w)) return false;

    // Consume wrapper header on a working copy to inspect the inner SEQUENCE…
    Span tmp = s;
    parse_tlv_header(tmp, w);            // tmp now at wrapper body
    TLV seq;
    if (!parse_tlv_header(tmp, seq)) return false;
    if (!(seq.cls == TC_Universal && seq.constructed && seq.tagnum == 0x10)) return false;

    if (!seq.indefinite) {
        out_body = limit_span(seq);
        // Now advance real 's' over the whole wrapper
        s.p = w.val + w.len;
        indefinite = false;
        return true;
    }
    else {
        // IMPORTANT: on the real 's', consume the wrapper header before skipping it
        TLV w2; if (!parse_tlv_header(s, w2)) return false; // s -> wrapper body
        // Provide a body span that starts right after the SEQUENCE header
        TLV seq2; if (!parse_tlv_header(s, seq2)) return false;
        if (!(seq2.cls == TC_Universal && seq2.constructed && seq2.tagnum == 0x10)) return false;
        out_body = s;                       // body begins after seq header
        // Advance 's' to the end of the wrapper efficiently
        if (!skip_constructed(s, w2)) return false;
        indefinite = true;
        return true;
    }
}

static bool read_ctx_integer(Span& s, uint8_t tagnum, uint64_t& value) {
    TLV w;
    if (!peek_ctx_constructed(s, tagnum, w)) return false;

    if (!w.indefinite) {
        Span body = limit_span(w);
        TLV inner;
        if (!parse_tlv_header(body, inner) || !is_integer_primitive(inner)) return false;
        if (!parse_integer(inner, value)) return false;
        s.p = w.val + w.len;
        return true;
    }
    else {
        TLV w2; if (!parse_tlv_header(s, w2)) return false; // consume wrapper header
        TLV inner;
        if (!parse_tlv_header(s, inner) || !is_integer_primitive(inner)) return false;
        if (!parse_integer(inner, value)) return false;
        s.p = inner.val + inner.len;
        if (is_eoc(s)) s.p += 2;
        if (is_eoc(s)) s.p += 2;
        return true;
    }
}

// Best-effort Seq-id stringify (lightweight)
static void seqid_choice_best_effort(Span sub, std::string& out_text) {
    std::string_view s1, s2; uint64_t intv = 0; bool have_int = false;
    while (!sub.empty()) {
        if (is_eoc(sub)) { sub.p += 2; break; }
        TLV c; if (!parse_tlv_header(sub, c)) break;
        if (is_string_primitive(c)) {
            if (c.len > 0 && s1.empty()) s1 = std::string_view((const char*)c.val, c.len);
            else if (c.len > 0 && s2.empty()) s2 = std::string_view((const char*)c.val, c.len);
            sub.p = c.val + c.len; continue;
        }
        if (is_integer_primitive(c) && !have_int) { if (parse_integer(c, intv)) have_int = true; sub.p = c.val + c.len; continue; }
        if (c.constructed) {
            if (!c.indefinite) { Span deeper = limit_span(c); seqid_choice_best_effort(deeper, out_text); sub.p = c.val + c.len; }
            else { Span deeper = sub; seqid_choice_best_effort(deeper, out_text); if (!skip_constructed(sub, c)) break; }
            if (!out_text.empty()) break;
            continue;
        }
        sub.p = c.val + (c.indefinite ? 0 : c.len);
    }
    if (out_text.empty()) {
        if (!s1.empty() && !s2.empty()) { out_text.append(s1.data(), s1.size()); out_text.push_back('|'); out_text.append(s2.data(), s2.size()); }
        else if (!s1.empty()) { out_text.append(s1.data(), s1.size()); }
        else if (have_int) { char buf[32]; std::snprintf(buf, sizeof(buf), "%llu", (unsigned long long)intv); out_text.append(buf); }
    }
}
static void parse_seqid_list(Span body, std::vector<std::string>& ids, size_t max_ids = 4) {
    while (!body.empty()) {
        if (is_eoc(body)) { body.p += 2; break; }
        TLV ch; if (!parse_tlv_header(body, ch)) break;
        if (ch.cls == TC_Context && ch.constructed) {
            if (!ch.indefinite) {
                Span sub = limit_span(ch); std::string out; seqid_choice_best_effort(sub, out);
                if (!out.empty()) { ids.emplace_back(std::move(out)); if (ids.size() >= max_ids) { body.p = ch.val + ch.len; break; } }
                body.p = ch.val + ch.len;
            }
            else {
                Span sub = body; std::string out; seqid_choice_best_effort(sub, out);
                if (!skip_constructed(body, ch)) break;
                if (!out.empty()) { ids.emplace_back(std::move(out)); if (ids.size() >= max_ids) break; }
            }
        }
        else {
            if (!skip_value(body, ch)) break;
        }
    }
}

// Parse a Blast-def-line whose fields ([0] title, [1] seqid, [2] taxid, …)
// appear directly in 'body'. If we encounter a nested SEQUENCE, we recurse.
static bool parse_blast_def_line_fields_inplace(
    Span& body,
    std::vector<std::pair<const uint8_t*, size_t>>& titles,
    std::vector<std::string>& ids,
    std::vector<uint64_t>& taxids)
{
    for (;;) {
        // Termination for definite/indefinite bodies
        if ((body.end - body.p) == 0) return true;
        if (is_eoc(body)) { body.p += 2; return true; }

        // Try known context-specific fields in any order
        const uint8_t* sp = nullptr; size_t sl = 0;
        if (read_ctx_string(body, 0, sp, sl)) {           // [0] title
            titles.emplace_back(sp, sl);
            continue;
        }
        Span seqids; bool indef = false;
        if (read_ctx_sequence(body, 1, seqids, indef)) {  // [1] seqid
            parse_seqid_list(seqids, ids);
            continue;
        }
        uint64_t tx = 0;
        if (read_ctx_integer(body, 2, tx)) {              // [2] taxid
            taxids.push_back(tx);
            continue;
        }

        // Peek next TLV to decide what to do
        Span peek = body;
        TLV t;
        if (!parse_tlv_header(peek, t)) return false;

        // If we see a universal SEQUENCE, it is an inner container; recurse into it.
        if (t.cls == TC_Universal && t.constructed && t.tagnum == 0x10) {
            if (!t.indefinite) {
                // definite-length inner SEQUENCE
                (void)parse_tlv_header(body, t);              // consume header on real stream
                Span inner = limit_span(t);
                if (!parse_blast_def_line_fields_inplace(inner, titles, ids, taxids)) return false;
                body.p = t.val + t.len;                       // advance past inner
            }
            else {
                // indefinite-length inner SEQUENCE
                TLV t2; if (!parse_tlv_header(body, t2)) return false; // consume header
                Span inner = body;                             // copy starting at first child
                if (!parse_blast_def_line_fields_inplace(inner, titles, ids, taxids)) return false;
                if (!skip_constructed(body, t2)) return false; // jump past the whole inner SEQUENCE
            }
            continue;
        }

        // Anything else: skip generically
        (void)parse_tlv_header(body, t);
        if (!skip_value(body, t)) return false;
    }
}

// One Blast-def-line (SEQUENCE with [0]=title [1]=seqid [2]=taxid)
static bool parse_blast_def_line(Span& s,
    std::vector<std::pair<const uint8_t*, size_t>>& titles,
    std::vector<std::string>& ids,
    std::vector<uint64_t>& taxids) {
    TLV t; if (!parse_tlv_header(s, t)) return false;
    if (!(t.cls == TC_Universal && t.constructed && t.tagnum == 0x10)) return false;

    if (!t.indefinite) {
        Span body = limit_span(t);

        // [0] title (optional)
        const uint8_t* sp = nullptr; size_t sl = 0;
        Span probe = body;
        if (read_ctx_string(probe, 0, sp, sl)) { titles.emplace_back(sp, sl); body = probe; }
        else {
            // fallback: bare string first field
            TLV maybe; Span tmp = body;
            if (parse_tlv_header(tmp, maybe) && is_string_primitive(maybe)) { titles.emplace_back(maybe.val, maybe.len); body.p = maybe.val + maybe.len; }
        }

        // [1] seqid
        Span seqids; bool indef = false;
        if (read_ctx_sequence(body, 1, seqids, indef)) { parse_seqid_list(seqids, ids); }
        else {
            TLV maybe; if (parse_tlv_header(body, maybe) && maybe.cls == TC_Universal && maybe.constructed && maybe.tagnum == 0x10) {
                Span seqids2 = limit_span(maybe); parse_seqid_list(seqids2, ids); body.p = maybe.val + maybe.len;
            }
        }

        // [2] taxid (optional)
        uint64_t tx = 0; if (read_ctx_integer(body, 2, tx)) taxids.push_back(tx);

        // skip rest
        while (!body.empty()) { TLV rest; if (!parse_tlv_header(body, rest)) break; if (!skip_value(body, rest)) break; }

        s.p = t.val + t.len; return true;
    }

    // indefinite body
    for (;;) {
        if (is_eoc(s)) { s.p += 2; break; }
        const uint8_t* sp = nullptr; size_t sl = 0;
        if (read_ctx_string(s, 0, sp, sl)) { titles.emplace_back(sp, sl); continue; }
        Span seqids; bool indef = false;
        if (read_ctx_sequence(s, 1, seqids, indef)) { parse_seqid_list(seqids, ids); continue; }
        uint64_t tx = 0; if (read_ctx_integer(s, 2, tx)) { taxids.push_back(tx); continue; }
        TLV any; if (!parse_tlv_header(s, any)) return false; if (!skip_value(s, any)) return false;
    }
    return true;
}

// Heuristic: many .phr records are [len_be][BER...] — return a Span that skips the prefix if it fits exactly.
static inline bool strip_len_prefix_if_fits(Span& s) {
    if ((size_t)(s.end - s.p) < 8) return false;
    uint32_t L = (uint32_t)s.p[0] << 24 | (uint32_t)s.p[1] << 16 | (uint32_t)s.p[2] << 8 | (uint32_t)s.p[3];
    const uint8_t* body = s.p + 4;
    const uint8_t* tail = body + (size_t)L;
    if (tail == s.end) { s.p = body; return true; }
    return false;
}

// Parse the BODY of a top-level SEQUENCE OF Blast-def-line (i.e., you've already consumed the SEQUENCE header).
static bool parse_defline_set_body(Span& body,
    std::vector<std::pair<const uint8_t*, size_t>>& titles,
    std::vector<std::string>& ids,
    std::vector<uint64_t>& taxids,
    bool indefinite) {
    if (!indefinite) {
        while (!body.empty()) {
            if (!parse_blast_def_line(body, titles, ids, taxids)) return false;
        }
        return true;
    }
    for (;;) {
        if (is_eoc(body)) { body.p += 2; break; }
        if (!parse_blast_def_line(body, titles, ids, taxids)) return false;
    }
    return true;
}

// Best-effort: if the next object is a context-specific constructed [X] that EXPLICIT-wraps a SEQUENCE,
// return a Span to the SEQUENCE'S BODY and advance 's' past the whole wrapper.
static bool unwrap_ctx_sequence(Span& s, Span& out_body, bool& indefinite) {
    TLV w;
    Span probe = s;
    if (!parse_tlv_header(probe, w)) return false;
    if (!(w.cls == TC_Context && w.constructed)) return false;

    TLV seq;
    if (!parse_tlv_header(probe, seq)) return false;
    if (!(seq.cls == TC_Universal && seq.constructed && seq.tagnum == 0x10)) return false;

    if (!seq.indefinite) {
        out_body = limit_span(seq);
        // Advance original span over whole wrapper:
        s.p = w.val + w.len;
        indefinite = false;
        return true;
    }
    else {
        // For indefinite inner SEQUENCE, give caller a body span starting at first child:
        out_body = probe; // after seq header
        // And advance original s by skipping the entire wrapper cheaply:
        if (!skip_constructed(s, w)) return false;
        indefinite = true;
        return true;
    }
}

// Top-level record: SEQUENCE whose contents may be the fields directly, or an inner SEQUENCE.
// We parse fields "in place" and recurse into any nested SEQUENCE nodes we encounter.
bool Phr::parse_record(OId i,
    std::vector<std::pair<const uint8_t*, size_t>>& titles,
    std::vector<std::string>& ids,
    std::vector<uint64_t>& taxids)
{
    uint64_t off = idx_.hdr_offsets[i];
    uint64_t off2 = idx_.hdr_offsets[i + 1];

    if (off2 < off || off2 > phr_.size) {
        std::fprintf(stderr,
            "warn: bad header offsets at i=%zu: [%llu,%llu) within phr.size=%zu; skipping\n",
            i, (unsigned long long)off, (unsigned long long)off2, phr_.size);
        throw out_of_range("Phr::parse_record");
    }

    Span s{ phr_.data + off, phr_.data + off2 };
    strip_len_prefix_if_fits(s);

    TLV top;
    if (!parse_tlv_header(s, top)) return false;
    if (!(top.cls == TC_Universal && top.constructed && top.tagnum == 0x10)) return false;

    if (!top.indefinite) {
        Span body = limit_span(top);
        if (!parse_blast_def_line_fields_inplace(body, titles, ids, taxids)) {
            s.debug_print(i);
            return false;
        }
    }
    else {
        // Indefinite: work on a local span (the slice already bounds the record)
        Span body = s;
        if (!parse_blast_def_line_fields_inplace(body, titles, ids, taxids)) {
            s.debug_print(i);
            return false;
        }
    }
}

static inline uint32_t be32u(const uint8_t* p) {
    return (uint32_t)p[0] << 24 | (uint32_t)p[1] << 16 | (uint32_t)p[2] << 8 | (uint32_t)p[3];
}
static inline uint64_t be64u(const uint8_t* p) {
    return ((uint64_t)be32u(p) << 32) | (uint64_t)be32u(p + 4);
}
static inline uint64_t le64u(const uint8_t* p) {
    return (uint64_t)p[0] | (uint64_t)p[1] << 8 | (uint64_t)p[2] << 16 | (uint64_t)p[3] << 24 |
        (uint64_t)p[4] << 32 | (uint64_t)p[5] << 40 | (uint64_t)p[6] << 48 | (uint64_t)p[7] << 56;
}

static bool parse_pin_db_any(const MappedFile& pin, size_t phr_size, PinIndex& out) {
    const uint8_t* const base = pin.data;
    const uint8_t* const end = pin.data + pin.size;
    if (!base || pin.size < 16) return false;

    // Read version & dbtype at fixed positions (present in both v4/v5)
    out.version = be32u(base + 0);
    out.dbtype = be32u(base + 4);

    // Scan from after version/dbtype; v5 packs a variable "title area" we ignore.
    const uint8_t* p = base + 8;

    auto valid_monotone = [](const std::vector<uint64_t>& v) -> bool {
        if (v.empty() || v.front() != 0) return false;
        for (size_t i = 1; i < v.size(); ++i) if (v[i] < v[i - 1]) return false;
        return true;
        };

    auto try_table32 = [&](const uint8_t* t, uint32_t count, std::vector<uint64_t>& outv) -> bool {
        size_t need = (size_t)count * 4;
        if ((size_t)(end - t) < need) return false;
        if (be32u(t + 0) != 0) return false;           // first offset must be 0
        outv.resize(count);
        for (uint32_t i = 0; i < count; ++i) outv[i] = be32u(t + 4ull * i);
        return valid_monotone(outv) && outv.back() == (uint64_t)phr_size;
        };

    auto try_table64 = [&](const uint8_t* t, uint32_t count, std::vector<uint64_t>& outv) -> bool {
        size_t need = (size_t)count * 8;
        if ((size_t)(end - t) < need) return false;
        if (be64u(t + 0) != 0) return false;           // first offset must be 0
        outv.resize(count);
        for (uint32_t i = 0; i < count; ++i) outv[i] = be64u(t + 8ull * i);
        return valid_monotone(outv) && outv.back() == (uint64_t)phr_size;
        };

    // Minimal bytes needed to test a candidate header (nseq+res+maxlen) and at least one offset
    const size_t MIN_CANDIDATE = 4 + 8 + 4 + 4;

    for (; p + MIN_CANDIDATE <= end; ++p) {
        uint32_t nseq = be32u(p + 0);
        if (nseq == 0 || nseq > 500000000u) continue;  // sanity

        uint64_t residues = le64u(p + 4);
        if (residues == 0 || residues > (1ull << 48)) continue;

        uint32_t maxlen = be32u(p + 12);
        if (maxlen == 0 || maxlen > 20000000u) continue;

        // Candidate found; offsets should start immediately after this triple.
        const uint8_t* toff = p + 16;
        std::vector<uint64_t> hdr;

        // Try 32-bit, then 64-bit.
        if (try_table32(toff, nseq + 1, hdr) || try_table64(toff, nseq + 1, hdr)) {
            out.nseq = nseq;
            out.residue_count = residues;   // informational
            out.max_seq_len = maxlen;       // informational
            out.hdr_offsets = std::move(hdr);
            return true;
        }
    }

    return false; // no matching header+offset table found
}

Phr::Phr(const string &phr_path, const string& pin_path) {
    MappedFile pin;
    if (!map_file(pin_path.c_str(), pin)) die("cannot map .pin");
    if (!map_file(phr_path.c_str(), phr_)) die("cannot map .phr");
    if (!parse_pin_db_any(pin, phr_.size, idx_))
        die("failed to parse .pin (dbV4 or dbV5 expected)");
    unmap_file(pin);
}

Phr::~Phr() {
    unmap_file(phr_);
}