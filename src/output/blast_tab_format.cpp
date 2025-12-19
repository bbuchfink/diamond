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

#include <sstream>
#include <set>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <functional>
#include <map>
#include "basic/match.h"
#include "output_format.h"
#include "util/sequence/sequence.h"
#include "data/sequence_file.h"
#include "data/taxonomy_nodes.h"
#include "util/system/system.h"

using namespace Output;
using std::array;
using std::endl;
using std::set;
using std::string;
using std::transform;
using std::back_inserter;
using std::vector;
using std::runtime_error;
using std::map;

map<FieldId, OutputField> TabularFormat::field_def = {
{ FieldId::QSeqId, { FieldId::QSeqId, "qseqid", "cseqid", "Query Seq - id", HspValues::NONE, Flags::IS_STRING } },
{ FieldId::QLen, { FieldId::QLen, "qlen", "clen", "Query sequence length", HspValues::NONE, Flags::NONE } },
{ FieldId::SSeqId, { FieldId::SSeqId, "sseqid",    "mseqid", "Subject Seq - id", HspValues::NONE, Flags::IS_STRING | Flags::SSEQID } },
{ FieldId::SAllSeqId, { FieldId::SAllSeqId, "sallseqid", "", "All subject Seq - id(s), separated by a ';'", HspValues::NONE, Flags::ALL_SEQIDS | Flags::IS_ARRAY | Flags::SSEQID } },
{ FieldId::SLen, { FieldId::SLen, "slen", "mlen", "Subject sequence length", HspValues::NONE, Flags::NONE } },
{ FieldId::QStart, { FieldId::QStart, "qstart", "cstart", "Start of alignment in query", HspValues::QUERY_START, Flags::NONE } },
{ FieldId::QEnd, { FieldId::QEnd, "qend", "cend", "End of alignment in query", HspValues::QUERY_END, Flags::NONE } },
{ FieldId::SStart, { FieldId::SStart, "sstart", "mstart", "Start of alignment in subject", HspValues::TARGET_START, Flags::NONE } },
{ FieldId::SEnd, { FieldId::SEnd, "send", "mend", "End of alignment in subject", HspValues::TARGET_END, Flags::NONE } },
{ FieldId::QSeq, { FieldId::QSeq, "qseq", "", "Aligned part of query sequence", HspValues::QUERY_COORDS, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::SSeq, { FieldId::SSeq, "sseq", "", "Aligned part of subject sequence", HspValues::TRANSCRIPT, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::EValue, { FieldId::EValue, "evalue", "evalue", "Expect value", HspValues::NONE, Flags::NONE } },
{ FieldId::BitScore, { FieldId::BitScore, "bitscore", "Bitscore", "Bit score", HspValues::NONE, Flags::NONE } },
{ FieldId::Score, { FieldId::Score, "score", "score", "Raw score", HspValues::NONE, Flags::NONE } },
{ FieldId::Length, { FieldId::Length, "length", "length", "Alignment length", HspValues::LENGTH, Flags::NONE } },
{ FieldId::PIdent, { FieldId::PIdent, "pident", "pident" ,"Percentage of identical matches", HspValues::IDENT | HspValues::LENGTH, Flags::NONE } },
{ FieldId::NIdent, { FieldId::NIdent, "nident", "nident", "Number of identical matches", HspValues::IDENT, Flags::NONE } },
{ FieldId::Mismatch, { FieldId::Mismatch, "mismatch", "mismatch", "Number of mismatches", HspValues::MISMATCHES, Flags::NONE } },
{ FieldId::Positive, { FieldId::Positive, "positive", "positive", "Number of positive - scoring matches", HspValues::TRANSCRIPT, Flags::NONE } },
{ FieldId::GapOpen, { FieldId::GapOpen, "gapopen", "gapopen", "Number of gap openings", HspValues::GAP_OPENINGS, Flags::NONE } },
{ FieldId::Gaps, { FieldId::Gaps, "gaps", "gaps", "Total number of gaps", HspValues::GAPS, Flags::NONE } },
{ FieldId::PPos, { FieldId::PPos, "ppos", "ppos", "Percentage of positive - scoring matches", HspValues::TRANSCRIPT, Flags::NONE } },
{ FieldId::QFrame, { FieldId::QFrame, "qframe", "", "Query frame", HspValues::NONE, Flags::NO_REALIGN  } },
{ FieldId::BTop, { FieldId::BTop, "btop", "", "Blast traceback operations (BTOP)", HspValues::TRANSCRIPT, Flags::IS_STRING | Flags::NO_REALIGN  } },
{ FieldId::STaxIds, { FieldId::STaxIds, "staxids", "", "Unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)", HspValues::NONE, Flags::IS_ARRAY | Flags::NO_REALIGN  } },
{ FieldId::SSciNames, { FieldId::SSciNames, "sscinames", "", "Unique Subject Scientific Name(s), separated by a ';'", HspValues::NONE, Flags::IS_ARRAY | Flags::NO_REALIGN  } },
{ FieldId::SSKingdoms, { FieldId::SSKingdoms, "sskingdoms",    "", "Unique Subject Super Kingdom(s), separated by a ';'", HspValues::NONE, Flags::IS_ARRAY | Flags::NO_REALIGN  } },
{ FieldId::STitle, { FieldId::STitle, "stitle", "", "Subject Title", HspValues::NONE, Flags::FULL_TITLES | Flags::IS_STRING  | Flags::SSEQID } },
{ FieldId::SAllTitles, { FieldId::SAllTitles, "salltitles", "", "All Subject Title(s), separated by a '<>'", HspValues::NONE, Flags::ALL_SEQIDS | Flags::FULL_TITLES | Flags::IS_ARRAY | Flags::SSEQID } },
{ FieldId::QCovHsp, { FieldId::QCovHsp, "qcovhsp", "ccovhsp", "Query coverage per HSP", HspValues::QUERY_COORDS, Flags::NONE } },
{ FieldId::QTitle, { FieldId::QTitle, "qtitle", "", "Query title", HspValues::NONE, Flags::IS_STRING } },
{ FieldId::FullSSeq, { FieldId::FullSSeq, "full_sseq", "", "Subject sequence", HspValues::NONE, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::QQual, { FieldId::QQual, "qqual", "", "Query quality values for the aligned part of the query", HspValues::QUERY_COORDS, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::QNum, { FieldId::QNum, "qnum", "", "Query ordinal id", HspValues::NONE, Flags::NONE } },
{ FieldId::SNum, { FieldId::SNum, "snum", "", "Subject ordinal id", HspValues::NONE, Flags::NONE } },
{ FieldId::SCovHsp, { FieldId::SCovHsp, "scovhsp", "mcovhsp", "Subject coverage per HSP", HspValues::TARGET_COORDS, Flags::NONE } },
{ FieldId::FullQQual, { FieldId::FullQQual, "full_qqual", "", "Query quality values", HspValues::NONE, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::FullQSeq, { FieldId::FullQSeq, "full_qseq", "", "Query sequence", HspValues::NONE, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::QSeqGapped, { FieldId::QSeqGapped, "qseq_gapped", "", "Aligned part of query sequence (with gaps)", HspValues::TRANSCRIPT, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::SSeqGapped, { FieldId::SSeqGapped, "sseq_gapped", "", "Aligned part of subject sequence (with gaps)", HspValues::TRANSCRIPT, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::QStrand, { FieldId::QStrand, "qstrand", "", "Query strand", HspValues::NONE, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::Cigar, { FieldId::Cigar, "cigar", "", "CIGAR string", HspValues::TRANSCRIPT, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::SKingdoms, { FieldId::SKingdoms, "skingdoms", "", "Unique Subject Kingdom(s), separated by a ';'", HspValues::NONE, Flags::IS_ARRAY | Flags::NO_REALIGN } },
{ FieldId::SPhylums, { FieldId::SPhylums, "sphylums", "", "Unique Subject Phylum(s), separated by a ';'", HspValues::NONE, Flags::IS_ARRAY | Flags::NO_REALIGN } },
{ FieldId::FullQSeqMate, { FieldId::FullQSeqMate, "full_qseq_mate", "", "Query sequence of the mate", HspValues::NONE, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::QSeqTranslated, { FieldId::QSeqTranslated, "qseq_translated", "", "Aligned part of query sequence (translated)", HspValues::TRANSCRIPT, Flags::IS_STRING | Flags::NO_REALIGN } },
{ FieldId::HspNum, { FieldId::HspNum, "hspnum", "", "Number of HSP within the subject", HspValues::NONE, Flags::NONE } },
{ FieldId::NormalizedBitscore, { FieldId::NormalizedBitscore, "normalized_bitscore", "", "Bitscore normalized by maximum self alignment score", HspValues::NONE, Flags::SELF_ALN_SCORES | Flags::NO_REALIGN } },
{ FieldId::NORMALIZED_NIDENT, { FieldId::NORMALIZED_NIDENT, "normalized_nident", "normalized_nident", "Number of identical matches normalized by maximum length", HspValues::IDENT | HspValues::LENGTH, Flags::NONE } },
{ FieldId::ApproxPIdent, { FieldId::ApproxPIdent, "approx_pident", "approx_pident", "Approximate percentage of identical matches", HspValues::COORDS, Flags::NONE } },
{ FieldId::CorrectedBitScore, { FieldId::CorrectedBitScore, "corrected_bitscore", "corrected_bitscore", "Bit score corrected for edge effects", HspValues::NONE, Flags::NONE } },
{ FieldId::SLineages, { FieldId::SLineages, "slineages", "", "Unique Subject Lineage(s), separated by a '<>'", HspValues::NONE, Flags::NO_REALIGN } }
#ifdef DP_STAT
{ FieldId::Reserved1, { FieldId::Reserved1, "reserved1", "reserved1", "", HspValues::NONE, Flags::NO_REALIGN } },
{ FieldId::Reserved2, { FieldId::Reserved2, "reserved2", "reserved2", "", HspValues::NONE, Flags::NO_REALIGN } }
#endif
};

template<>
struct EnumTraits<Header> {
	static const SEMap<Header> from_string;
};

const SEMap<Header> EnumTraits<Header>::from_string = { {"0", Header::NONE}, {"simple", Header::SIMPLE}, {"verbose", Header::VERBOSE} };

static std::function<void(const TabularFormat&, const HspContext&, Output::Info&)> make_invalid_match_handler(const string& key)
{
    return [key](const TabularFormat&, const HspContext&, Output::Info&) {
        throw runtime_error(string("Invalid output field: ") + key);
        };
}

static std::function<void(const TabularFormat&, Output::Info&)> make_invalid_intro_handler(const string& key)
{
    return [key](const TabularFormat&, Output::Info&) {
        throw runtime_error(string("Invalid output field: ") + key);
        };
}

static void print_staxids(TextBuffer& out, int64_t subject_global_id, const SequenceFile& db, bool json = false)
{
    out.print(db.taxids(subject_global_id), json ? ',' : ';');
}

static vector<TaxId> lineage(TaxId taxid, SequenceFile& db) {
    vector<TaxId> l;
    TaxId i = taxid;
    int n = 0;
    while (true) {
        if (i <= 0)
            return {};
        if (i == 1)
            break;
        l.push_back(i);
        i = db.get_parent(i);
        if (l.size() >= MAX_LINEAGE)
            throw runtime_error("Lineage too long for taxid " + std::to_string(taxid));
    }
    return l;
}

static void print_lineage(int64_t target_oid, SequenceFile& db, TextBuffer& out, bool json = false) {
    const vector<TaxId> taxids = db.taxids(target_oid);
    if (taxids.empty()) {
        out << (json ? " []" : "N/A");
        return;
    }
    set<vector<TaxId>> lineages;
    for (TaxId i : taxids) {
        const vector<TaxId> l = lineage(i, db);
        if (!l.empty())
            lineages.insert(l);
    }
    if (lineages.empty()) {
        out << (json ? " []" : "N/A");
        return;
    }
    if (json) {
        out << " [" DEFAULT_LINE_DELIMITER;
    }
    for (auto i = lineages.begin(); i != lineages.end(); ++i) {
        if (i != lineages.begin())
            out << (json ? "," DEFAULT_LINE_DELIMITER : "<>");
        if (json) {
            out << "\t\t[";
        }
        const vector<TaxId>& lin = *i;
        for (auto j = lin.rbegin(); j != lin.rend(); ++j) {
            if (j != lin.rbegin())
                out << (json ? ", " : "; ");
            if (json)
                out << '\"';
            out << db.taxon_scientific_name(*j);
            if (json)
                out << '\"';
        }
        if (json) {
            out << ']';
        }
    }
    if (json) {
        out << DEFAULT_LINE_DELIMITER "\t]";
    }
}

map<FieldId, FieldCallbacks> TabularFormat::field_callbacks = [] {
    map<FieldId, FieldCallbacks> callbacks;
    for (const auto& [id, field] : TabularFormat::field_def) {
        callbacks.emplace(id, FieldCallbacks{ make_invalid_match_handler(field.key), make_invalid_intro_handler(field.key) });
    }

    callbacks[FieldId::QSeqId].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out.write_until(r.query_title.c_str(), Util::Seq::id_delimiters);
        };
    callbacks[FieldId::QSeqId].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out.write_until(info.query.title, Util::Seq::id_delimiters);
        };

    callbacks[FieldId::QLen].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.query_len;
        };
    callbacks[FieldId::QLen].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << info.query.len;
        };

    callbacks[FieldId::SSeqId].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        OutputFormat::print_title(info.out, r.target_title.c_str(), false, false, "");
        };
    callbacks[FieldId::SSeqId].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::SAllSeqId].match = [](const TabularFormat& format, const HspContext& r, Output::Info& info) {
        OutputFormat::print_title(info.out, r.target_title.c_str(), false, true, format.is_json ? "," : ";", 0, format.is_json);
        };
    callbacks[FieldId::SAllSeqId].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::SLen].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.subject_len;
        };
    callbacks[FieldId::SLen].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::QStart].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.oriented_query_range().begin_ + 1;
        };
    callbacks[FieldId::QStart].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::QEnd].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.oriented_query_range().end_ + 1;
        };
    callbacks[FieldId::QEnd].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::SStart].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.subject_source_range().begin_ + 1;
        };
    callbacks[FieldId::SStart].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::SEnd].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.subject_source_range().end_;
        };
    callbacks[FieldId::SEnd].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::QSeq].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        r.query.source().print(info.out, r.query_source_range().begin_, r.query_source_range().end_, input_value_traits);
        };
    callbacks[FieldId::QSeq].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::SSeq].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        vector<Letter> seq;
        seq.reserve(r.subject_range().length());
        for (HspContext::Iterator j = r.begin(); j.good(); ++j)
            if (!(j.op() == op_insertion))
                seq.push_back(j.subject());
        info.out << Sequence(seq);
        };
    callbacks[FieldId::SSeq].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::EValue].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out.print_e(r.evalue());
        };
    callbacks[FieldId::EValue].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::BitScore].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.bit_score();
        };
    callbacks[FieldId::BitScore].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::Score].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.score();
        };
    callbacks[FieldId::Score].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::Length].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.length();
        };
    callbacks[FieldId::Length].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::PIdent].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.id_percent();
        };
    callbacks[FieldId::PIdent].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::NIdent].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.identities();
        };
    callbacks[FieldId::NIdent].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::Mismatch].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.mismatches();
        };
    callbacks[FieldId::Mismatch].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::Positive].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.positives();
        };
    callbacks[FieldId::Positive].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::GapOpen].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.gap_openings();
        };
    callbacks[FieldId::GapOpen].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::Gaps].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.gaps();
        };
    callbacks[FieldId::Gaps].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::PPos].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << (double)r.positives() * 100.0 / r.length();
        };
    callbacks[FieldId::PPos].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::QFrame].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.blast_query_frame();
        };
    callbacks[FieldId::QFrame].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '0';
        };

    callbacks[FieldId::BTop].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        unsigned n_matches = 0;
        for (HspContext::Iterator i = r.begin(); i.good(); ++i) {
            switch (i.op()) {
            case op_match:
                ++n_matches;
                break;
            case op_substitution:
            case op_frameshift_forward:
            case op_frameshift_reverse:
                if (n_matches > 0) {
                    info.out << n_matches;
                    n_matches = 0;
                }
                info.out << i.query_char() << i.subject_char();
                break;
            case op_insertion:
                if (n_matches > 0) {
                    info.out << n_matches;
                    n_matches = 0;
                }
                info.out << i.query_char() << '-';
                break;
            case op_deletion:
                if (n_matches > 0) {
                    info.out << n_matches;
                    n_matches = 0;
                }
                info.out << '-' << i.subject_char();
                break;
            }
        }
        if (n_matches > 0)
            info.out << n_matches;
        };
    callbacks[FieldId::BTop].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::STaxIds].match = [](const TabularFormat& format, const HspContext& r, Output::Info& info) {
        print_staxids(info.out, r.subject_oid, *info.db, format.is_json);
        };
    callbacks[FieldId::STaxIds].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::SSciNames].match = [](const TabularFormat& format, const HspContext& r, Output::Info& info) {
        const vector<TaxId> tax_id = info.db->taxids(r.subject_oid);
        print_taxon_names(tax_id.begin(), tax_id.end(), *info.db, info.out, format.is_json);
        };
    callbacks[FieldId::SSciNames].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::SSKingdoms].match = [](const TabularFormat& format, const HspContext& r, Output::Info& info) {
        const set<TaxId> tax_id = info.db->rank_taxid(info.db->taxids(r.subject_oid), Rank::superkingdom);
        print_taxon_names(tax_id.begin(), tax_id.end(), *info.db, info.out, format.is_json);
        };
    callbacks[FieldId::SSKingdoms].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::STitle].match = [](const TabularFormat& format, const HspContext& r, Output::Info& info) {
        OutputFormat::print_title(info.out, r.target_title.c_str(), true, false, format.is_json ? "," : "<>", 0, false);
        };
    callbacks[FieldId::STitle].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::SAllTitles].match = [](const TabularFormat& format, const HspContext& r, Output::Info& info) {
        OutputFormat::print_title(info.out, r.target_title.c_str(), true, true, format.is_json ? "," : "<>", 0, format.is_json);
        };
    callbacks[FieldId::SAllTitles].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::QCovHsp].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.qcovhsp();
        };
    callbacks[FieldId::QCovHsp].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::QTitle].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.query_title;
        };
    callbacks[FieldId::QTitle].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << info.query.title;
        };

    callbacks[FieldId::FullSSeq].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.subject_seq;
        };
    callbacks[FieldId::FullSSeq].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::QQual].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        if (strlen(info.query.qual) == 0) {
            info.out << '*';
            return;
        }
        info.out << string(info.query.qual + r.query_source_range().begin_, info.query.qual + r.query_source_range().end_).c_str();
        };
    callbacks[FieldId::QQual].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::QNum].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.query_oid;
        };
    callbacks[FieldId::QNum].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::SNum].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.subject_oid;
        };
    callbacks[FieldId::SNum].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::SCovHsp].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.scovhsp();
        };
    callbacks[FieldId::SCovHsp].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::FullQQual].match = [](const TabularFormat&, const HspContext&, Output::Info& info) {
        info.out << (strlen(info.query.qual) ? info.query.qual : "*");
        };
    callbacks[FieldId::FullQQual].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << (strlen(info.query.qual) ? info.query.qual : "*");
        };

    callbacks[FieldId::FullQSeq].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        r.query.source().print(info.out, input_value_traits);
        };
    callbacks[FieldId::FullQSeq].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.query.source_seq.print(info.out, input_value_traits);
        };

    callbacks[FieldId::QSeqGapped].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        for (HspContext::Iterator i = r.begin(); i.good(); ++i)
            info.out << i.query_char();
        };
    callbacks[FieldId::QSeqGapped].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::SSeqGapped].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        for (HspContext::Iterator i = r.begin(); i.good(); ++i)
            info.out << i.subject_char();
        };
    callbacks[FieldId::SSeqGapped].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::QStrand].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        if (align_mode.query_translated)
            info.out << ((r.blast_query_frame() > 0) ? '+' : '-');
        else
            info.out << '+';
        };
    callbacks[FieldId::QStrand].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::Cigar].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        print_cigar(r, info.out);
        };
    callbacks[FieldId::Cigar].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::SKingdoms].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        const set<TaxId> tax_id = info.db->rank_taxid(info.db->taxids(r.subject_oid), Rank::kingdom);
        print_taxon_names(tax_id.begin(), tax_id.end(), *info.db, info.out);
        };
    callbacks[FieldId::SKingdoms].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::SPhylums].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        const set<TaxId> tax_id = info.db->rank_taxid(info.db->taxids(r.subject_oid), Rank::phylum);
        print_taxon_names(tax_id.begin(), tax_id.end(), *info.db, info.out);
        };
    callbacks[FieldId::SPhylums].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::FullQSeqMate].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        if (config.query_file.size() == 2) {
            info.query.mate_seq.print(info.out, input_value_traits);
            return;
        }
        info.out << '*';
        };

    callbacks[FieldId::QSeqTranslated].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        if (config.frame_shift) {
            vector<Letter> seq;
            seq.reserve(r.query_range().length());
            for (HspContext::Iterator j = r.begin(); j.good(); ++j)
                if (j.op() != op_deletion && j.op() != op_frameshift_forward && j.op() != op_frameshift_reverse)
                    seq.push_back(j.query());
            info.out << Sequence(seq);
            return;
        }
        r.query.index(r.frame()).print(info.out, r.query_range().begin_, r.query_range().end_, amino_acid_traits);
        };
    callbacks[FieldId::QSeqTranslated].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::ApproxPIdent].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.approx_id();
        };
    callbacks[FieldId::ApproxPIdent].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << "-1";
        };

    callbacks[FieldId::CorrectedBitScore].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.corrected_bit_score();
        };

    callbacks[FieldId::SLineages].match = [](const TabularFormat& format, const HspContext& r, Output::Info& info) {
        print_lineage(r.subject_oid, *info.db, info.out, format.is_json);
        };
    callbacks[FieldId::SLineages].query_intro = [](const TabularFormat&, Output::Info& info) {
        info.out << '*';
        };

    callbacks[FieldId::HspNum].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.hsp_num;
        };

    callbacks[FieldId::NormalizedBitscore].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out.print_d(r.bit_score() / std::max(r.query_self_aln_score, r.target_self_aln_score));
        };

    callbacks[FieldId::NORMALIZED_NIDENT].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out.print_d((double)r.identities() / std::max(r.query.index(r.frame()).length(), r.subject_len));
        };

#ifdef DP_STAT
    callbacks[FieldId::Reserved1].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.reserved1();
        };
    callbacks[FieldId::Reserved2].match = [](const TabularFormat&, const HspContext& r, Output::Info& info) {
        info.out << r.reserved2();
        };
#endif

    return callbacks;
    }();

Header TabularFormat::header_format(unsigned workflow) {
	const bool cluster = workflow == Config::cluster || workflow == Config::DEEPCLUST;
	if (workflow != Config::blastp && !cluster)
		throw runtime_error("header_format");
	if (!config.output_header.present())
		return Header::NONE;
	if (config.output_header.empty())
		return cluster ? Header::SIMPLE : Header::VERBOSE;
	if (config.output_header.size() > 1)
		throw runtime_error("Invalid header format: " + join(" ", config.output_header.cbegin(), config.output_header.cend()));
	const Header h = from_string<Header>(config.output_header.front());
	if (h == Header::VERBOSE && cluster)
		throw runtime_error("Verbose header format is not supported for cluster workflow.");
	return h;
}

TabularFormat::TabularFormat(bool json) :
    OutputFormat((json ? OutputFormat::json : OutputFormat::blast_tab), HspValues::NONE, Flags::NONE, json ? ',' : '\0'),
    is_json(json)
{
    static const FieldId stdf[] = { FieldId::QSeqId, FieldId::SSeqId, FieldId::PIdent, FieldId::Length, FieldId::Mismatch, FieldId::GapOpen, FieldId::QStart, FieldId::QEnd, FieldId::SStart, FieldId::SEnd, FieldId::EValue, FieldId::BitScore };
    const vector<string>& f = config.output_format;
    if (f.size() <= 1) {
        fields = vector<FieldId>(stdf, stdf + 12);
        if (config.frame_shift == 0)
            hsp_values = HspValues::QUERY_COORDS | HspValues::TARGET_COORDS | HspValues::LENGTH | HspValues::IDENT | HspValues::MISMATCHES | HspValues::GAP_OPENINGS;
        else
            hsp_values = HspValues::TRANSCRIPT;
        flags |= Output::Flags::SSEQID;
        return;
    }
    for (vector<string>::const_iterator i = f.begin() + 1; i != f.end(); ++i) {
        auto it = std::find_if(field_def.begin(), field_def.end(), [i](const auto& f) {return f.second.key == *i; });
        if (it == field_def.end())
            throw std::runtime_error(string("Invalid output field: ") + *i);
        const FieldId id = it->first;
        if (id == FieldId::STaxIds)
            needs_taxon_id_lists = true;
        if (id == FieldId::SSciNames || id == FieldId::SSKingdoms || id == FieldId::SKingdoms || id == FieldId::SPhylums || id == FieldId::SLineages || id >= FieldId::COUNT) {
            needs_taxon_scientific_names = true;
            needs_taxon_id_lists = true;
        }
        if (id == FieldId::SSKingdoms || id == FieldId::SKingdoms || id == FieldId::SPhylums || id >= FieldId::COUNT) {
            needs_taxon_nodes = true;
            needs_taxon_ranks = true;
        }
        if (id == FieldId::SLineages) {
            needs_taxon_nodes = true;
        }
        fields.push_back(id);
        if (id == FieldId::FullSSeq || id == FieldId::ApproxPIdent)
            flags |= Flags::TARGET_SEQS;
        if (id == FieldId::QQual || id == FieldId::FullQQual)
            config.store_query_quality = true;
        if (id == FieldId::FullQSeqMate)
            needs_paired_end_info = true;
        if (id == FieldId::ApproxPIdent && score_matrix.name() != "BLOSUM62")
            throw std::runtime_error("Approximate identity is only supported for the BLOSUM62 scoring matrix.");
        if ((id == FieldId::FullQSeqMate || id == FieldId::QSeqTranslated) && !align_mode.query_translated)
            throw std::runtime_error("Output field only supported for translated search.");
        hsp_values |= field_def.at(id).hsp_values;
        flags |= field_def.at(id).flags;
    }
    //if (config.traceback_mode == TracebackMode::NONE && config.max_hsps == 1 && !needs_transcript && !needs_stats && !config.query_range_culling && config.min_id == 0.0 && config.query_cover == 0.0 && config.subject_cover == 0.0)
        //config.traceback_mode = TracebackMode::SCORE_ONLY;
}

void TabularFormat::print_match(const HspContext& r, Output::Info& info)
{
    TextBuffer& out = info.out;
    const char* prepos = "\t";
    if (is_json) {
        if (r.hit_num != 0)
            out << ",";
        out << "\n\t{\n";
    }
    for (auto i = fields.cbegin(); i != fields.cend(); ++i) {
        const OutputField& field = field_def.at(*i);
        if (is_json) {
            out << prepos << "\"" << field.key << "\":";
            if (flag_any(field.flags, Flags::IS_STRING))
                out << "\"";
            if (flag_any(field.flags, Flags::IS_ARRAY))
                out << "[";
        }
        const FieldCallbacks& callbacks = field_callbacks.at(*i);
        callbacks.match(*this, r, info);
        if (is_json) {
            if (flag_any(field.flags, Flags::IS_STRING))
                out << "\"";
            if (flag_any(field.flags, Flags::IS_ARRAY))
                out << "]";
            if (i < fields.end() - 1)
                out << ",\n";
            else
                out << "\n";
        }
        else if (i < fields.end() - 1)
            out << '\t';
    }
    out << (is_json ? "\t}" : "\n");
}

void TabularFormat::print_query_intro(Output::Info& info) const
{
    TextBuffer& out = info.out;
    if (info.unaligned && config.report_unaligned == 1) {
        for (auto i = fields.cbegin(); i != fields.cend(); ++i) {
            const FieldCallbacks& callbacks = field_callbacks.at(*i);
            callbacks.query_intro(*this, info);
            if (i < fields.end() - 1)
                out << '\t';
        }
        out << '\n';
    }
}

void TabularFormat::output_header(Consumer& f, bool cluster) const {
    vector<string> headers;
    transform(fields.begin(), fields.end(), back_inserter(headers), [cluster](FieldId i) {
        const OutputField& field = TabularFormat::field_def.at(i);
        const string& key = cluster ? field.clust_key : field.key;
        if (cluster && key.empty())
            throw runtime_error("Output field not supported for clustering: " + field.key);
        return key;
        });
    const string s = join("\t", headers.begin(), headers.end()) + '\n';
    f.consume(s.data(), s.length());
    return;
}

void TabularFormat::print_header(Consumer& f, int mode, const char* matrix, int gap_open, int gap_extend, double evalue, const char* first_query_name, unsigned first_query_len) const {
    const Header h = header_format(Config::blastp);
    if (h == Header::VERBOSE) {
        std::ostringstream ss;
        ss << "# DIAMOND v" << Const::version_string << ". http://github.com/bbuchfink/diamond" << endl;
        ss << "# Invocation: " << config.invocation << endl;
        ss << "# Fields: ";
		for (auto i = fields.begin(); i != fields.end(); ++i) {
            if(i != fields.begin())
				ss << ", ";
            ss << field_def.at(*i).description;
        }
        const string s(ss.str());
        f.consume(s.data(), s.length());
    }
    else if (h == Header::SIMPLE) {
        output_header(f, false);
    }
    if (is_json) {
        const char* hdr = "[";
        f.consume(hdr, strlen(hdr));
    }
}

void TabularFormat::print_footer(Consumer &f) const {
    if(is_json) {
        const char* s = "\n]";
		f.consume(s, strlen(s));
    }
}