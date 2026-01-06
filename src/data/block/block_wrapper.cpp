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

#include "block_wrapper.h"

using std::vector;
using std::string;

BlockWrapper::BlockWrapper(const Block& block, Flags flags, const ValueTraits& value_traits) :
	SequenceFile(SequenceFile::Type::BLOCK, flags, FormatFlags::LENGTH_LOOKUP | FormatFlags::TITLES_LAZY | FormatFlags::SEEKABLE, value_traits),
	block_(block),
	oid_(0)
{
}

int64_t BlockWrapper::file_count() const {
	return 1;
}

bool BlockWrapper::files_synced() {
	return true;
}

SequenceFile::SeqInfo BlockWrapper::read_seqinfo() {
	if (oid_ >= block_.seqs().size()) {
		++oid_;
		return SeqInfo(0, 0);
	}
	const Loc l = block_.seqs().length(oid_);
	if (l == 0)
		throw std::runtime_error("Database with sequence length 0 is not supported");
	return SeqInfo(oid_++, l);
}

void BlockWrapper::putback_seqinfo() {
	--oid_;
}

void BlockWrapper::close() {
}

void BlockWrapper::set_seqinfo_ptr(OId i) {
	oid_ = i;
}

OId BlockWrapper::tell_seq() const {
	return oid_;
}

bool BlockWrapper::eof() const {
	return oid_ >= block_.seqs().size();
}

void BlockWrapper::init_seq_access() {
	set_seqinfo_ptr(0);
}

bool BlockWrapper::read_seq(vector<Letter>& seq, string& id, std::vector<char>* quals)
{
	throw OperationNotSupported();
}

void BlockWrapper::create_partition_balanced(int64_t max_letters) {
	throw OperationNotSupported();
}

void BlockWrapper::save_partition(const string& partition_file_name, const string& annotation) {
	throw OperationNotSupported();
}

int BlockWrapper::get_n_partition_chunks() {
	throw OperationNotSupported();
}

void BlockWrapper::init_seqinfo_access() {
}

void BlockWrapper::seek_chunk(const Chunk& chunk) {
	throw OperationNotSupported();
}

std::string BlockWrapper::seqid(OId oid, bool all, bool full_titles) {
	return block_.ids()[oid];
}

size_t BlockWrapper::id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) {
	return block_.ids().length(seq_info.pos);
}

void BlockWrapper::seek_offset(size_t p) {
}

void BlockWrapper::read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) {
	*(dst - 1) = Sequence::DELIMITER;
	*(dst + len) = Sequence::DELIMITER;
	std::copy(block_.seqs().ptr(pos), block_.seqs().end(pos), dst);
	++pos;
}

void BlockWrapper::read_id_data(const int64_t oid, char* dst, size_t len, bool all, bool full_titles) {
	std::copy(block_.ids().ptr(oid), block_.ids().end(oid), dst);
	dst[len] = '\0';
}

void BlockWrapper::skip_id_data() {
}

uint64_t BlockWrapper::sequence_count() const {
	return block_.seqs().size();
}

uint64_t BlockWrapper::letters() const {
	return block_.seqs().letters();
}

int BlockWrapper::db_version() const {
	throw OperationNotSupported();
}

int BlockWrapper::program_build_version() const {
	throw OperationNotSupported();
}

int BlockWrapper::build_version() {
	throw OperationNotSupported();
}

BlockWrapper::~BlockWrapper()
{
}

DbFilter* BlockWrapper::filter_by_accession(const std::string& file_name)
{
	throw OperationNotSupported();
}

std::string BlockWrapper::file_name()
{
	return {};
}

std::vector<TaxId> BlockWrapper::taxids(size_t oid) const
{
	throw OperationNotSupported();
}

void BlockWrapper::seq_data(size_t oid, std::vector<Letter>& dst)
{
	throw OperationNotSupported();
}

Loc BlockWrapper::seq_length(size_t oid)
{
	throw OperationNotSupported();
}

void BlockWrapper::end_random_access(bool dictionary)
{
	if (!dictionary)
		return;
	free_dictionary();
}
