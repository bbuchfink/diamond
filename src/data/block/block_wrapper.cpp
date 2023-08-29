/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "block_wrapper.h"

using std::vector;
using std::string;

BlockWrapper::BlockWrapper(const Block& block, Metadata metadata, Flags flags, const ValueTraits& value_traits) :
	SequenceFile(SequenceFile::Type::BLOCK, Alphabet::STD, flags, FormatFlags::LENGTH_LOOKUP | FormatFlags::TITLES_LAZY | FormatFlags::SEEKABLE, value_traits),
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

void BlockWrapper::create_partition_balanced(size_t max_letters) {
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

std::string BlockWrapper::seqid(OId oid) const {
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

void BlockWrapper::read_id_data(const int64_t oid, char* dst, size_t len) {
	std::copy(block_.ids().ptr(oid), block_.ids().end(oid), dst);
	dst[len] = '\0';
}

void BlockWrapper::skip_id_data() {
}

int64_t BlockWrapper::sequence_count() const {
	return block_.seqs().size();
}

size_t BlockWrapper::letters() const {
	return block_.seqs().letters();
}

int BlockWrapper::db_version() const {
	throw OperationNotSupported();
}

int BlockWrapper::program_build_version() const {
	throw OperationNotSupported();
}

SequenceFile::Metadata BlockWrapper::metadata() const {
	return Metadata();
}

int BlockWrapper::build_version() {
	throw OperationNotSupported();
}

BlockWrapper::~BlockWrapper()
{
}

void BlockWrapper::close_weakly()
{
}

void BlockWrapper::reopen()
{
}

BitVector* BlockWrapper::filter_by_accession(const std::string& file_name)
{
	throw OperationNotSupported();
}

const BitVector* BlockWrapper::builtin_filter()
{
	return nullptr;
}

std::string BlockWrapper::file_name()
{
	return {};
}

int64_t BlockWrapper::sparse_sequence_count() const
{
	throw OperationNotSupported();
}

std::vector<TaxId> BlockWrapper::taxids(size_t oid) const
{
	throw OperationNotSupported();
}

void BlockWrapper::seq_data(size_t oid, std::vector<Letter>& dst) const
{
	throw OperationNotSupported();
}

size_t BlockWrapper::seq_length(size_t oid) const
{
	throw OperationNotSupported();
}

void BlockWrapper::init_random_access(const size_t query_block, const size_t ref_blocks, bool dictionary)
{
	if (dictionary)
		load_dictionary(query_block, ref_blocks);
}

void BlockWrapper::end_random_access(bool dictionary)
{
	if (!dictionary)
		return;
	free_dictionary();
}
