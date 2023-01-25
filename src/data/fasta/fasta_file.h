#include <vector>
#include <list>
#include "../sequence_file.h"

struct FastaFile : public SequenceFile
{

	struct WriteAccess {};

	FastaFile(const std::vector<std::string> &file_name, Metadata metadata = Metadata(), Flags flags = Flags::NONE, const ValueTraits& value_traits = amino_acid_traits);
	FastaFile(const std::string& file_name, bool overwrite, const WriteAccess&);

	virtual int64_t file_count() const override;
	virtual void create_partition_balanced(size_t max_letters) override;
	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") override;
	virtual int get_n_partition_chunks() override;
	virtual void close() override;
	virtual void set_seqinfo_ptr(OId i) override;
	virtual OId tell_seq() const override;
	virtual bool eof() const override;
	virtual bool files_synced() override;
	virtual void init_seq_access() override;
	virtual void init_seqinfo_access() override;
	virtual void seek_chunk(const Chunk& chunk) override;
	virtual SeqInfo read_seqinfo() override;
	virtual void putback_seqinfo() override;
	virtual size_t id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) override;
	virtual void seek_offset(size_t p) override;
	virtual void read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) override;
	virtual void read_id_data(const int64_t oid, char* dst, size_t len) override;
	virtual void skip_id_data() override;
	virtual int64_t sequence_count() const override;
	virtual bool read_seq(std::vector<Letter>& seq, std::string& id, std::vector<char>* quals = nullptr) override;
	virtual size_t letters() const override;
	virtual int db_version() const override;
	virtual int program_build_version() const override;
	virtual Metadata metadata() const override;
	virtual int build_version() override;
	virtual ~FastaFile();
	virtual void close_weakly() override;
	virtual void reopen() override;
	virtual BitVector* filter_by_accession(const std::string& file_name) override;
	virtual const BitVector* builtin_filter() override;
	virtual std::string file_name() override;
	virtual int64_t sparse_sequence_count() const override;
	virtual std::vector<TaxId> taxids(size_t oid) const override;
	virtual void seq_data(size_t oid, std::vector<Letter>& dst) const override;
	virtual size_t seq_length(size_t oid) const override;
	virtual void init_random_access(const size_t query_block, const size_t ref_blocks, bool dictionary = true) override;
	virtual void end_random_access(bool dictionary = true) override;
	virtual void init_write() override;
	virtual void write_seq(const Sequence& seq, const std::string& id) override;

	int64_t line_count() const;

	static void prep_db(const std::string& path);

private:

	std::pair<int64_t, int64_t> init_read();

	std::list<TextInputFile> file_;
	std::list<TextInputFile>::iterator file_ptr_;
	std::unique_ptr<OutputFile> out_file_;
	std::unique_ptr<const Sequence_file_format> format_;
	OId oid_;
	int64_t seqs_, letters_;
	
};
