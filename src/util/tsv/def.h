#pragma once
#include <vector>
#include <stdexcept>

namespace Util { namespace Tsv {

enum class Type {
	STRING, INT64
};

using Schema = std::vector<Type>;
using RecordId = int64_t;

struct InvalidType : public std::runtime_error {
	InvalidType() :
		std::runtime_error("Invalid type in schema.")
	{}
};

struct SchemaMismatch : public std::runtime_error {
	SchemaMismatch() :
		std::runtime_error("Mismatching schema.")
	{}
};

}}
