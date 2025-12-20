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

#include <assert.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <set>
#include "util/system/system.h"
#include "pal.h"
#include "volume.h"
#include "util/string/string.h"

using std::string;
using std::runtime_error;
using std::set;
using std::tie;
using std::vector;
using std::advance;

static string trim(const string& text)
{
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == string::npos) {
        return string();
    }
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

static vector<string> split_whitespace(const std::string& text)
{
    std::istringstream in(text);
    std::vector<string> parts;
    string token;
    while (in >> token) {
        parts.push_back(token);
    }
    return parts;
}

[[nodiscard]] std::expected<vector<string>::iterator, Diamond::Error> Pal::recurse(const string& path, vector<string>::iterator volume_it) {
	auto pal = Pal::open(path);
    if(!pal)
		return std::unexpected(pal.error());
	volume_it = volumes.insert(volume_it, pal->volumes.begin(), pal->volumes.end());
	advance(volume_it, pal->volumes.size());
    
    for(const auto& kv : pal->metadata) {
        if (metadata.find(kv.first) != metadata.end()) {
            if(kv.first == "TITLE" || kv.first == "NSEQ" || kv.first == "LENGTH")
				continue;
            return std::unexpected(Diamond::Error{ Diamond::ErrorCode::PAL_PARSE_ERROR,
                "Duplicate key '" + kv.first + "' in nested PAL file: " + path });
        }
        else
            metadata[kv.first] = kv.second;
	}
    const OId base = oid_index.back();
    for(vector<OId>::const_iterator it = pal->oid_index.begin() + 1; it != pal->oid_index.end(); ++it)
		oid_index.push_back(*it + base);
    sequence_count += pal->sequence_count;
    letters += pal->letters;
	version = pal->version;
    return volume_it;
}

[[nodiscard]] std::expected<Pal, Diamond::Error> Pal::open(const string& path) {
	const set<string> supported_keys = { "TITLE", "MEMB_BIT", "SEQIDLIST", "NSEQ", "LENGTH", "TAXIDLIST" };
    Pal pal;
    string db_dir, file;
    tie(db_dir, file) = absolute_path(path);
    if (!exists(path + ".pal") && !ends_with(path, ".pal")) {
		pal.volumes.push_back(db_dir + PATH_SEPARATOR + file);
    }
    else {
        const string pal_path = ends_with(path, ".pal") ? path : path + ".pal";
        std::ifstream in(pal_path);
        if (!in)
            return std::unexpected(Diamond::Error{ Diamond::ErrorCode::PAL_OPEN_ERROR, "Unable to open PAL file: " + pal_path });
        string line;
        size_t line_number = 0;
        while (std::getline(in, line)) {
            ++line_number;
            const auto comment = line.find('#');
            if (comment != string::npos) {
                line = line.substr(0, comment);
            }

            line = trim(line);
            if (line.empty()) {
                continue;
            }

            const auto key_end = line.find_first_of(" \t");
            if (key_end == string::npos) {
                return std::unexpected(Diamond::Error{ Diamond::ErrorCode::PAL_PARSE_ERROR,
                    "Line " + std::to_string(line_number) + " is missing a value: " + line });
            }

            const string key = line.substr(0, key_end);
            const string value = trim(line.substr(key_end + 1));
            if (value.empty()) {
                throw runtime_error("Line " + std::to_string(line_number) + " has an empty value: " + line);
            }

            if (key == "DBLIST") {
                const auto volumes = split_whitespace(value);
                if (volumes.empty()) {
                    return std::unexpected(Diamond::Error{ Diamond::ErrorCode::PAL_PARSE_ERROR,
                    "DBLIST on line " + std::to_string(line_number) + " does not list any volumes" });
                }
                pal.volumes.insert(pal.volumes.end(), volumes.begin(), volumes.end());
                for (string& s : pal.volumes)
                    if (!is_absolute_path(s) && !s.empty() && s[0] != '"')
                        s = db_dir + PATH_SEPARATOR + s;
                continue;
            }

			if (supported_keys.find(key) == supported_keys.end())
                return std::unexpected(Diamond::Error{ Diamond::ErrorCode::BLASTDB_UNSUPPORTED_FEATURE,
                    "Unsupported PAL key '" + key + "' on line " + std::to_string(line_number) });

            const auto duplicate = pal.metadata.find(key);
            if (duplicate != pal.metadata.end()) {
                return std::unexpected(Diamond::Error{ Diamond::ErrorCode::PAL_PARSE_ERROR,
                    "Duplicate key '" + key + "' on line " + std::to_string(line_number) });
            }

            pal.metadata[key] = value;
        }
    }
    pal.sequence_count = 0;
    pal.letters = 0;
    pal.oid_index.push_back(0);
    
	for (vector<string>::iterator it = pal.volumes.begin(); it != pal.volumes.end();) {
        if (it->length() >= 2 && (*it)[0] == '"' && (*it)[it->size() - 1] == '"') {
            const string nested = it->substr(1, it->length() - 2);
            it = pal.volumes.erase(it);
            auto e = pal.recurse(is_absolute_path(nested) ? nested : db_dir + PATH_SEPARATOR + nested, it);
            if(!e)
				return std::unexpected(e.error());
            it = e.value();
        }
        else {
            Volume vol = Volume::Open(*it, 0, 0, 0);
            pal.sequence_count += vol.index().num_oids;
            pal.oid_index.push_back(vol.index().num_oids + pal.oid_index.back());
            pal.letters += vol.index().total_length;
            pal.version = vol.index().version;
            ++it;
        }
    }
    if (pal.metadata.find("SEQIDLIST") != pal.metadata.end()) {
        if(ends_with(pal.metadata["SEQIDLIST"], ".bsl"))
            return std::unexpected(Diamond::Error{ Diamond::ErrorCode::BLASTDB_UNSUPPORTED_FEATURE,
				"Binary SEQIDLIST files (.bsl) are not supported, use text file instead: " + pal.metadata["SEQIDLIST"] });
        if (!is_absolute_path(pal.metadata["SEQIDLIST"])) {
            pal.metadata["SEQIDLIST"] = db_dir + PATH_SEPARATOR + pal.metadata["SEQIDLIST"];
        }
    }
    if (pal.metadata.find("TAXIDLIST") != pal.metadata.end() && !is_absolute_path(pal.metadata["TAXIDLIST"])) {
        pal.metadata["TAXIDLIST"] = db_dir + PATH_SEPARATOR + pal.metadata["TAXIDLIST"];
    }
    assert(pal.sequence_count > 0);
    return pal;
}

int Pal::volume(OId oid) const noexcept {
	assert(oid < sequence_count);
	auto it = std::upper_bound(oid_index.begin(), oid_index.end(), oid);
	return it - oid_index.begin() - 1;
}