/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

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

vector<string>::iterator Pal::recurse(const string& path, vector<string>::iterator volume_it) {
	Pal pal(path);
	const auto index = std::distance(volumes.begin(), volume_it);
    volumes.insert(volume_it, pal.volumes.begin(), pal.volumes.end());
	    
    for(const auto& kv : pal.metadata) {
        if (metadata.find(kv.first) != metadata.end()) {
            if(kv.first == "TITLE" || kv.first == "NSEQ" || kv.first == "LENGTH")
				continue;
            throw runtime_error("Duplicate key '" + kv.first + "' in nested PAL file: " + path);
        }
        else
            metadata[kv.first] = kv.second;
	}
    const OId base = oid_index.back();
    for(vector<OId>::const_iterator it = pal.oid_index.begin() + 1; it != pal.oid_index.end(); ++it)
		oid_index.push_back(*it + base);
    sequence_count += pal.sequence_count;
    letters += pal.letters;
	version = pal.version;
	return volumes.begin() + index + pal.volumes.size();
}

Pal::Pal(const string& path) {
	const set<string> supported_keys = { "TITLE", "MEMB_BIT", "SEQIDLIST", "NSEQ", "LENGTH", "TAXIDLIST" };
    string db_dir, file;
    tie(db_dir, file) = absolute_path(path);
    if (!exists(path + ".pal") && !ends_with(path, ".pal")) {
		volumes.push_back(db_dir + PATH_SEPARATOR + file);
    }
    else {
        const string pal_path = ends_with(path, ".pal") ? path : path + ".pal";
        std::ifstream in(pal_path);
        if (!in)
            throw runtime_error("Unable to open PAL file: " + pal_path);
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
                throw runtime_error("Error parsing PAL file: line " + std::to_string(line_number) + " is missing a value: " + line);
            }

            const string key = line.substr(0, key_end);
            const string value = trim(line.substr(key_end + 1));
            if (value.empty())
                throw runtime_error("Error parsing PAL file: line " + std::to_string(line_number) + " has an empty value: " + line);

            if (key == "DBLIST") {
                const auto vls = split_whitespace(value);
                if (vls.empty()) {
                    throw runtime_error("Error parsing PAL file: DBLIST on line " + std::to_string(line_number) + " does not list any volumes");
                }
                volumes.insert(volumes.end(), vls.begin(), vls.end());
                for (string& s : volumes)
                    if (!is_absolute_path(s) && !s.empty() && s[0] != '"')
                        s = db_dir + PATH_SEPARATOR + s;
                continue;
            }

            if (supported_keys.find(key) == supported_keys.end())
                throw runtime_error("Error parsing PAL file: Unsupported PAL key '" + key + "' on line " + std::to_string(line_number));

            const auto duplicate = metadata.find(key);
            if (duplicate != metadata.end()) {
                throw runtime_error("Error parsing PAL file: Duplicate key '" + key + "' on line " + std::to_string(line_number));
            }

            metadata[key] = value;
        }
    }
    sequence_count = 0;
    letters = 0;
    oid_index.push_back(0);
    
	for (vector<string>::iterator it = volumes.begin(); it != volumes.end();) {
        if (it->length() >= 2 && (*it)[0] == '"' && (*it)[it->size() - 1] == '"') {
            const string nested = it->substr(1, it->length() - 2);
            it = volumes.erase(it);
            it = recurse(is_absolute_path(nested) ? nested : db_dir + PATH_SEPARATOR + nested, it);
        }
        else {
            BlastVolume vol(*it, 0, 0, 0, false);
            sequence_count += vol.index().num_oids;
            oid_index.push_back(vol.index().num_oids + oid_index.back());
            letters += vol.index().total_length;
            version = vol.index().version;
            ++it;
        }
    }
    if (metadata.find("SEQIDLIST") != metadata.end()) {
        if (ends_with(metadata["SEQIDLIST"], ".bsl"))
            throw runtime_error("Binary SEQIDLIST files(.bsl) are not supported, use text file instead : " + metadata["SEQIDLIST"]);
        if (!is_absolute_path(metadata["SEQIDLIST"])) {
            metadata["SEQIDLIST"] = db_dir + PATH_SEPARATOR + metadata["SEQIDLIST"];
        }
    }
    if (metadata.find("TAXIDLIST") != metadata.end() && !is_absolute_path(metadata["TAXIDLIST"])) {
        metadata["TAXIDLIST"] = db_dir + PATH_SEPARATOR + metadata["TAXIDLIST"];
    }
    assert(sequence_count > 0);
}

int Pal::volume(OId oid) const noexcept {
	assert(oid < sequence_count);
	auto it = std::upper_bound(oid_index.begin(), oid_index.end(), oid);
    return int(it - oid_index.begin() - 1);
}