/****
DIAMOND protein aligner
Copyright (C) 2019-2024 Max Planck Society for the Advancement of Science e.V.

Code developed by Klaus Reuter

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

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include "filestack.h"

#define AUTOTAG std::string(__FUNCTION__)+"_"+std::to_string(__LINE__)

class Parallelizer {
    public:
        Parallelizer();
        ~Parallelizer();

        void init(const std::string & tempdir = std::string());
        void clear();

        int get_rank();
        std::string get_id();
        std::string get_work_directory();
        int get_n_registered();

        bool is_master();
        bool register_workers(const double sleep_s=5.0);

        bool barrier(const std::string & tag);

        bool create_stack_from_file(const std::string & tag, const std::string & file_name);

        bool create_stack(const std::string & tag, std::string sfx="");
        bool delete_stack(const std::string & tag);
        std::shared_ptr<FileStack> get_stack(const std::string & tag);

        static void sleep(const double sleep_s);

        void log(const std::string & buf);

        const std::string LOG = "log";
        const std::string COMMAND = "command";
        const std::string WORKERS = "workers";
        const std::string REGISTER = "register";

        static std::shared_ptr<Parallelizer> get();

        void list_filestacks();

    private:
	    static std::shared_ptr<Parallelizer> instance_ptr;

        const std::string BARRIER = "barrier";
        std::string work_directory;
        std::string barrier_file;

        int rank;
        std::string id;
        int n_registered;
        bool master_flag;
        int i_barrier;
        bool initialized;

        std::string get_barrier_file_name(const std::string & step, const std::string & tag, int i);

        bool clean(std::vector<std::string> & file_list);

        std::vector<std::string> continuous_cleanup_list;
        std::vector<std::string> final_cleanup_list;

        std::unordered_map<std::string, std::shared_ptr<FileStack>> fs_map;
};