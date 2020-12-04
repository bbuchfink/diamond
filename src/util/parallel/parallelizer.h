#ifndef _PARALLELIZER_H_
#define _PARALLELIZER_H_

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <algorithm>
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

#endif
