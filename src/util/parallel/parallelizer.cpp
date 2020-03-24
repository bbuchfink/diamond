#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <stdexcept>
#include <memory>
#include <unistd.h>
#include <sys/stat.h>

#include "filestack.h"
#include "parallelizer.h"

using namespace std;


#undef DEBUG
#ifdef DEBUG
#define DBG(x) (std::cerr << x << std::endl)
#else
#define DBG(x)
#endif


string join_path(const string & path_1, const string & path_2) {
// #ifdef WINDOWS
//     const string sep = "\";
// #else
    const string sep = "/";
// #endif
    return path_1 + sep + path_2;
}



Parallelizer::Parallelizer() : work_directory("libworkstack"), n_registered(0), i_barrier(0) {
    {
        char* env_str = std::getenv("TMPDIR");
        if (env_str) {
            work_directory = join_path(string(env_str), work_directory);
        }
        env_str = std::getenv("SLURM_JOBID");
        if (env_str) {
            work_directory = work_directory + "_" + string(env_str);
        }
    }

    {
        errno = 0;
        int s = mkdir(work_directory.c_str(), 00770);
        if (s != 0) {
            if (errno == EEXIST) {
                // directory did already exist
            } else {
                throw(runtime_error("could not create working directory " + work_directory + " for parallelizer"));
            }
        }
    }

    {
        int rank = -1;
        const vector<string> env_opts = {"SLURM_PROCID", "PARALLEL_RANK"};
        for (auto env : env_opts) {
            const char* env_str = std::getenv(env.c_str());
            if (env_str) {
                rank = std::stoi(env_str);
                break;
            }
        }
        if (rank < 0) {
            const string msg = "parallel: Could not determine the parallel rank. Please set it via one of the environment variables "
                + env_opts[0] + ", " + env_opts[1] + ".";
            throw std::runtime_error(msg);
        }
        if (rank == 0) {
            master_flag = true;
        } else {
            master_flag = false;
        }
        // TODO: use hostname and process id, alternatively
        id = "rank_" + to_string(rank);
    }

    create_stack(LOG, id);
    create_stack(COMMAND);
    create_stack(WORKERS);
    create_stack(REGISTER);

    barrier_file = join_path(work_directory, BARRIER);

    get_stack(LOG)->clear();
    if (is_master()) {
        get_stack(COMMAND)->clear();
        get_stack(WORKERS)->clear();
        get_stack(REGISTER)->clear();
        sleep(1.0);
    } else {
        sleep(1.0);
    }

    register_workers();
}

Parallelizer::~Parallelizer() {
    clean(continuous_cleanup_list);
    // clean(final_cleanup_list);
}

string Parallelizer::get_id() {
    return id;
}

int Parallelizer::get_n_registered() {
    return n_registered;
}

bool Parallelizer::is_master() {
    return master_flag;
}

string Parallelizer::get_barrier_file_name(const string & step, const string & tag, int i) {
    return barrier_file + "_" + step + "_" + tag + "_" + to_string(i);
}

bool Parallelizer::barrier(const string & tag) {
    auto cmd_file_name = get_barrier_file_name("cmd", tag, i_barrier);
    auto cmd_fs = FileStack(cmd_file_name);
    auto ack_file_name = get_barrier_file_name("ack", tag, i_barrier);
    auto ack_fs = FileStack(ack_file_name);

    static const string msg = "WAIT";
    if (is_master()) {
        ack_fs.clear();
        cmd_fs.push(msg);
    }
    cmd_fs.poll_query(msg);
    ack_fs.push(id);

    static const string msg_ok = "GOON";
    if (is_master()) {
        const size_t n_workers = get_stack(WORKERS)->size();
        ack_fs.poll_size(n_workers);
        cmd_fs.push(msg_ok);
    }
    cmd_fs.poll_query(msg_ok);

    if (is_master()) {
        clean(continuous_cleanup_list);
        continuous_cleanup_list.push_back(cmd_file_name);
        continuous_cleanup_list.push_back(ack_file_name);
    }

    i_barrier++;
    return true;
}


bool Parallelizer::register_workers(const double sleep_s) {
    get_stack(REGISTER)->push(id);
    sleep(sleep_s);
    if (is_master()) {
        string line;
        while (get_stack(REGISTER)->pop(line)) {
            get_stack(WORKERS)->push(line);
            n_registered++;
        }
    }
    return true;
}


bool Parallelizer::create_stack(const std::string & tag, std::string sfx) {
    if (fs_map.find(tag) == fs_map.end()) {
        if (sfx.size() > 0) {
            sfx = "_" + sfx;
        }
        const string file_name = join_path(work_directory, tag + sfx);
        fs_map.emplace(tag, shared_ptr<FileStack>(new FileStack(file_name)));
        final_cleanup_list.push_back(file_name);
        return true;
    } else {
        return false;
    }
}


std::shared_ptr<FileStack> Parallelizer::get_stack(const std::string & tag) {
    return fs_map.at(tag);
}


bool Parallelizer::delete_stack(const std::string & tag) {
    if (fs_map.find(tag) != fs_map.end()) {
        fs_map.erase(tag);
        return true;
    } else {
        return false;
    }
}


void Parallelizer::sleep(const double sleep_s) {
    const chrono::duration<double> sleep_time(sleep_s);
    this_thread::sleep_for(sleep_time);
}


bool Parallelizer::clean(vector<string> & file_list) {
    for (auto s : file_list) {
        errno = 0;
        unlink(s.c_str());
    }
    file_list.clear();
    return true;
}
