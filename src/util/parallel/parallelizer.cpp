#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <stdexcept>
#include <memory>
#include "../system.h"
#ifndef WIN32
#include <unistd.h>
#include <sys/stat.h>
#endif

#include "multiprocessing.h"
// #define DEBUG
#undef DEBUG
#include "filestack.h"
#include "parallelizer.h"

using namespace std;


std::shared_ptr<Parallelizer> Parallelizer::instance_ptr = nullptr;

std::shared_ptr<Parallelizer> Parallelizer::get() {
    DBG("");
    if (! instance_ptr) {
        instance_ptr = shared_ptr<Parallelizer>(new Parallelizer());
    }
    return instance_ptr;
}


Parallelizer::Parallelizer() : work_directory("parallelizer"), n_registered(0), master_flag(true), i_barrier(0), initialized(false) {
    DBG("");
    // call init() later explicitly for final inizialization including worker registration
}


void Parallelizer::init(const string & tempdir) {
    DBG("");
    if (tempdir.size() > 0) {
        work_directory = join_path(tempdir, work_directory);
    }
    DBG("work_directory = " + work_directory);

    errno = 0;
#ifndef WIN32
    int s = mkdir(work_directory.c_str(), 00770);
#else
    int s = 0;
#endif
    if (s != 0) {
        if (errno == EEXIST) {
            // directory did already exist
        } else {
            throw(runtime_error("could not create working directory " + work_directory + " for parallelizer"));
        }
    }

    char hostname[1024];
    hostname[1023] = '\0';
#ifdef WIN32    
#else
    gethostname(hostname, 1023);
	id = string(hostname) + "_" + to_string(getpid());
#endif
    DBG("id = " + id);

    create_stack(LOG, id);
    create_stack(COMMAND);
    create_stack(WORKERS);
    create_stack(REGISTER);

    barrier_file = join_path(work_directory, BARRIER);

    log("PARALLELIZER BEGIN");

    initialized = true;
}


void Parallelizer::clear() {
}

Parallelizer::~Parallelizer() {
    DBG("");
    if (initialized) {
        log("PARALLELIZER END");
        clean(continuous_cleanup_list);
        clean(final_cleanup_list);
    }
}

string Parallelizer::get_id() {
    return id;
}

int Parallelizer::get_rank() {
    return rank;
}

string Parallelizer::get_work_directory() {
    return work_directory;
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
    if (! initialized) {
        return false;
    }

    auto cmd_file_name = get_barrier_file_name("cmd", tag, i_barrier);
    DBG(cmd_file_name);
    auto cmd_fs = FileStack(cmd_file_name);
    auto ack_file_name = get_barrier_file_name("ack", tag, i_barrier);
    DBG(ack_file_name);
    auto ack_fs = FileStack(ack_file_name);

    static const string msg = "WAIT";
    if (is_master()) {
        ack_fs.clear();
        cmd_fs.push(msg);
    }
    cmd_fs.poll_query(msg);
    ack_fs.push(id);

    DBG(msg);

    static const string msg_ok = "GOON";
    if (is_master()) {
        const size_t n_workers = get_stack(WORKERS)->size();
        ack_fs.poll_size(n_workers);
        cmd_fs.push(msg_ok);
    }
    cmd_fs.poll_query(msg_ok);

    DBG(msg_ok);

    if (is_master()) {
        clean(continuous_cleanup_list);
        continuous_cleanup_list.push_back(cmd_file_name);
        continuous_cleanup_list.push_back(ack_file_name);
    }

    i_barrier++;
    return true;
}


bool Parallelizer::register_workers(const double sleep_s) {
    DBG(id);
    get_stack(REGISTER)->push(id);
    sleep(sleep_s);
    if (is_master()) {
        string line;
        while (get_stack(REGISTER)->pop(line)) {
            get_stack(WORKERS)->push(line);
            n_registered++;
        }
        DBG("n_registered = " + to_string(n_registered));
    }
    return true;
}


bool Parallelizer::create_stack(const std::string & tag, std::string sfx) {
    if (fs_map.find(tag) == fs_map.end()) {
        if (sfx.size() > 0) {
            sfx = "_" + sfx;
        }
        const string file_name = join_path(work_directory, tag + sfx);
        DBG(file_name);
        return create_stack_from_file(tag, file_name);
    } else {
        return false;
    }
}


bool Parallelizer::create_stack_from_file(const std::string & tag, const std::string & file_name) {
    delete_stack(tag);
    fs_map.emplace(tag, shared_ptr<FileStack>(new FileStack(file_name)));
    DBG(file_name);
    return true;
}


std::shared_ptr<FileStack> Parallelizer::get_stack(const std::string & tag) {
	// cerr << __PRETTY_FUNCTION__ << endl;
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
        UNLINK(s.c_str());
    }
    file_list.clear();
    return true;
}


void Parallelizer::list_filestacks() {
	// cerr << __PRETTY_FUNCTION__ << endl;
    for (auto item : fs_map)
        cerr << item.first << " : " << item.second << endl;
}


void Parallelizer::log(const string & buf) {
	auto log_stack = get_stack(LOG);

    // we use ms since the Epoch as the universal timestamp
    const auto ms = chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now() - chrono::time_point<chrono::system_clock>{}
        ).count();
    const string tagged_buf = to_string(ms) + ' ' + buf + '\n';

    log_stack->push_non_locked(tagged_buf);
}
