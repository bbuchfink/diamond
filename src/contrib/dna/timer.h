/****
DIAMOND protein aligner
Copyright (C) 2022 Dimitrios Koutsogiannis

Code developed by Dimitrios Koutsogiannis <dimitrios.koutsogiannis@tue.mpg.de>

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

#include "chrono"
#include "mutex"
#include "iostream"


namespace Dna{
    using Duration = std::chrono::duration<long, std::ratio<1, 1000000000>>;

struct ExtensionTimer {
    ExtensionTimer():
            total_time(),
            preprocessing_time(),
            postprocessing_time(),
            extension(),
            next(),
            chaining(),
            extension_filter()
    {};
    ExtensionTimer& operator+=(ExtensionTimer& other)noexcept{
        this->total_time += (other.total_time);
        this->extension += (other.extension);
        this->preprocessing_time += (other.preprocessing_time);
        this->postprocessing_time += (other.postprocessing_time);
        this->next += (other.next);
        this->chaining += (other.chaining);
        this->extension_filter += (other.extension_filter);

        other.total_time =  std::chrono::duration<int>(0);
        other.extension =  std::chrono::duration<int>(0);
        other.preprocessing_time =  std::chrono::duration<int>(0);
        other.postprocessing_time =  std::chrono::duration<int>(0);
        other.next =  std::chrono::duration<int>(0);
        other.chaining = std::chrono::duration<int>(0);
        other.extension_filter = std::chrono::duration<int>(0);

        return *this;
    };
    void update(int operation , std::chrono::duration<long, std::ratio<1, 1000000000>> duration){
        switch(operation){
            case 0:
                total_time.operator+=(duration);
                break;
            case 1:
                preprocessing_time.operator+=(duration);
                break;
            case 2:
                postprocessing_time.operator+=(duration);
                break;
            case 3:
                extension_filter.operator+=(duration);
                break;
            case 4:
                extension.operator+=(duration);
                break;
            case 5:
                next.operator+=(duration);
                break;
            case 6:
                chaining.operator+=(duration);
                break;
        }
    };
    Duration total_time;
    Duration preprocessing_time;
    Duration postprocessing_time;
    Duration extension;
    Duration next;
    Duration chaining;
    Duration extension_filter;

    std::mutex mtx;
};

struct TotalTime: ExtensionTimer{
    void print() {
        //std::cerr << "Total Time: " << total_time.count() / 1'000'000'000.0 <<  " seconds\n";
        //std::cerr << "Pre-Processing: " << preprocessing_time.count() / 1'000'000'000.0 <<  " seconds\n";
        std::cerr << "Chaining: " << chaining.count() / 1'000'000'000.0 <<  " seconds\n";
        std::cerr << "Extension-Time: " <<  extension.count() / 1'000'000'000.0 << " seconds\n";
        std::cerr << "seed-lookup: " <<  next.count() / 1'000'000'000.0 << " seconds\n";
        std::cerr << "build Hsp from cigar (part of extension): " << postprocessing_time.count() / 1'000'000'000.0 << " seconds\n";
        std::cerr << "ungapped: " << preprocessing_time.count() / 1'000'000'000.0 << " seconds\n";
        //std::cerr << "extension_filter: " << extension_filter.count() / 1'000'000'000.0 << " seconds\n";
    }
    ~TotalTime(){
        print();

    };
};
}
