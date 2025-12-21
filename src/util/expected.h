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

#pragma once
#include <string>
#include <exception>

namespace Diamond {

    enum struct ErrorCode {
        PAL_OPEN_ERROR, PAL_PARSE_ERROR, BLASTDB_UNSUPPORTED_FEATURE
    };

    struct Error {
        ErrorCode code;
        std::string message;
    };

    template<typename E>
    class Unexpected
    {
    public:
        explicit Unexpected(const E& error) : error_(error)
        {
        }

        explicit Unexpected(E&& error) : error_(static_cast<E&&>(error))
        {
        }

        E& error()
        {
            return error_;
        }

        const E& error() const
        {
            return error_;
        }

    private:
        E error_;
    };

    class BadExpectedAccess : public std::exception
    {
    public:
        const char* what() const noexcept override
        {
            return "bad expected access";
        }
    };

    template<typename T, typename E>
    class Expected
    {
    public:
        Expected() : has_value_(true)
        {
            new (&storage_.value_) T();
        }

        Expected(const T& value) : has_value_(true)
        {
            new (&storage_.value_) T(value);
        }

        Expected(T&& value) : has_value_(true)
        {
            new (&storage_.value_) T(static_cast<T&&>(value));
        }

        Expected(const Unexpected<E>& error) : has_value_(false)
        {
            new (&storage_.error_) E(error.error());
        }

        Expected(Unexpected<E>&& error) : has_value_(false)
        {
            new (&storage_.error_) E(static_cast<E&&>(error.error()));
        }

        Expected(const Expected& other) : has_value_(other.has_value_)
        {
            if (has_value_)
                new (&storage_.value_) T(other.storage_.value_);
            else
                new (&storage_.error_) E(other.storage_.error_);
        }

        Expected(Expected&& other) : has_value_(other.has_value_)
        {
            if (has_value_)
                new (&storage_.value_) T(static_cast<T&&>(other.storage_.value_));
            else
                new (&storage_.error_) E(static_cast<E&&>(other.storage_.error_));
        }

        ~Expected()
        {
            destroy();
        }

        Expected& operator=(const Expected& other)
        {
            if (this == &other)
                return *this;
            assign_copy(other);
            return *this;
        }

        Expected& operator=(Expected&& other)
        {
            if (this == &other)
                return *this;
            assign_move(static_cast<Expected&&>(other));
            return *this;
        }

        bool has_value() const
        {
            return has_value_;
        }

        explicit operator bool() const
        {
            return has_value_;
        }

        T& value()
        {
            if (!has_value_)
                throw BadExpectedAccess();
            return storage_.value_;
        }

        const T& value() const
        {
            if (!has_value_)
                throw BadExpectedAccess();
            return storage_.value_;
        }

        constexpr const T* operator->() const noexcept {
            return &storage_.value_;
        }

        constexpr T* operator->() noexcept {
            return &storage_.value_;
        }

        E& error()
        {
            if (has_value_)
                throw BadExpectedAccess();
            return storage_.error_;
        }

        const E& error() const
        {
            if (has_value_)
                throw BadExpectedAccess();
            return storage_.error_;
        }

        template<typename U>
        T value_or(U&& default_value) const
        {
            if (has_value_)
                return storage_.value_;
            return T(static_cast<U&&>(default_value));
        }

        void swap(Expected& other)
        {
            if (has_value_ && other.has_value_)
            {
                using std::swap;
                swap(storage_.value_, other.storage_.value_);
            }
            else if (has_value_ && !other.has_value_)
            {
                T temp(static_cast<T&&>(storage_.value_));
                storage_.value_.~T();
                new (&storage_.error_) E(static_cast<E&&>(other.storage_.error_));
                other.storage_.error_.~E();
                new (&other.storage_.value_) T(static_cast<T&&>(temp));
                other.has_value_ = true;
                has_value_ = false;
            }
            else if (!has_value_ && other.has_value_)
            {
                other.swap(*this);
            }
            else
            {
                using std::swap;
                swap(storage_.error_, other.storage_.error_);
            }
        }

    private:
        struct storage
        {
            storage()
            {
            }

            ~storage()
            {
            }

            union
            {
                T value_;
                E error_;
            };
        } storage_;

        bool has_value_;

        void destroy()
        {
            if (has_value_)
                storage_.value_.~T();
            else
                storage_.error_.~E();
        }

        void assign_copy(const Expected& other)
        {
            if (has_value_ && other.has_value_)
            {
                storage_.value_ = other.storage_.value_;
            }
            else if (has_value_ && !other.has_value_)
            {
                storage_.value_.~T();
                new (&storage_.error_) E(other.storage_.error_);
                has_value_ = false;
            }
            else if (!has_value_ && other.has_value_)
            {
                storage_.error_.~E();
                new (&storage_.value_) T(other.storage_.value_);
                has_value_ = true;
            }
            else
            {
                storage_.error_ = other.storage_.error_;
            }
        }

        void assign_move(Expected&& other)
        {
            if (has_value_ && other.has_value_)
            {
                storage_.value_ = static_cast<T&&>(other.storage_.value_);
            }
            else if (has_value_ && !other.has_value_)
            {
                storage_.value_.~T();
                new (&storage_.error_) E(static_cast<E&&>(other.storage_.error_));
                has_value_ = false;
            }
            else if (!has_value_ && other.has_value_)
            {
                storage_.error_.~E();
                new (&storage_.value_) T(static_cast<T&&>(other.storage_.value_));
                has_value_ = true;
            }
            else
            {
                storage_.error_ = static_cast<E&&>(other.storage_.error_);
            }
        }
    };

    template<typename E>
    Unexpected<E> make_unexpected(E error)
    {
        return Unexpected<E>(static_cast<E&&>(error));
    }

}