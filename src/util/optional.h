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

#pragma once
#include <new>
#include <utility>
#include <type_traits>
#include <stdexcept>

struct in_place_t { explicit in_place_t() {} };
static const in_place_t in_place{};

struct nullopt_t { explicit nullopt_t(int) {} };
static const nullopt_t nullopt{ 0 };

struct bad_optional_access : std::logic_error {
  bad_optional_access() : std::logic_error("bad optional access") {}
};

template <class T>
class optional {
  static_assert(!std::is_reference<T>::value, "optional<T&> is not supported");
  static_assert(!std::is_void<T>::value,      "optional<void> is not supported");

public:
  typedef T value_type;

  optional() noexcept : engaged_(false) {}
  optional(nullopt_t) noexcept : engaged_(false) {}

  optional(const T& v) : engaged_(false) { emplace(v); }
  optional(T&& v)      : engaged_(false) { emplace(std::move(v)); }

  optional(in_place_t) : engaged_(false) { emplace(); }

  template <class... Args>
  explicit optional(in_place_t, Args&&... args) : engaged_(false) {
    emplace(std::forward<Args>(args)...);
  }

  optional(const optional& other) : engaged_(false) {
    if (other.engaged_) emplace(*other.ptr());
  }

  optional(optional&& other) noexcept(std::is_nothrow_move_constructible<T>::value)
    : engaged_(false) {
    if (other.engaged_) emplace(std::move(*other.ptr()));
  }

  ~optional() { reset(); }

  optional& operator=(nullopt_t) noexcept {
    reset();
    return *this;
  }

  optional& operator=(const optional& other) {
    if (this == &other) return *this;

    if (engaged_ && other.engaged_) {
      *ptr() = *other.ptr();
    } else if (engaged_ && !other.engaged_) {
      reset();
    } else if (!engaged_ && other.engaged_) {
      emplace(*other.ptr());
    }
    return *this;
  }

  optional& operator=(optional&& other) noexcept(
      std::is_nothrow_move_constructible<T>::value &&
      std::is_nothrow_move_assignable<T>::value) {
    if (this == &other) return *this;

    if (engaged_ && other.engaged_) {
      *ptr() = std::move(*other.ptr());
    } else if (engaged_ && !other.engaged_) {
      reset();
    } else if (!engaged_ && other.engaged_) {
      emplace(std::move(*other.ptr()));
    }
    return *this;
  }

  optional& operator=(const T& v) {
    if (engaged_) *ptr() = v;
    else emplace(v);
    return *this;
  }

  optional& operator=(T&& v) {
    if (engaged_) *ptr() = std::move(v);
    else emplace(std::move(v));
    return *this;
  }

  bool has_value() const noexcept { return engaged_; }
  explicit operator bool() const noexcept { return engaged_; }

  T& value() & {
    if (!engaged_) throw bad_optional_access();
    return *ptr();
  }
  const T& value() const & {
    if (!engaged_) throw bad_optional_access();
    return *ptr();
  }
  T&& value() && {
    if (!engaged_) throw bad_optional_access();
    return std::move(*ptr());
  }

  T* operator->() { return &value(); }
  const T* operator->() const { return &value(); }

  T& operator*() & { return value(); }
  const T& operator*() const & { return value(); }
  T&& operator*() && { return std::move(value()); }

  template <class U>
  T value_or(U&& default_value) const {
    return engaged_ ? *ptr() : static_cast<T>(std::forward<U>(default_value));
  }

  void reset() noexcept {
    if (engaged_) {
      ptr()->~T();
      engaged_ = false;
    }
  }

  template <class... Args>
  T& emplace(Args&&... args) {
    reset();
    new (&storage_) T(std::forward<Args>(args)...);
    engaged_ = true;
    return *ptr();
  }

  void swap(optional& other) noexcept(
      std::is_nothrow_move_constructible<T>::value &&
      std::is_nothrow_move_assignable<T>::value) {
    if (engaged_ && other.engaged_) {
      using std::swap;
      swap(*ptr(), *other.ptr());
    } else if (engaged_ && !other.engaged_) {
      other.emplace(std::move(*ptr()));
      reset();
    } else if (!engaged_ && other.engaged_) {
      emplace(std::move(*other.ptr()));
      other.reset();
    }
  }

private:
  T* ptr() { return reinterpret_cast<T*>(&storage_); }
  const T* ptr() const { return reinterpret_cast<const T*>(&storage_); }
  struct Storage {
      alignas(T) unsigned char storage[sizeof(T)];
  };
  Storage storage_;
  bool engaged_;
};

template <class T>
inline void swap(optional<T>& a, optional<T>& b) noexcept(noexcept(a.swap(b))) {
  a.swap(b);
}

template <class T>
inline bool operator==(const optional<T>& a, const optional<T>& b) {
  if (a.has_value() != b.has_value()) return false;
  return a.has_value() ? (*a == *b) : true;
}
template <class T>
inline bool operator!=(const optional<T>& a, const optional<T>& b) {
  return !(a == b);
}

template <class T>
inline bool operator==(const optional<T>& a, nullopt_t) { return !a.has_value(); }
template <class T>
inline bool operator==(nullopt_t, const optional<T>& a) { return !a.has_value(); }
template <class T>
inline bool operator!=(const optional<T>& a, nullopt_t) { return a.has_value(); }
template <class T>
inline bool operator!=(nullopt_t, const optional<T>& a) { return a.has_value(); }

template <class T>
inline optional<typename std::decay<T>::type> make_optional(T&& v) {
    return optional<typename std::decay<T>::type>(std::forward<T>(v));
}