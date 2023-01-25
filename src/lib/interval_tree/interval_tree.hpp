#pragma once

#include "interval_tree_fwd.hpp"
#include "interval_types.hpp"

#include <string>
#include <memory>
#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <iterator>
#include <type_traits>

#include <iostream>

namespace lib_interval_tree
{
//############################################################################################################
    enum class rb_color
    {
        fail,
        red,
        black,
        double_black
    };
//############################################################################################################
    using default_interval_value_type = int;
//############################################################################################################
    template <typename numerical_type, typename interval_kind_ = closed>
    struct interval
    {
    public:
        using value_type = numerical_type;
        using interval_kind = interval_kind_;

        /**
         *  Constructs an interval. low MUST be smaller than high.
         */
#ifndef INTERVAL_TREE_SAFE_INTERVALS
#if __cplusplus >= 201703L
        constexpr
#endif
        interval(value_type low, value_type high)
            : low_{low}
            , high_{high}
        {
            assert(low <= high);
        }
#else
#if __cplusplus >= 201703L
        constexpr
#endif
        interval(value_type low, value_type high)
            : low_{std::min(low, high)}
            , high_{std::max(low, high)}
        {
        }
#endif

        /**
         *  Returns if both intervals equal.
         */
        friend bool operator==(interval const& lhs, interval const& other)
        {
            return lhs.low_ == other.low_ && lhs.high_ == other.high_;
        }

        /**
         *  Returns if both intervals are different.
         */
        friend bool operator!=(interval const& lhs, interval const& other)
        {
            return lhs.low_ != other.low_ || lhs.high_ != other.high_;
        }

        /**
         *  Returns the lower bound of the interval
         */
        value_type low() const
        {
            return low_;
        }

        /**
         *  Returns the upper bound of the interval
         */
        value_type high() const
        {
            return high_;
        }

        /**
         *  Returns whether the intervals overlap.
         *  For when both intervals are closed.
         */
        bool overlaps(value_type l, value_type h) const
        {
            return low_ <= h && l <= high_;
        }

        /**
         *  Returns whether the intervals overlap, excluding border.
         *  For when at least one interval is open (l&r).
         */
        bool overlaps_exclusive(value_type l, value_type h) const
        {
            return low_ < h && l < high_;
        }

        /**
         *  Returns whether the intervals overlap
         */
        bool overlaps(interval const& other) const
        {
            return overlaps(other.low_, other.high_);
        }

        /**
         *  Returns whether the intervals overlap, excluding border.
         */
        bool overlaps_exclusive(interval const& other) const
        {
            return overlaps_exclusive(other.low_, other.high_);
        }

        /**
         *  Returns whether the given value is in this.
         */
        bool within(value_type value) const
        {
            return interval_kind::within(low_, high_, value);
        }

        /**
         *  Returns whether the given interval is in this.
         */
        bool within(interval const& other) const
        {
            return low_ <= other.low_ && high_ >= other.high_;
        }

        /**
         *  Calculates the distance between the two intervals.
         *  Overlapping intervals have 0 distance.
         */
        value_type operator-(interval const& other) const
        {
            if (overlaps(other))
                return 0;
            if (high_ < other.low_)
                return other.low_ - high_;
            else
                return low_ - other.high_;
        }

        /**
         *  Returns the size of the interval.
         */
        value_type size() const
        {
            return high_ - low_;
        }

        /**
         *  Creates a new interval from this and other, that contains both intervals and whatever
         *  is between.
         */
        interval join(interval const& other) const
        {
            return {std::min(low_, other.low_), std::max(high_, other.high_)};
        }

    private:
        value_type low_;
        value_type high_;
    };
//############################################################################################################
    /**
     *  Creates a safe interval that puts the lower bound left automatically.
     */
    template <typename numerical_type, typename interval_kind_ = closed>
#if __cplusplus >= 201703L
    constexpr
#endif
    interval <numerical_type, interval_kind_> make_safe_interval(numerical_type lhs, numerical_type rhs)
    {
        return interval <numerical_type, interval_kind_>{std::min(lhs, rhs), std::max(lhs, rhs)};
    }
//############################################################################################################
    template <typename numerical_type = default_interval_value_type, typename interval_type_ = interval <numerical_type, closed>>
    class node
    {
    private:
        using node_type = node <numerical_type, interval_type_>;

    public:
        using interval_type = interval_type_;
        using value_type = numerical_type;

    public:
        friend lib_interval_tree::interval_tree <interval_type>;
        friend lib_interval_tree::const_interval_tree_iterator <node <numerical_type, interval_type> >;
        friend lib_interval_tree::interval_tree_iterator <node <numerical_type, interval_type> >;

    public:
        node(node* parent, interval_type interval)
            : interval_{std::move(interval)}
            , max_{interval.high()}
            , parent_{parent}
            , left_{}
            , right_{}
            , color_{rb_color::fail}
        {
        }

        ~node()
        {
        }

        interval_type interval() const
        {
            return interval_;
        }

        value_type max() const
        {
            return max_;
        }

        bool is_left() const noexcept
        {
            return this == parent_->left_;
        }

        bool is_right() const noexcept
        {
            return this == parent_->right_;
        }

        bool is_root() const noexcept
        {
            return !parent_;
        }

        /**
         *  Returns the color of the node.
         */
        rb_color color() const
        {
            return color_;
        }

        /**
         *  Returns the parent node up the tree.
         */
        node const* parent() const
        {
            return parent_;
        }

        /**
         *  Returns the left node (readonly).
         */
        node const* left() const
        {
            return left_;
        }

        /**
         *  Returns the right node (readonly).
         */
        node const* right() const
        {
            return right_;
        }

        /**
         *  Returns the height of the node in the tree. Where height = how many parents does it have.
         *  The root has no parents and is therefor has height 0.
         */
        int height() const
        {
            int counter{0};
            for (auto* p = parent_; p != nullptr; p = p->parent_)
                ++counter;
            return counter;
        }

        /**
         *  Returns the lower bound of the interval of this node
         */
        value_type low() const
        {
            return interval_.low();
        }

        /**
         *  Returns the upper bound of the interval of this node
         */
        value_type high() const
        {
            return interval_.high();
        }

private:
        void set_interval(interval_type const& ival)
        {
            interval_ = ival;
        }

        void kill() const noexcept
        {
            auto* parent = parent_;
            if (is_left())
            {
                delete parent_->left_;
                parent->left_ = nullptr;
            }
            else
            {
                delete parent_->right_;
                parent->right_ = nullptr;
            }
        }

    private:
        interval_type interval_;
        value_type max_;
        node* parent_;
        node* left_;
        node* right_;
        rb_color color_;
    };
//############################################################################################################
    template <typename node_type, typename owner_type>
    class basic_interval_tree_iterator : public std::forward_iterator_tag
    {
    public:
        friend interval_tree <typename node_type::interval_type>;

        using tree_type = interval_tree <typename node_type::interval_type>;
        using value_type = node_type;

        using node_ptr_t = typename std::conditional <
            std::is_const <typename std::remove_pointer <owner_type>::type>::value,
            node_type const*,
            node_type*
        >::type;

    public:
        constexpr basic_interval_tree_iterator(basic_interval_tree_iterator const&) = default;
        basic_interval_tree_iterator& operator=(basic_interval_tree_iterator const&) = default;

        bool operator!=(basic_interval_tree_iterator const& other) const
        {
            return node_ != other.node_;
        }

        bool operator==(basic_interval_tree_iterator const& other) const
        {
            return node_ == other.node_;
        }

        /**
         *  Returns the max property of the node.
         */
        typename node_type::interval_type::value_type max() const
        {
            return node_->max();
        }

        /**
         *  Returns the color of the node.
         */
        rb_color color() const
        {
            return node_->color();
        }

        typename tree_type::interval_type interval() const
        {
            return node_->interval();
        }

        virtual ~basic_interval_tree_iterator() = default;

    protected:
        basic_interval_tree_iterator(node_ptr_t node, owner_type owner)
            : node_{node}
            , owner_{owner}
        {
        }

    protected:
        node_ptr_t node_;
        owner_type owner_;
    };
//############################################################################################################
    template <typename node_type>
    class const_interval_tree_iterator
        : public basic_interval_tree_iterator <node_type,
                                               interval_tree <typename node_type::interval_type> const*>
    {
    public:
        using tree_type = interval_tree <typename node_type::interval_type>;
        using iterator_base = basic_interval_tree_iterator <node_type, tree_type const*>;
        using value_type = typename iterator_base::value_type;
        using iterator_base::node_;
        using iterator_base::owner_;

        friend tree_type;

    public:
        const_interval_tree_iterator& operator++()
        {
            if (!node_)
            {
                node_ = owner_->root_;

                if (!node_)
                    return *this;

                while(node_->left_)
                    node_ = node_->left_;
            }

            if (node_->right_)
            {
                node_ = node_->right_;

                while (node_->left_)
                    node_ = node_->left_;
            }
            else
            {
                auto* parent = node_->parent_;
                while (parent != nullptr && node_ == parent->right_)
                {
                    node_ = parent;
                    parent = parent->parent_;
                }
                node_ = parent;
            }

            return *this;
        }

        const_interval_tree_iterator operator++(int)
        {
            const_interval_tree_iterator cpy = *this;
            operator++();
            return cpy;
        }

        typename value_type::interval_type operator*() const
        {
            if (node_)
                return node_->interval();
            else
                throw std::out_of_range("dereferencing interval_tree_iterator out of bounds");
        }

        /**
         *  Returns an iterator to the parent of this node.
         *  will equal std::end(tree) if there is no parent node.
         */
        const_interval_tree_iterator parent() const
        {
            if (node_)
                return {node_->parent_, owner_};
            else
                throw std::out_of_range("interval_tree_iterator out of bounds");
        }

        /**
         *  Continues down the left side of this node.
         *  will equal std::end(tree) if there is no left node.
         */
        const_interval_tree_iterator left() const
        {
            if (node_)
                return {node_->left_, owner_};
            else
                throw std::out_of_range("interval_tree_iterator out of bounds");
        }

        /**
         *  Continues down the right side of this node.
         *  will equal std::end(tree) if there is no right node.
         */
        const_interval_tree_iterator right() const
        {
            if (node_)
                return {node_->right_, owner_};
            else
                throw std::out_of_range("interval_tree_iterator out of bounds");
        }

        value_type const* operator->() const
        {
            return node_;
        }

    private:
        const_interval_tree_iterator(node_type const* node, tree_type const* owner)
            : basic_interval_tree_iterator <node_type, tree_type const*> {node, owner}
        {
        }
    };
//############################################################################################################
    template <typename node_type>
    class interval_tree_iterator
        : public basic_interval_tree_iterator <node_type,
                                               interval_tree <typename node_type::interval_type>*>
    {
    public:
        using tree_type = interval_tree <typename node_type::interval_type>;
        using iterator_base = basic_interval_tree_iterator <node_type, tree_type*>;
        using value_type = typename iterator_base::value_type;
        using iterator_base::node_;
        using iterator_base::owner_;

        friend tree_type;

    public:
        interval_tree_iterator& operator++()
        {
            if (!node_)
            {
                node_ = owner_->root_;

                if (!node_)
                    return *this;

                while(node_->left_)
                    node_ = node_->left_;
            }

            if (node_->right_)
            {
                node_ = node_->right_;

                while (node_->left_)
                    node_ = node_->left_;
            }
            else
            {
                auto* parent = node_->parent_;
                while (parent != nullptr && node_ == parent->right_)
                {
                    node_ = parent;
                    parent = parent->parent_;
                }
                node_ = parent;
            }

            return *this;
        }

        interval_tree_iterator operator++(int)
        {
            interval_tree_iterator cpy = *this;
            operator++();
            return cpy;
        }

        /**
         *  Returns an iterator to the parent of this node.
         *  will equal std::end(tree) if there is no parent node.
         */
        interval_tree_iterator parent()
        {
            if (node_)
                return {node_->parent_, owner_};
            else
                throw std::out_of_range("interval_tree_iterator out of bounds");
        }

        /**
         *  Continues down the left side of this node.
         *  will equal std::end(tree) if there is no left node.
         */
        interval_tree_iterator left()
        {
            if (node_)
                return {node_->left_, owner_};
            else
                throw std::out_of_range("interval_tree_iterator out of bounds");
        }

        /**
         *  Continues down the right side of this node.
         *  will equal std::end(tree) if there is no right node.
         */
        interval_tree_iterator right()
        {
            if (node_)
                return {node_->right_, owner_};
            else
                throw std::out_of_range("interval_tree_iterator out of bounds");
        }

        typename value_type::interval_type operator*() const
        {
            if (node_)
                return node_->interval();
            else
                throw std::out_of_range("interval_tree_iterator out of bounds");
        }

        value_type* operator->()
        {
            return node_;
        }

    private:
        interval_tree_iterator(node_type* node, tree_type* owner)
            : basic_interval_tree_iterator <node_type, tree_type*> {node, owner}
        {
        }
    };
//############################################################################################################
    template <typename IntervalT = interval <int, closed>>
    class interval_tree
    {
    public:
        using interval_type = IntervalT;
        using value_type = typename interval_type::value_type;
        using node_type = node <value_type, interval_type>;
        using iterator = interval_tree_iterator <node_type>;
        using const_iterator = const_interval_tree_iterator <node_type>;
        using size_type = long long;
        using this_type = interval_tree<interval_type>;

    public:
        friend const_interval_tree_iterator <node_type>;
        friend interval_tree_iterator <node_type>;

        interval_tree()
            : root_{nullptr}
            , size_{0}
        {
        }

        ~interval_tree()
        {
            clear();
        }

        interval_tree(interval_tree const& other)
            : root_{nullptr}
            , size_{0}
        {
            operator=(other);
        }

    public:
        interval_tree& operator=(interval_tree const& other)
        {
            if (!empty())
                clear();

            if (other.root_ != nullptr)
                root_ = copyTreeImpl(other.root_, nullptr);

            size_ = other.size_;

            return *this;
        }

        /**
         *  Removes all from this tree.
         */
        void clear()
        {
            for (auto i = std::begin(*this); i != std::end(*this);)
                i = erase(i);
        }

        /**
         *  Returns the root node from this tree.
         */
        iterator root()
        {
            return {root_, this};
        }

        /**
         *  Returns the root node from this tree.
         */
        const_iterator root() const
        {
            return {root_, this};
        }

        /**
         *  Inserts an interval into the tree.
         */
        template <typename IntervalType = interval_type>
        iterator insert(IntervalType&& ival)
        {
            node_type* z = new node_type(nullptr, std::forward <IntervalType&&> (ival));
            node_type* y = nullptr;
            node_type* x = root_;
            while (x)
            {
                y = x;
                if (z->interval_.low() < x->interval_.low())
                    x = x->left_;
                else
                    x = x->right_;
            }
            z->parent_ = y;
            if (!y)
                root_ = z;
            else if (z->interval_.low() < y->interval_.low())
                y->left_ = z;
            else
                y->right_ = z;
            z->color_ = rb_color::red;

            insert_fixup(z);
            recalculate_max(z);
            ++size_;
            return {z, this};
        }

        /**
         *  Inserts an interval into the tree if no other interval overlaps it.
         *  Otherwise merge the interval with the being overlapped.
         *
         *  @param ival The interval
         *  @param exclusive Exclude borders.
         */
        iterator insert_overlap(interval_type const& ival, bool exclusive = false)
        {
            auto iter = overlap_find(ival, exclusive);
            if (iter == end())
                return insert(ival);
            else
            {
                auto merged = iter->interval().join(ival);
                erase(iter);
                return insert_overlap(merged);
            }
        }

        /**
         *  Erases the element pointed to be iter.
         */
        iterator erase(iterator iter)
        {
            if (!iter.node_)
                throw std::out_of_range("cannot erase end iterator");

            auto next = iter;
            ++next;

            node_type* y;
            if (!iter.node_->left_ || !iter.node_->right_)
                y = iter.node_;
            else
                y = successor(iter.node_);

            node_type* x;
            if (y->left_)
                x = y->left_;
            else
                x = y->right_;

            if (x)
                x->parent_ = y->parent_;

            auto* x_parent = y->parent_;

            if (!y->parent_)
                root_ = x;
            else if (y->is_left())
                y->parent_->left_ = x;
            else
                y->parent_->right_ = x;

            if (y != iter.node_)
            {
                iter.node_->interval_ = y->interval_;
                iter.node_->max_ = y->max_;
                recalculate_max(iter.node_);
            }

            if (x && x->color_ == rb_color::red)
            {
                if (x_parent)
                    erase_fixup(x, x_parent, y->is_left());
                else
                    x->color_ = rb_color::black;
            }

            delete y;

            --size_;
            return next;
        }

        /**
         *  Returns the size of the object.
         */
        size_type size() const
        {
            return size_;
        }

        /**
         *  Finds the first exact match.
         *
         *  @param ival The interval to find an exact match for within the tree.
         *  @param compare A comparison function to use.
         */
        template <typename CompareFunctionT>
        iterator find(interval_type const& ival, CompareFunctionT const& compare)
        {
            if (root_ == nullptr)
                return end();
            return iterator{find_i(root_, ival, compare), this};
        }

        /**
         *  Finds the first exact match.
         *
         *  @param ival The interval to find an exact match for within the tree.
         *  @param compare A comparison function to use.
         */
        template <typename CompareFunctionT>
        const_iterator find(interval_type const& ival, CompareFunctionT const& compare) const
        {
            if (root_ == nullptr)
                return end();
            return const_iterator{find_i(root_, ival, compare), this};
        }

        /**
         *  Finds the first exact match.
         *
         *  @param ival The interval to find an exact match for within the tree.
         */
        iterator find(interval_type const& ival)
        {
            return find(ival, [](auto const& lhs, auto const& rhs){return lhs == rhs;});
        }
        /**
         *  Finds the first exact match.
         *
         *  @param ival The interval to find an exact match for within the tree.
         */
        const_iterator find(interval_type const& ival) const
        {
            return find(ival, [](auto const& lhs, auto const& rhs){return lhs == rhs;});
        }

        /**
         *  Finds all exact matches and returns the amount of intervals found.
         */
        template <typename FunctionT, typename CompareFunctionT>
        void find_all(interval_type const& ival, FunctionT const& on_find, CompareFunctionT const& compare)
        {
            if (root_ == nullptr)
                return;
            find_all_i<this_type, iterator>(this, root_, ival, on_find, compare);
        }
        template <typename FunctionT, typename CompareFunctionT>
        void find_all(interval_type const& ival, FunctionT const& on_find, CompareFunctionT const& compare) const
        {
            if (root_ == nullptr)
                return;
            find_all_i<this_type, const_iterator>(this, root_, ival, on_find, compare);
        }

        template <typename FunctionT>
        void find_all(interval_type const& ival, FunctionT const& on_find)
        {
            find_all(ival, on_find, [](auto const& lhs, auto const& rhs){return lhs == rhs;});
        }
        template <typename FunctionT>
        void find_all(interval_type const& ival, FunctionT const& on_find) const
        {
            find_all(ival, on_find, [](auto const& lhs, auto const& rhs){return lhs == rhs;});
        }

        /**
         *  Finds the next exact match EXCLUDING from.
         *
         *  @param from The iterator to search from EXCLUSIVE!
         *  @param ival The interval to find an exact match for within the tree.
         *  @param compare A comparison function to use.
         */
        template <typename CompareFunctionT>
        iterator find_next_in_subtree(iterator from, interval_type const& ival, CompareFunctionT const& compare)
        {
            if (root_ == nullptr)
                return end();
            return iterator{find_i_ex(from.node_, ival, compare), this};
        }
        template <typename CompareFunctionT>
        const_iterator find_next_in_subtree(iterator from, interval_type const& ival, CompareFunctionT const& compare) const
        {
            if (root_ == nullptr)
                return end();
            return iterator{find_i_ex(from.node_, ival, compare), this};
        }

        /**
         *  Finds the next exact match EXCLUDING from.
         *
         *  @param from The iterator to search from, EXCLUSIVE!
         *  @param ival The interval to find an exact match for within the tree.
         *  @param compare A comparison function to use.
         */
        iterator find_next_in_subtree(iterator from, interval_type const& ival)
        {
            return find_next_in_subtree(from, ival, [](auto const& lhs, auto const& rhs){return lhs == rhs;});
        }
        const_iterator find_next_in_subtree(iterator from, interval_type const& ival) const
        {
            return find_next_in_subtree(from, ival, [](auto const& lhs, auto const& rhs){return lhs == rhs;});
        }

        /**
         *  Finds the first interval that overlaps with ival.
         *
         *  @param ival The interval to find an overlap for within the tree.
         *  @param exclusive Exclude edges?
         */
        iterator overlap_find(interval_type const& ival, bool exclusive = false)
        {
            if (root_ == nullptr)
                return end();
            if (exclusive)
                return iterator{overlap_find_i<true>(root_, ival), this};
            else
                return iterator{overlap_find_i<false>(root_, ival), this};
        }
        const_iterator overlap_find(interval_type const& ival, bool exclusive = false) const
        {
            if (root_ == nullptr)
                return end();
            if (exclusive)
                return const_iterator{overlap_find_i<true>(root_, ival), this};
            else
                return const_iterator{overlap_find_i<false>(root_, ival), this};
        }

        /**
         *  Finds all intervals that overlaps with ival.
         *
         *  @param ival The interval to find an overlap for within the tree.
         *  @param exclusive Exclude edges?
         */
        template <typename FunctionT>
        void overlap_find_all(interval_type const& ival, FunctionT& on_find, bool exclusive = false)
        {
            if (root_ == nullptr)
                return;
            if (exclusive)
                overlap_find_all_i<this_type, true, iterator>(this, root_, ival, on_find);
            else
                overlap_find_all_i<this_type, false, iterator>(this, root_, ival, on_find);
        }
        template <typename FunctionT>
        void overlap_find_all(interval_type const& ival, FunctionT& on_find, bool exclusive = false) const
        {
            if (root_ == nullptr)
                return;
            if (exclusive)
                overlap_find_all_i<this_type, true, const_iterator>(this, root_, ival, on_find);
            else
                overlap_find_all_i<this_type, false, const_iterator>(this, root_, ival, on_find);
        }

        /**
         *  Finds the next interval that overlaps with ival
         *
         *  @param from The iterator to start from, EXCLUSIVE!
         *  @param ival The interval to find an overlap for within the tree.
         *  @param exclusive Exclude edges?
         */
        iterator overlap_find_next_in_subtree(iterator from, interval_type const& ival, bool exclusive = false)
        {
            if (root_ == nullptr)
                return end();
            return iterator{overlap_find_i_ex(from.node_, ival, exclusive), this};
        }
        const_iterator overlap_find_next_in_subtree(const_iterator  from, interval_type const& ival, bool exclusive = false) const
        {
            if (root_ == nullptr)
                return end();
            return const_iterator {overlap_find_i_ex(from.node_, ival, exclusive), this};
        }

        /**
         *  Deoverlaps the tree but returns it as a copy.
         */
        interval_tree deoverlap_copy()
        {
            interval_tree fresh;
            for (auto i = begin(), e = end(); i != e; ++i)
                fresh.insert_overlap(*i);

            return fresh;
        }

        /**
         *  Merges all overlapping intervals by erasing overlapping intervals and reinserting the merged interval.
         */
        interval_tree& deoverlap()
        {
            *this = deoverlap_copy();
            return *this;
        }

        /**
         *  Only works with deoverlapped trees.
         *  Creates an interval tree that contains all gaps between the intervals as intervals.
         */
        interval_tree punch() const
        {
            if (empty())
                return {};
            auto min = std::begin(*this)->interval().low();
            auto max = root_->max_;
            return punch({min, max});
        }

        /**
         *  Only works with deoverlapped trees.
         *  Removes all intervals from the given interval and produces a tree that contains the remaining intervals.
         *  This is basically the other punch overload with ival = [tree_lowest, tree_highest]
         */
        interval_tree punch(interval_type const& ival) const
        {
            if (empty())
                return {};

            interval_tree result;
            auto i = std::begin(*this);
            if (ival.low() < i->interval().low())
                result.insert({ival.low(), i->interval().low()});

            for (auto e = end(); i != e; ++i)
            {
                auto next = i; ++next;
                if (next != e)
                    result.insert({i->interval().high(), next->interval().low()});
                else
                    break;
            }

            if (i != end() && i->interval().high() < ival.high())
                result.insert({i->interval().high(), ival.high()});

            return result;
        }

        iterator begin()
        {
            if (!root_)
                return {nullptr, this};

            auto* iter = root_;

            while (iter->left_)
                iter = iter->left_;

            return{iter, this};
        }
        iterator end()
        {
            return {nullptr, this};
        }

        const_iterator cbegin() const
        {
            if (!root_)
                return {nullptr, this};

            auto* iter = root_;

            while (iter->left_)
                iter = iter->left_;

            return const_iterator{iter, this};
        }
        const_iterator cend() const
        {
            return const_iterator{nullptr, this};
        }
        const_iterator begin() const
        {
            return cbegin();
        }
        const_iterator end() const
        {
            return cend();
        }

        /**
         *  Returns wether or not the tree is empty
         */
        bool empty() const noexcept
        {
            return root_ == nullptr;
        }

    private:
        node_type* copyTreeImpl(node_type* root, node_type* parent)
        {
            if (root)
            {
                auto* cpy = new node_type(parent, root->interval());
                cpy->color_ = root->color_;
                cpy->max_ = root->max_;
                cpy->left_ = copyTreeImpl(root->left_, cpy);
                cpy->right_ = copyTreeImpl(root->right_, cpy);
                return cpy;
            }
            return nullptr;
        };

        template <typename ThisType, typename IteratorT, typename FunctionT, typename ComparatorFunctionT>
        static bool find_all_i
        (
            typename std::conditional<std::is_same<IteratorT, iterator>::value, ThisType, ThisType const>::type* self,
            node_type* ptr, 
            interval_type const& ival, 
            FunctionT const& on_find, 
            ComparatorFunctionT const& compare
        )
        {
            if (compare(ptr->interval(), ival))
            {
                if (!on_find(IteratorT{ptr, self}))
                    return false;
            }
            if (ptr->left_ && ival.high() <= ptr->left_->max())
            {
                // no right? can only continue left
                if (!ptr->right_ || ival.low() > ptr->right_->max())
                    return find_all_i<ThisType, IteratorT>(self, ptr->left_, ival, on_find, compare);

                if (!find_all_i<ThisType, IteratorT>(self, ptr->left_, ival, on_find, compare))
                    return false;
            }
            if (ptr->right_ && ival.high() <= ptr->right_->max())
            {
                if (!ptr->left_ || ival.low() > ptr->left_->max())
                    return find_all_i<ThisType, IteratorT>(self, ptr->right_, ival, on_find, compare);

                if (!find_all_i<ThisType, IteratorT>(self, ptr->right_, ival, on_find, compare))
                    return false;
            }
            return true;
        }

        template <typename ComparatorFunctionT>
        node_type* find_i(node_type* ptr, interval_type const& ival, ComparatorFunctionT const& compare) const
        {
            if (compare(ptr->interval(), ival))
                return ptr;
            else
                return find_i_ex(ptr, ival, compare);
        }

        // excludes ptr
        template <typename ComparatorFunctionT>
        node_type* find_i_ex(node_type* ptr, interval_type const& ival, ComparatorFunctionT const& compare) const
        {
            if (ptr->left_ && ival.high() <= ptr->left_->max())
            {
                // no right? can only continue left
                if (!ptr->right_ || ival.low() > ptr->right_->max())
                    return find_i(ptr->left_, ival, compare);

                auto* res = find_i(ptr->left_, ival, compare);
                if (res != nullptr)
                    return res;
            }
            if (ptr->right_ && ival.high() <= ptr->right_->max())
            {
                if (!ptr->left_ || ival.low() > ptr->left_->max())
                    return find_i(ptr->right_, ival, compare);

                auto* res = find_i(ptr->right_, ival, compare);
                if (res != nullptr)
                    return res;
            }
            return nullptr;
        }

        template <bool Exclusive>
        node_type* overlap_find_i(node_type* ptr, interval_type const& ival) const
        {
#if __cplusplus >= 201703L
            if constexpr (Exclusive)
#else
            if (Exclusive)
#endif
            {
                if (ptr->interval().overlaps_exclusive(ival))
                    return ptr;
            }
            else
            {
                if (ptr->interval().overlaps(ival))
                    return ptr;
            }

            return overlap_find_i_ex<Exclusive>(ptr, ival);
        }

        template <typename ThisType, bool Exclusive, typename IteratorT, typename FunctionT>
        static bool overlap_find_all_i
        (
            typename std::conditional<std::is_same<IteratorT, iterator>::value, ThisType, ThisType const>::type* self,
            node_type* ptr, 
            interval_type const& ival, 
            FunctionT& on_find
        )
        {
#if __cplusplus >= 201703L
            if constexpr (Exclusive)
#else
            if (Exclusive)
#endif
            {
                if (ptr->interval().overlaps_exclusive(ival))
                {
                    if (!on_find(IteratorT{ptr, self}))
                    {
                        return false;
                    }
                }
            }
            else
            {
                if (ptr->interval().overlaps(ival))
                {
                    if (!on_find(IteratorT{ptr, self}))
                    {
                        return false;
                    }
                }
            }
            if (ptr->left_ && ptr->left_->max() >= ival.low())
            {
                // no right? can only continue left
                // or interval low is bigger than max of right branch.
                if (!ptr->right_ || ival.low() > ptr->right_->max())
                    return overlap_find_all_i<ThisType, Exclusive, IteratorT>(self, ptr->left_, ival, on_find);

                if (!overlap_find_all_i<ThisType, Exclusive, IteratorT>(self, ptr->left_, ival, on_find))
                    return false;
            }
            if (ptr->right_ && ptr->right_->max() >= ival.low())
            {
                if (!ptr->left_ || ival.low() > ptr->right_->max())
                    return overlap_find_all_i<ThisType, Exclusive, IteratorT>(self, ptr->right_, ival, on_find);

                if (!overlap_find_all_i<ThisType, Exclusive, IteratorT>(self, ptr->right_, ival, on_find))
                    return false;
            }
            return true;
        }

        // excludes ptr
        template <bool Exclusive>
        node_type* overlap_find_i_ex(node_type* ptr, interval_type const& ival) const
        {
            if (ptr->left_ && ptr->left_->max() >= ival.low())
            {
                // no right? can only continue left
                // or upper bounds higher than what is contained right? continue left.
                if (!ptr->right_ || ival.low() > ptr->right_->max())
                    return overlap_find_i<Exclusive>(ptr->left_, ival);

                auto* res = overlap_find_i<Exclusive>(ptr->left_, ival);
                if (res != nullptr)
                    return res;
            }
            if (ptr->right_ && ptr->right_->max() >= ival.low())
            {
                if (!ptr->left_ || ival.low() > ptr->left_->max())
                    return overlap_find_i<Exclusive>(ptr->right_, ival);

                auto* res = overlap_find_i<Exclusive>(ptr->right_, ival);
                if (res != nullptr)
                    return res;
            }
            return nullptr;
        }

        node_type* successor(node_type* node)
        {
            if (node->right_)
                return minimum(node->right_);
            auto* y = node->parent_;
            while (y && node == y->right_)
            {
                node = y;
                y = y->parent_;
            }
            return y;
        }

        bool is_descendant(iterator par, iterator desc)
        {
            auto p = desc->parent_;
            for (; p && p != par.node_; p = p->parent_) {}
            return p != nullptr;
        }

        /**
         *  Set v inplace of u. Does not delete u.
         *  Creates orphaned nodes. A transplant call must be succeeded by delete calls.
         */
        void transplant(node_type* u, node_type* v)
        {
            if (u->is_root())
                root_ = v;
            else if (u->is_left())
                u->parent_->left_ = v;
            else
                u->parent_->right_ = v;
            if (v)
                v->parent_ = u->parent_;
        }

        /**
         *  Get leftest of x.
         */
        node_type* minimum(node_type* x) const
        {
            while (x->left_)
                x = x->left_;
            return x;
        }

        void left_rotate(node_type* x)
        {
            auto* y = x->right_;
            x->right_ = y->left_;
            if (y->left_)
                y->left_->parent_ = x;

            y->parent_ = x->parent_;
            if (!x->parent_)
                root_ = y;
            else if (x->is_left())
                x->parent_->left_ = y;
            else
                x->parent_->right_ = y;

            y->left_ = x;
            x->parent_ = y;

            // max fixup
            if (x->left_ && x->right_)
                x->max_ = std::max(x->interval_.high(), std::max(x->left_->max_, x->right_->max_));
            else if (x->left_)
                x->max_ = std::max(x->interval_.high(), x->left_->max_);
            else if (x->right_)
                x->max_ = std::max(x->interval_.high(), x->right_->max_);
            else
                x->max_ = x->interval_.high();

            if (y->right_)
                y->max_ = std::max(y->interval_.high(), std::max(y->right_->max_, x->max_));
            else
                y->max_ = std::max(y->interval_.high(), x->max_);
        }

        void right_rotate(node_type* y)
        {
            auto* x = y->left_;
            y->left_ = x->right_;

            if (x->right_)
                x->right_->parent_ = y;

            x->parent_= y->parent_;
            if (!y->parent_)
                root_ = x;
            else if (y->is_left())
                y->parent_->left_ = x;
            else
                y->parent_->right_ = x;

            x->right_ = y;
            y->parent_ = x;

            // max fixup
            if (y->left_ && y->right_)
                y->max_ = std::max(y->interval_.high(), std::max(y->left_->max_, y->right_->max_));
            else if (y->left_)
                y->max_ = std::max(y->interval_.high(), y->left_->max_);
            else if (y->right_)
                y->max_ = std::max(y->interval_.high(), y->right_->max_);
            else
                y->max_ = y->interval_.high();

            if (x->left_)
                x->max_ = std::max(x->interval_.high(), std::max(x->left_->max_, y->max_));
            else
                x->max_ = std::max(x->interval_.high(), y->max_);
        }

        void recalculate_max(node_type* reacalculation_root)
        {
            auto* p = reacalculation_root;
            while (p && p->max_ <= reacalculation_root->max_)
            {
                if (p->left_ && p->left_->max_ > p->max_)
                    p->max_ = p->left_->max_;
                if (p->right_ && p->right_->max_ > p->max_)
                    p->max_ = p->right_->max_;
                p = p->parent_;
            }
        }

        void insert_fixup(node_type* z)
        {
            while (z->parent_ && z->parent_->color_ == rb_color::red)
            {
                if (!z->parent_->parent_)
                    break;
                if (z->parent_ == z->parent_->parent_->left_)
                {
                    node_type* y = z->parent_->parent_->right_;
                    if (y && y->color_ == rb_color::red)
                    {
                        z->parent_->color_ = rb_color::black;
                        y->color_ = rb_color::black;
                        z->parent_->parent_->color_ = rb_color::red;
                        z = z->parent_->parent_;
                    }
                    else
                    {
                        if (z == z->parent_->right_)
                        {
                            z = z->parent_;
                            left_rotate(z);
                        }
                        z->parent_->color_ = rb_color::black;
                        z->parent_->parent_->color_ = rb_color::red;
                        right_rotate(z->parent_->parent_);
                    }
                }
                else
                {
                    node_type* y = z->parent_->parent_->left_;
                    if (y && y->color_ == rb_color::red)
                    {
                        z->parent_->color_ = rb_color::black;
                        y->color_ = rb_color::black;
                        z->parent_->parent_->color_ = rb_color::red;
                        z = z->parent_->parent_;
                    }
                    else
                    {
                        if (z->is_left())
                        {
                            z = z->parent_;
                            right_rotate(z);
                        }
                        z->parent_->color_ = rb_color::black;
                        z->parent_->parent_->color_ = rb_color::red;
                        left_rotate(z->parent_->parent_);
                    }
                }
            }
            root_->color_ = rb_color::black;
        }

        void erase_fixup(node_type* x, node_type* x_parent, bool y_is_left)
        {
            while (x != root_ && x->color_ == rb_color::black)
            {
                node_type* w;
                if (y_is_left)
                {
                    w = x_parent->right_;
                    if (w->color_ == rb_color::red)
                    {
                        w->color_ = rb_color::black;
                        x_parent->color_ = rb_color::red;
                        left_rotate(x_parent);
                        w = x_parent->right_;
                    }

                    if (w->left_->color_ == rb_color::black && w->right_->color_ == rb_color::black)
                    {
                        w->color_ = rb_color::red;
                        x = x_parent;
                        x_parent = x->parent_;
                        y_is_left = (x == x_parent->left_);
                    }
                    else
                    {
                        if (w->right_->color_ == rb_color::black)
                        {
                            w->left_->color_ = rb_color::black;
                            w->color_ = rb_color::red;
                            right_rotate(w);
                            w = x_parent->right_;
                        }

                        w->color_ = x_parent->color_;
                        x_parent->color_ = rb_color::black;
                        if (w->right_)
                            w->right_->color_ = rb_color::black;

                        left_rotate(x_parent);
                        x = root_;
                        x_parent = nullptr;
                    }
                }
                else
                {
                    w = x_parent->left_;
                    if (w->color_ == rb_color::red)
                    {
                        w->color_ = rb_color::black;
                        x_parent->color_ = rb_color::red;
                        right_rotate(x_parent);
                        w = x_parent->left_;
                    }

                    if (w->right_->color_ == rb_color::black && w->left_->color_ == rb_color::black)
                    {
                        w->color_ = rb_color::red;
                        x = x_parent;
                        x_parent = x->parent_;
                        y_is_left = (x == x_parent->left_);
                    }
                    else
                    {
                        if (w->left_->color_ == rb_color::black)
                        {
                            w->right_->color_ = rb_color::black;
                            w->color_ = rb_color::red;
                            left_rotate(w);
                            w = x_parent->left_;
                        }

                        w->color_ = x_parent->color_;
                        x_parent->color_ = rb_color::black;
                        if (w->left_)
                            w->left_->color_ = rb_color::black;

                        right_rotate(x_parent);
                        x = root_;
                        x_parent = nullptr;
                    }
                }
            }

            x->color_ = rb_color::black;
        }

    private:
        node_type* root_;
        size_type size_;
    };
//############################################################################################################
	template <typename T, typename Kind = closed>
	using interval_tree_t = interval_tree <interval <T, Kind>>;
//############################################################################################################
}
