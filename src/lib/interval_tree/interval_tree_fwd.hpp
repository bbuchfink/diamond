#pragma once

namespace lib_interval_tree
{
    template <typename numerical_type, typename interval_kind_>
    struct interval;

    template <typename IntervalT>
    class interval_tree;

    template <typename numerical_type, typename interval_type>
    class node;

    template <typename node_type, typename owner_type>
    class basic_interval_tree_iterator;

    template <typename node_type>
    class const_interval_tree_iterator;

    template <typename node_type>
    class interval_tree_iterator;
}
