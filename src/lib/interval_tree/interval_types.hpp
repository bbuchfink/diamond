#pragma once

namespace lib_interval_tree
{
    // (]
    struct left_open
    {
        template <typename numerical_type>
        static inline bool within(numerical_type b, numerical_type e, numerical_type p)
        {
            return (b < p) && (p <= e);
        }
    };
    // [)
    struct right_open
    {
        template <typename numerical_type>
        static inline bool within(numerical_type b, numerical_type e, numerical_type p)
        {
            return (b <= p) && (p < e);
        }
    };
    // []
    struct closed
    {
        template <typename numerical_type>
        static inline bool within(numerical_type b, numerical_type e, numerical_type p)
        {
            return (b <= p) && (p <= e);
        }
    };
    // ()
    struct open
    {
        template <typename numerical_type>
        static inline bool within(numerical_type b, numerical_type e, numerical_type p)
        {
            return (b < p) && (p < e);
        }
    };
}
