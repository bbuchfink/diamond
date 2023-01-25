#pragma once
#include "../dp/flags.h"

namespace Extension {

enum class Mode {
	BANDED_FAST, BANDED_SLOW, FULL, GLOBAL
};

int band(int len, const Mode mode);
HspValues filter_hspvalues();

}

template<>
struct EnumTraits<Extension::Mode> {
	static const SEMap<Extension::Mode> from_string;
};
