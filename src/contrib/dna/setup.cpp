namespace Dna {

	const map<Sensitivity, SensitivityTraits> sensitivity_traits = {
		{{Sensitivity::FASTER, {true,  false, 20.0,  9,    0,      0,      0,     1,        16,   nullptr,  1.0,     2.0,       dna,      12, 0.0, 40, 0.1 }},
		{Sensitivity::FAST, {true,  false, 20.0,  9,    0,      0,      0,     1,        16,   nullptr,  1.0,     2.0,       dna,      10, 0.0, 40, 0.5}},
		{Sensitivity::DEFAULT, {true,  false, 20.0,  9,    0,      0,      0,     1,        16,   nullptr,  1.0,     2.0,       dna,      10, 0.0, 20, 0.5 }},
		{Sensitivity::SENSITIVE, {true,  false, 20.0,  9,    0,      0,      0,     1,        16,   nullptr,  1.0,     2.0,       dna,      6, 0.0, 20, 0.5 }},
		{Sensitivity::VERY_SENSITIVE, {true,  false, 20.0,  9,    0,      0,      0,     1,        16,   nullptr,  1.0,     2.0,       dna,      5, 0.0, 17, 0.5 }},
		{Sensitivity::ULTRA_SENSITIVE, {true,  false, 20.0,  9,    0,      0,      0,     1,        16,   nullptr,  1.0,     2.0,       dna,      4, 0.0, 15, 0.5 }},
		} };

    const map<Sensitivity, vector<string>> shape_codes = {
        { Sensitivity::ULTRA_SENSITIVE,
      {"111111111111"} },
      { Sensitivity::VERY_SENSITIVE,
       { "1111111111111"} },
    { Sensitivity::SENSITIVE,
     { "11111111111111" } },
    { Sensitivity::DEFAULT,
      { "111111111111111"} },
    { Sensitivity::FAST,
      { "111111111111111" } },
   { Sensitivity::FASTER,
      { "111111111111111111" } } };

    void setup_search(Sensitivity sens, Search::Config& cfg)
{
	const SensitivityTraits& traits = sensitivity_traits[(int)align_mode.sequence_type].at(sens);
	config.sensitivity = sens;
	::Config::set_option(cfg.chain_fraction_align, config.chain_fraction_align_, 0.0, traits.chain_fraction_align);
    ::Config::set_option(cfg.min_chain_score, config.min_chain_score_, 0, traits.min_chain_score);
    ::Config::set_option(cfg.max_overlap_extension, config.max_overlap_extension_, 0.0, traits.max_overlap_extension);

    if(align_mode.mode == AlignMode::blastn)
       Reduction::reduction = dna;
    cfg.chain_pen_gap = config.chain_pen_gap_scale * 0.01 * (double)shapes[0].length_;
    cfg.chain_pen_skip = config.chain_pen_skip_scale * 0.01 * (double)shapes[0].length_;

}


}